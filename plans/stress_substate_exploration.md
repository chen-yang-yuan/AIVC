# Tumor stress vs target-gene substates (`code/1_exploration/`)

Current as of 2026-11-14. Design (2026-11-07), first run (2026-11-08), the revised run after the user's review
(2026-11-08: two gene-list variants and two correspondence approaches as peers, outputs split into
`stress_axes/`, `cytoplasmic/`, `nuclear/`, `comparison/`) and the streamlined run (2026-11-08: canonical gene
lists only, two rule-based stress labellings, no k-means / UMAP), followed by the Sankey-style flow figures
(2026-11-08) and an audit of code, outputs and this document (2026-11-08). This document holds the latest state of
the analysis design and its results. Rewrite the affected sections when something changes
and append a dated entry to the Log (`CLAUDE.md` maintenance rule).

## 1. Purpose and scope

The sample-corrected pooled subtyping (`plans/subtype_exploration.md` §5) found tumor substates in the
cytoplasmic expression of the 634 target genes (S-phase, G2/M, Interferon, Hypoxia, ER stress, Stress response,
Complement, Stroma border, Myeloid border, Baseline) and a different, weaker organisation in the nuclear expression
of the same genes (ARI 0.05 with the cytoplasmic labels). The model will predict cytoplasmic activity from the
cell's identity and environment. Two questions follow:

1. Do the target-gene substates correspond to tumor stress, i.e. do cells in a substate sit in a recognisable
   stress condition (hypoxic, crowded, under immune attack, at a stromal border, heat-shocked, granule-forming)?
2. Does the cytoplasmic expression of the target genes carry more information about stress than the nuclear
   expression of the same genes? If so, the choice of cytoplasmic activity as the model's outcome is supported
   by the environment-to-outcome link, not only by the substate structure.

Scope: the five subtyped samples (BC, OC, CC, LC, Prostate; Skin excluded as in the subtyping), the 623,838
pooled tumor cells (`cell_type_merged == "Malignant cell"`), the 634-gene set unchanged. The pooled cell is the
unit of analysis; sample is a covariate and per-sample breakdowns are robustness checks. The design has to make
a cytoplasmic advantage defensible against three known biases (depth, circularity, spatial autocorrelation) rather
than assume it; either outcome of question 2 is recorded.

## 2. The previous stress analysis and what is kept

`code/2_outcome/5_score.ipynb` scored the six-sample merged tumor object (`data/merged_data/adata_tumor_scored.h5ad`,
672,964 cells, whole-cell `normalize_total(1e4)` + `log1p`):

| Old score | Definition | Group rule | Kept? |
|---|---|---|---|
| `hypoxia_score` | 2nd-highest log expression among HIF1A, EPAS1, NFE2L2, CREB1, RELA, RELB, NFKB1, NFKB2 | `<= median -> Low` (522k Low / 151k High: ties) | column kept for continuity; no longer the hypoxia measure |
| `heat_shock_score` | same rule over 11 HSP genes (HSP90AA1, HSP90AB1, HSPA1A, HSPA1B, HSPA6, HSPA8, HSPH1, DNAJB1, HSPB1, HSPD1, HSPE1) | same | genes kept, scored with `score_genes` |
| `mechanical` | inverse mean distance to the 10 nearest **tumor** cells, z-scored | median | kept, over all cells |
| `immune_cell_proximity`, `tcell_proximity` | count of immune / T cells within 50 um | Away (0) / Close | kept |
| `tcell_attack_weighted` | T cells within 100 um, Gaussian kernel (sigma 50), weighted by their GZMA/GZMB/GZMK, z-scored | median | kept |

Refinements in this analysis:

- Five samples, pooled thresholds, every score and group written to a CSV (`output/1_exploration/stress/
  stress_axes/stress_scores.csv`), not only into an h5ad. The old notebook cannot be re-run as it stands (the BC ssGSEA file
  changed shape), so the scores are recomputed, not read from the old object.
- Crowding over all cells of the tissue (the old code fitted the neighbour search on tumor cells only, despite
  its comment).
- The 8-TF hypoxia rule is retired as the hypoxia measure: it is near zero for most cells at Xenium depth, and
  three of the eight TFs are in the 634 set. Hypoxia gets an expression-free proxy (vascular distance) and a
  canonical expression signature (HALLMARK_HYPOXIA).
- Away/Close is kept for counts, and every continuous score is split by pooled tertiles (top third vs bottom
  third, middle excluded for binary contrasts) rather than `<= median`, which sent tied cells to Low.
- New axes: core vs boundary, myeloid contact, stromal contact, morphology, SG marker expression.

## 3. Stress axes

All axes are computed per sample (memory) and then concatenated; thresholds are set on the pooled cells. Spatial
scores are in physical units (um, counts) and are directly comparable across samples. Expression scores use the
whole-cell log-normalised counts (all 5,001 genes), which is independent of the nuclear/cytoplasmic split used to
define the substates.

**Primary tier: environmental, expression-free.** This tier is the fair test of question 2, since it uses no
transcript of the tumor cell itself.

| Axis | Score(s) | Group rule (pooled) |
|---|---|---|
| Vascular distance (hypoxia proxy) | distance to the nearest Endothelial cell; endothelial count within 50 and 100 um | tertiles: Far vs Near |
| Crowding (mechanical) | inverse mean distance to the 10 nearest cells of any type; `n_cells_30um` | tertiles: High vs Low |
| Core vs boundary | tumor fraction among cells within 50 um; distance to the nearest non-tumor cell | tertiles: Core vs Boundary |
| Immune attack | T-cell (any: CD4+ T, CD8+ T, T cell) count within 50 um and distance; GZM-weighted Gaussian kernel (100 um, sigma 50, weights = mean GZMA/GZMB/GZMK of the neighbouring T cells) and its label-free version over all neighbouring non-tumor cells (§3.1); immune (any) count within 50 um | Away/Close on the counts; tertiles on the distance and the kernel score |
| Myeloid contact | Myeloid + Dendritic count within 50 um; distance | Away/Close; tertiles |
| Stromal contact | Fibroblast (CAF) count within 50 um; distance | Away/Close; tertiles |
| Morphology (mechanical, secondary) | cell area; nucleus-to-cell area ratio | tertiles |

Cells whose sample has no cell of a reference type get the count 0 and the distance capped at the pooled 99th
percentile; the cap is recorded.

**Secondary tier: expression-based stress** (`sc.tl.score_genes` on whole-cell log-normalised counts, z-scored
on the pooled cells, tertile split), with the canonical gene lists as they stand on the panel:

- hypoxia: HALLMARK_HYPOXIA (`data/_utils/hallmark_pathways_filtered.gmt`), plus the old 8-TF score as a column;
- heat shock: the 11 HSP genes above;
- TNFA/NF-kB and interferon-gamma response hallmarks (inflammatory / immune stress);
- **stress-granule (SG) marker expression**: the SG transcriptome list of `data/_utils/SG_markers.xlsx` with the
  rule of `code/2_outcome/3_detection.py` (`Fraction of RNA molecules in SGs > 0.4`: 406 genes, 129 on the
  panel), scored on whole-cell expression; where the SG detection was run, the granule burden per cell
  (granule count and SG-marker transcripts in granules, from `processed_data/granule_adata_0.25_1.h5ad`, log1p)
  is added as a companion score.

**Gene lists.** Every expression signature is scored with its canonical gene list as it stands on the panel.
Overlaps with the 634 set: HALLMARK_HYPOXIA 26 of 96 (SLC2A1, LDHA, NDRG1, BNIP3, P4HA1, ALDOA among others),
heat shock 1 of 11 (HSPA8), TNFA/NFKB 17 of 110, IFNG response 24 of 106, SG markers 21 of 129; the 8-TF rule
shares HIF1A, EPAS1 and NFE2L2, GZMA/B/K none. Whole-cell counts are nuclear plus cytoplasmic, and for a
cytoplasm-enriched gene the cytoplasmic part dominates, so a signature that shares genes with the 634 set is
mechanically closer to the cytoplasmic matrix than to the nuclear one (the circularity caveat of §6).

*Why the minus-634 arm was retired (2026-11-08).* The review run of 2026-11-08 scored every signature twice,
canonical and minus the 634 target genes (`no634`), and carried both through every table, stress labelling and
figure as peers, so that the user could pick the arm with the stronger signal. The canonical lists won on every
metric and the minus-634 lists showed the same associations at lower strength: approach 1 on the cytoplasmic
substates max AMI 0.059 vs 0.026, max Cramer's V 0.398 vs 0.270, max R2 beyond sample 0.115 vs 0.052,
Interferon x IFNG log2 OR +3.30 vs +2.01, Hypoxia x hypoxia +2.78 vs +1.28; the expression archetypes AMI 0.045
vs 0.021; for the cytoplasmic-vs-nuclear comparison the canonical rows kept a small cytoplasmic lead (AMI ratio
1.20-1.23) where the minus-634 rows were level (1.04-1.13). Removing a quarter of a signature's genes (hypoxia,
TNFA/NFKB, IFNG) lowered its signal-to-noise at Xenium depth without changing what it detected, so the second
arm added cost (double axis tables, labellings and figures) but no decision. The arm was therefore retired on
the user's instruction and the analysis runs on the canonical lists only; its purpose, guarding against the
shared-gene circularity, is kept as an explicit caveat (§6) and is covered by the expression-free axes and the
environmental archetypes, which use no transcript of the tumor cell. The retired numbers are recorded in the
Log entry of the review run.

### 3.1 Cell-type annotation across samples

The five samples share the 18-category `cell_type_merged` vocabulary but not the annotation route
(`code/1_preprocessing/1_clean.ipynb`): BC, OC and CC carry 10x fine labels (31-33 types such as Macrophage, cDC2,
CD8+ Tex, iCAF_1, myCAF, VenECs) mapped to the merged vocabulary by `data/_utils/cell_type_dict.pkl`; LC and
Prostate (and Skin) were annotated by hand at the level of 10x graph clusters (33 and 22 clusters,
`data/_utils/cell_type_dict_manual.pkl`), unmapped clusters became Unknown and ambiguous ones Mixed. Share of all
cells (%), checked 2026-11-07:

| | BC | OC | CC | LC | Prostate |
|---|---|---|---|---|---|
| B cell | 10.0 | 0.3 | 2.1 | 6.4 | 0 |
| CD4+ T cell | 11.8 | 0.6 | 4.4 | 8.4 | 0 |
| CD8+ T cell | 7.4 | 0.8 | 2.5 | 5.8 | 3.4 |
| Dendritic cell | 5.7 | 0.5 | 2.7 | 7.7 | 0 |
| Myeloid cell | 8.2 | 4.8 | 3.9 | 8.1 | 3.9 |
| Endothelial cell | 4.8 | 5.0 | 4.3 | 13.7 | 5.1 |
| Fibroblast (CAF) | 13.0 | 17.9 | 29.4 | 12.9 | 10.9 |
| Epithelial (non-malignant) | 8.5 | 2.5 | 1.9 | 7.7 | 0 |
| Mixed / Unknown | 0 / 0 | 0 / 0 | 0 / 0 | 2.1 / 0.7 | 0 / 11.0 |

Prostate has no B-cell, CD4+ T, dendritic or non-malignant epithelial label and 11 % Unknown: its single T
cluster is labelled CD8+ T cell and rare types were absorbed into Unknown or neighbouring clusters. The cluster
route cannot resolve rare types, so fine immune counts are not comparable between the two routes. Endothelial,
fibroblast and myeloid labels exist in every sample at plausible shares, so the vascular, stromal and myeloid
axes are the least affected; crowding, morphology and core/boundary do not use fine labels. The substate labels
already inherit these annotations (`spatial_context` distances, Stroma border and Myeloid border overrides).

Rules that follow:

1. **Coarse reference sets only** for label-based axes: T cell (any) = CD4+ T + CD8+ T + T cell; immune (any) =
   B + T (any) + Dendritic + Myeloid + Mast; myeloid = Myeloid + Dendritic; endothelial; fibroblast. No
   CD8-specific or B-cell axis anywhere (a CD8 distance for BC, OC, CC was considered in the design and not
   implemented, since Prostate's single T cluster makes it incomparable).
2. **Away/Close is the primary group rule for label-based counts** (0 vs >= 1 within 50 um): it depends on the
   presence of a labelled cell, not on how many were resolved. Tertiles stay for the continuous scores.
3. **Unknown and Mixed cells** are excluded from every reference set and from the denominator of the tumor
   fraction (core/boundary uses labelled non-tumor cells only). Prostate's 11 % Unknown is reported next to its
   results, and Prostate is the first sample dropped in leave-one-sample-out.
4. **Label-free companion for every contact axis** (the harmonisation check): a Gaussian-kernel score (100 um,
   sigma 50) of lineage-marker expression in the neighbouring non-tumor cells, whatever their label. Panel
   markers: vascular PECAM1, CDH5, CLDN5, KDR, ENG; cytotoxic CD3E, CD2, CD8A, GZMA, GZMB, GZMK, PRF1; myeloid
   CD68, CD163, CD14, ITGAM; stromal only PDGFRB and FAP (no collagens, DCN or LUM on the panel; flagged as weak).
   A label-based association counts as robust only when its label-free companion shows the same sign, and the
   information-content comparison (analysis C) is reported on both versions. The GZM-weighted T-cell-attack score
   therefore exists in two forms: restricted to labelled T cells and over all neighbouring non-tumor cells.
5. **Optional follow-up** (not in the first pass): transfer the fine labels of BC, OC and CC to LC and Prostate on
   the shared panel (ingest or a kNN classifier on log-normalised expression) so that all five samples sit on one
   route (§9).

## 4. Substate labellings compared

Both labellings exist for the same pooled cells:

| Labelling | Source | Granularity |
|---|---|---|
| Cytoplasmic | `output/1_exploration/pooled_sample_corrected/annotation/cell_subtypes_pooled.csv` | 16 Louvain clusters, 11 program labels |
| Nuclear | `output/1_exploration/pooled_sample_corrected_nuclear/annotation/cell_subtypes_pooled.csv` | 21 Louvain clusters, 12 program labels |

Per-cell program scores (23 `score_*` columns per compartment) come from `characteristics/cell_annotations_
{cytoplasmic,nuclear}.csv`; the 50 sample-corrected PCs and the 30-NN expression graphs from the two embedded
objects (`adata_pooled_cytoplasmic_embedded.h5ad`, `adata_pooled_nuclear_embedded.h5ad`).

Labels that are spatial by construction are flagged: Stroma border (fibroblast within 25 um) is excluded from the
stromal-contact axis and Myeloid border from the myeloid-contact axis. Low depth cells are excluded everywhere.
Analyses run at both granularities; the label level is the one reported.

### 4.1 Explicit stress subtypes (approach 2)

Besides the per-axis groups (approach 1: tertiles / Away-Close of one axis at a time), every pooled cell gets
two explicit stress labels so that the correspondence can be read label against label. Composite scores are
the mean of the sign-aligned pooled z of their member axes, re-standardised (`stress_axes/composites.csv`):

| Composite | Members (sign) | Tier |
|---|---|---|
| Vascular deprivation | vascular_dist (+), vascular_n100 (-), k_vascular (-) | env |
| Crowding | crowding (+), n_cells_30um (+) | env |
| Core isolation | tumor_frac (+), nontumor_dist (+) | env |
| Immune attack | tcell_n50, tcell_gzm_kernel, immune_n50, k_cytotoxic (+), tcell_dist (-) | env |
| Myeloid contact | myeloid_n50, k_myeloid (+), myeloid_dist (-) | env |
| Stromal contact | fibroblast_n50, k_stromal (+), fibroblast_dist (-) | env |
| Hypoxia (expr), Heat shock, Inflammatory (TNFA/NFKB), Interferon (IFNG) | the canonical signature score | expr |

Two labellings (`stress_axes/stress_subtypes.csv`): `stress_env` = the environmental composite with the largest
score if that score exceeds 0.5 (standardised units), else "Unstressed" (expression-free, the fair labelling for
question 2); `stress_expr` = the same rule over the four expression composites. The SG-marker score and the two
granule-burden axes remain stress axes of approach 1 (profiles, boxplots, dot plot) but are not subtyping
features (user decision 2026-11-08: granule-level scores are to be read across the substates, not used to define
stress classes). Sizes per sample, the environment x response cross-tab (`stress_env_vs_expr.csv/png`) and
spatial maps of both labellings are in `stress_axes/`.

**Why two separate labellings rather than one pooled labelling.** The review run also clustered the cells on
the pooled environmental + expression composites (MiniBatchKMeans, k by silhouette); that labelling was dropped
on 2026-11-08 and the two rule-based labellings are kept apart, for these reasons:

1. *The two feature sets measure different things.* The environmental composites are exposure, measured from
   the neighbourhood without any transcript of the cell; the expression composites are the cell's own
   transcriptional response. A pooled cluster such as "Core isolation + Vascular deprivation" or "Inflammatory
   (TNFA/NFKB) + Hypoxia (expr)" (the names the k-means centres received) cannot say whether a substate tracks
   the environment or the response, and that distinction is exactly what question 1 asks and what the virtual
   cell model is built to predict (environment -> response).
2. *Question 2 needs an expression-free labelling.* The expression composites are whole-cell scores of the same
   cell and share genes with the 634 set (§3), so any labelling that includes them inherits the shared-gene bias
   towards the cytoplasmic matrix. Only a labelling built from the environmental composites alone makes the
   cytoplasmic-vs-nuclear comparison fair; pooling would remove that control.
3. *Pooling drowns the expression information.* The environmental composites are smooth over hundreds of
   micrometres and dominate the variance of the pooled space; the k-means solution had a silhouette of 0.12 at
   the chosen k = 5 (clusters are arbitrary cuts of continuous gradients) and four of its five clusters were
   named by environmental composites, with the expression programmes appearing only as secondary members. The
   separate expression archetype recovers IFNG and hypoxia responses that the pooled clustering blurred.
4. *Separate labels make the environment x response table an output.* With `stress_env` and `stress_expr` kept
   apart, their cross-tab (which environments co-occur with which responses, pooled and per sample) is a
   direct, interpretable result; a pooled label hides it inside the cluster definitions. The dominant-archetype
   rule is also transparent and reproducible across samples, whereas a k-means partition changes with k, the
   subsample and the scaling.
5. *Limitation.* The rule assigns one stress per cell, so cells under several stresses are under-represented.
   The composite scores are saved per cell in `stress_subtypes.csv` for any multi-label reading.

## 5. Analyses and results (streamlined run 2026-11-08, `6_stress_substates.ipynb`, `output/1_exploration/stress/`)

Everything runs on the pooled cells with spatial block bootstrap (200 um tiles within sample, B = 100) for every
CI and block folds for every cross-validation. 28 axes: 20 expression-free (`env`), 5 canonical signature scores
(`signature`), 3 other (8-TF rule, granule count, SG transcripts in granules). Low depth cells are dropped from
the label analyses (cytoplasmic n = 606,329; nuclear n = 570,991); Stroma border and Myeloid border are dropped
from the axes that define them. Two correspondence approaches run side by side:

- **Approach 1, stress groups across substates**: A median z per substate and axis, R2 of the axis explained by
  the substate beyond sample (bootstrap CI), epsilon2; B pooled tertile (Away/Close) groups vs substates,
  Cramer's V, AMI, log2 odds ratios pooled / Mantel-Haenszel / per sample, reproducible = |log2 OR| > 0.5, CI
  excludes 0, same sign in >= 3 of 5 samples (a table column; no longer drawn on the figures).
- **Approach 2, explicit stress subtypes**: the two labellings of §4.1 vs the substates, contingency tables,
  log2 enrichment (observed / expected) per cell with the same reproducibility rule, Cramer's V and AMI with CIs.

### 5.1 Stress subtypes (`stress_axes/`)

The environmental archetypes split the 623,838 cells into seven groups (Core isolation 107k, Stromal contact
106k, Vascular deprivation 93k, Crowding 92k, Myeloid contact 81k, Unstressed 80k, Immune attack 64k); the
expression archetypes into Unstressed 193k, Heat shock 150k, Hypoxia 101k, Interferon 91k, Inflammatory 89k
(`stress_subtype_sizes.csv`, also per sample: CC supplies most of the Core isolation / Crowding / Hypoxia cells,
LC and BC most of the Immune attack cells, OC most of the Heat shock cells). Without the Stress granule composite
the 112k cells it used to claim went mostly to Unstressed (+52k) and Heat shock (+24k). The environment x
response cross-tab (`stress_env_vs_expr.csv/png`; log2 observed / expected, pooled, with the per-sample sign
count) is the direct reading of which environment elicits which response: Core isolation -> Hypoxia (expr)
+0.79 (28 % of its cells; same sign in 4/5 samples), Vascular deprivation -> Hypoxia +0.44 (22 %; 2/5, OC and
CC drive it), Immune attack -> Interferon (IFNG) +0.86 (26 %; 4/5), Myeloid contact -> Interferon +0.40 and
Inflammatory +0.32 (5/5), Stromal contact -> Heat shock +0.30 and away from Hypoxia / Interferon (-0.58 each,
4/5 and 5/5), environmentally Unstressed -> expression Unstressed +0.36 and away from Hypoxia (-1.20, 4/5).
Crowding has no expression partner (Unstressed +0.20). Every environment still holds 25-40 % expression-
Unstressed cells, so the environment shifts the response distribution rather than fixing it. The same table as
a flow figure: `sankey_env_vs_expr.png` (ribbon width = cells, ribbons with log2 enrichment > 0.5 opaque, the
rest faint). Per-sample
composition (`stress_subtype_composition_by_sample.png`) and the maps (`spatial_stress_subtypes_<sample>.png`)
show that CC is the crowded, core-isolated sample and LC the vessel-rich, immune-rich one, as in the per-sample
axis medians (`stress_axes_summary.csv`).

### 5.2 Cytoplasmic substates vs stress (`cytoplasmic/`, question 1)

**Approach 1.** Median / max over the axis sets (Cramer's V; AMI; R2 beyond sample):

| Axis set | median V / max V | median AMI / max AMI | median R2 / max R2 | reproducible pairs |
|---|---|---|---|---|
| env (20) | 0.233 / 0.338 | 0.020 / 0.042 | 0.035 / 0.082 | 80 |
| signature (5) | 0.299 / 0.398 | 0.034 / 0.059 | 0.072 / 0.115 | 22 |
| other (3) | 0.129 / 0.131 | 0.006 / 0.006 | 0.006 / 0.008 | 0 |

Strongest reproducible pairs (pooled log2 OR, Mantel-Haenszel in brackets): Interferon x IFNG +3.30 [+3.46],
Hypoxia x hypoxia +2.78 [+3.07], Stroma border x tumor core -2.43 [-3.18] and x non-tumor distance -2.29,
Stress response x TNFA/NFKB +2.05 [+1.51], x non-tumor distance +2.01, x vascular kernel -1.96, x vascular
distance +1.95 (5/5), x fibroblast distance +1.93, x tumor core +1.86, x cytotoxic kernel -1.77, x crowding
-1.57 [-2.50] (5/5); Baseline x hypoxia -1.73 (5/5); Complement x TNFA/NFKB +1.71 [+2.02] (5/5). The granule
axes (SG-marker score, granule count, SG transcripts in granules) show no reproducible association with any
substate (max V 0.13, R2 <= 0.008; `boxplot_axes_by_substate.png` panels 13-15: the pooled medians are flat
across substates and the per-sample medians sit on top of each other). Per substate, with the expression-free
axes only:

| Substate | Reproducible environmental associations (log2 OR) |
|---|---|
| Hypoxia | vessel-poor -1.4, vascular kernel -1.3, core +1.3, far from vessels +1.3, far from fibroblasts +1.2, far from non-tumor cells +1.2, cytotoxic kernel -1.1, far from T cells +1.1 |
| Interferon | cytotoxic kernel +1.2, close to T cells (distance -1.1, count +0.8), close to myeloid cells (-1.1, kernel +0.9, count +0.8), immune count +0.8 |
| Stress response | far from non-tumor cells +2.0, vascular kernel -2.0, far from vessels +1.9, far from fibroblasts +1.9, core +1.9, cytotoxic kernel -1.8, stromal kernel -1.7, crowding -1.6, vessel-poor -1.5, myeloid kernel -1.5, sparse (n_cells_30um -1.3), large cells +1.1 |
| Complement / acute phase | vascular kernel +1.5, close to T cells -1.4, myeloid kernel +1.3, cytotoxic kernel +1.3, crowding +1.2, small cells -1.2, immune / T / myeloid counts +1.0 to +1.1 |
| Baseline | close to fibroblasts (-1.1, count +0.9), close to vessels (-0.8, count +0.8) |
| S-phase, G2/M | nucleus ratio +0.8; core +0.7, stromal kernel +0.6 (weak) |
| ER stress | sparse (n_cells_30um -0.9) only |
| Stroma border, Myeloid border | boundary -2.4 / -2.3, immune and myeloid contact + (positive controls) |

**Approach 2** (Cramer's V [CI]; AMI [CI]; reproducible cells of the table; strongest enrichment):

| Stress labelling | V | AMI | repro. | Enriched pairs (substate -> stress subtype, log2 enrichment, share of the substate's cells) |
|---|---|---|---|---|
| stress_env (7) | 0.152 [0.149, 0.156] | 0.032 [0.030, 0.033] | 17 | Stroma border -> Stromal contact +1.27 (41 %); Stress response -> Core isolation +1.24 (41 %), -> Vascular deprivation +0.52 (22 %); Interferon -> Immune attack +0.71 (17 %); Hypoxia -> Vascular deprivation +0.63 (23 %), -> Core isolation +0.50 (24 %); Complement -> Immune attack +0.58; ER stress -> Unstressed +0.62 |
| stress_expr (5) | 0.233 [0.227, 0.240] | 0.048 [0.046, 0.051] | 20 | Interferon -> Interferon (IFNG) +1.85 (53 %); Hypoxia -> Hypoxia (expr) +1.33 (41 %); Complement -> Inflammatory +1.27 (34 %); Stroma border -> Inflammatory +1.03 (29 %); Stress response -> Inflammatory +0.91 (27 %) |

The explicit stress subtypes tell the same story as the per-axis groups with one addition: the share of a
substate's cells that actually sits in the matching stress subtype is modest (Interferon 53 % in the IFNG
archetype, Hypoxia 41 % in the hypoxia archetype and 23 % in Vascular deprivation; Stress response 41 % in Core
isolation), so the substates are stress-enriched populations, not stress classes.

**Signal strength for question 1** (`cytoplasmic/signal_strength.csv`): the strongest signal comes from
approach 1 with the signature axes (max AMI 0.059, V 0.398, R2 0.115, Interferon x IFNG +3.30), followed by
approach 2 with the expression archetype (AMI 0.048, V 0.233); the expression-free side gives approach 1 env
axes max AMI 0.042 / V 0.338 (core / boundary and stromal axes) and the environmental archetype AMI 0.032 /
V 0.152.

Intuitive figures: `boxplot_axes_by_substate.png` (raw units; pooled medians: distance to the nearest endothelial
cell 98 um for Stress response and 92 um for Hypoxia vs 49-70 um for the other substates; cells within 30 um 20
for Stress response vs 25-26 for Hypoxia and the proliferating substates; hallmark hypoxia score 0.20 for Hypoxia
vs 0.05-0.14 elsewhere; IFNG score 0.10 for Interferon vs <= 0.02 elsewhere; granule-level scores flat),
`stacked_groups_by_substate.png`, `dotplot_substate_axes.png`, `stress_subtype_enrichment.png` (log2 enrichment
with cell counts), `stress_subtype_composition.png`, `sankey_substates_vs_stress.png` (three-column flow:
environmental archetype | substates | expression archetype; ribbon width = cells, coloured by substate, the
reproducibly enriched pairs opaque and everything else faint, so the Interferon -> Immune attack / IFNG,
Hypoxia -> Vascular deprivation / Hypoxia, Stress response -> Core isolation / Inflammatory and Stroma border ->
Stromal contact flows stand out against the bulk of unspecific flow), `spatial_Hypoxia_vs_vascular_dist.png`,
`spatial_Interferon_vs_k_cytotoxic.png`, `spatial_Stress_response_vs_crowding.png` (the Stress response
substate fills the interior of the CC tumor nests where crowding is lowest, and the same niche in BC and OC),
`spatial_labels_vs_stress_<sample>.png` (substates | environmental archetype | expression archetype),
`forest_headline.png` (per-sample log2 ORs: the Hypoxia-vessel association is OC- and CC-driven with LC
opposite; the Stress response niche is CC-driven with the same sign in BC, OC and LC but not Prostate;
Complement-immune is Prostate-driven).

### 5.3 Nuclear substates vs stress (`nuclear/`, same analysis)

| Axis set | median V / max V | median AMI / max AMI | median R2 / max R2 | reproducible pairs |
|---|---|---|---|---|
| env (20) | 0.175 / 0.256 | 0.011 / 0.025 | 0.019 / 0.048 | 83 |
| signature (5) | 0.266 / 0.327 | 0.027 / 0.042 | 0.069 / 0.089 | 26 |
| other (3) | 0.084 / 0.123 | 0.003 / 0.006 | 0.004 / 0.007 | 0 |

Strongest pairs: Hypoxia x hypoxia +2.73, Stress response x TNFA/NFKB +2.62, Interferon x IFNG +2.49, Stress
response x close to fibroblasts -2.07, x small cells -1.89, x boundary -1.84, x close to vessels -1.82,
Complement x IFNG +1.72, EMT/stroma x close to fibroblasts -1.42. The nuclear Hypoxia substate has the same
vascular / core profile as the cytoplasmic one (vessel-poor -1.4, far from vessels +1.3, core +1.3); the
nuclear Interferon substate is **not** immune-adjacent (far from fibroblasts +1.1, far from vessels +0.9, core
+0.8; no cytotoxic or T-cell association) although it expresses the IFNG programme; the nuclear Stress response
substate is a small (about 12k cells) boundary, vessel-near, immune-adjacent population with high TNFA/NFKB, the
opposite environment of the cytoplasmic Stress response substate (consistent with the 0.02 label agreement in
`plans/subtype_exploration.md`). Approach 2: stress_env V 0.106 [0.103, 0.109], AMI 0.017 [0.016, 0.018] (13
reproducible cells; strongest Hypoxia x Stromal contact -1.13; enriched EMT/stroma -> Stromal contact +0.83,
Interferon -> Core isolation +0.67, Hypoxia -> Vascular deprivation +0.65); stress_expr V 0.205 [0.200, 0.211],
AMI 0.038 [0.037, 0.041] (Stress response -> Inflammatory +1.74 (48 %), Interferon -> Interferon +1.38 (38 %),
Hypoxia -> Hypoxia +1.32 (42 %)). The nuclear labelling reaches about the same strength as the cytoplasmic one
on the signature axes and the expression archetype but clearly less on the expression-free axes and on the
environmental archetype (`nuclear/sankey_substates_vs_stress.png`: on the expression side the opaque flows match
the cytoplasmic ones, on the environmental side the nuclear Interferon substate flows to Core isolation rather
than Immune attack, Hypoxia to Core isolation / Vascular deprivation as in the cytoplasmic labelling, and the
nuclear Stress response and EMT/stroma substates to Stromal contact).

### 5.4 Cytoplasmic vs nuclear (`comparison/`, question 2)

**Correspondence statistics side by side** (`AB_comparison.csv`, `stress_subtype_comparison.csv`,
`signal_strength.csv`; medians over the axis set for approach 1, the labelling for approach 2):

| Approach / axis set | AMI cyto vs nuc (ratio) | V cyto vs nuc | R2 beyond sample cyto vs nuc (diff, CI) |
|---|---|---|---|
| 1 / env axes (20) | 0.020 vs 0.011 (1.76), cyto higher on 19 | 0.233 vs 0.175 | 0.035 vs 0.019 (+0.016 [+0.014, +0.021]; CI excludes 0 on 20/20) |
| 1 / signature axes (5) | 0.034 vs 0.027 (1.23), 4 of 5 | 0.299 vs 0.266 | 0.072 vs 0.069 (+0.003 [+0.001, +0.004]) |
| 1 / other axes (3) | 0.006 vs 0.003, 3 of 3 | 0.129 vs 0.084 | 0.006 vs 0.004 |
| 2 / stress_env | 0.032 vs 0.017 (1.91; diff +0.015 [+0.014, +0.016]) | 0.152 vs 0.106 (+0.046 [+0.044, +0.049]) | |
| 2 / stress_expr | 0.048 vs 0.038 (1.25; +0.010 [+0.009, +0.011]) | 0.233 vs 0.205 (+0.028 [+0.027, +0.031]) | |

The cytoplasmic labelling is more informative about stress on every row, and the margin is largest where the
stress measure is expression-free (env axes ratio 1.76, environmental archetype 1.91) and smallest where it is
a signature score of the same cell (1.23-1.25; the retired minus-634 arm had shown that the signature-axis
margin is mostly the shared genes, §3). Cluster-level AMI gives the same ordering (stress_env 0.037 vs 0.022,
stress_expr 0.048 vs 0.038), so the difference is not an artefact of the number of labels.

**Analysis C, information in the matched PCs** (`C_summary.csv`, `information_cv.png`,
`information_within_sample.png`, `information_by_depth.png`, `morans_i.png`): ridge regression of each axis on
the 50 sample-corrected PCs, block-CV R2 (pooled z and axes standardised within sample), LOSO, AUROC, Moran's I
on the 30-NN graphs.

| Expression-free axes (20) | Cytoplasmic PCs | Nuclear PCs | Both |
|---|---|---|---|
| median block-CV R2, pooled z | 0.072 | 0.052 | 0.098 |
| median block-CV R2, within-sample z | 0.088 | 0.060 | |
| cyto better (pooled z / within-sample z / LOSO / AUROC / Moran's I) | 19 / 19 / 17 / 18 / 18 of 20 | | |

Largest within-sample R2 (cyto / nuc): fibroblast distance 0.155 / 0.117, stromal kernel 0.145 / 0.109,
fibroblast count 0.134 / 0.096, non-tumor distance 0.129 / 0.099, tumor fraction 0.129 / 0.091, vascular
distance 0.121 / 0.092, vascular kernel 0.117 / 0.090, crowding 0.091 / 0.055; T-cell and immune axes
0.035-0.089 vs 0.026-0.066 (cyto better in 5/5 samples for tumor fraction, non-tumor distance, T-cell count /
distance, the cytotoxic and vascular kernels and the fibroblast axes). Signature axes: hypoxia 0.305 / 0.275,
IFNG 0.305 / 0.264, TNFA/NFKB 0.231 / 0.217, heat shock 0.086 / 0.056 (cyto better), SG markers 0.042 / 0.065
(nuclear better). Nuclear PCs remain better for cell area (0.079 vs 0.102), SG marker expression, granule burden
and the 8-TF rule; both PC sets together beat either alone everywhere. LOSO R2 across samples is
negative for most expression-free axes (the environment-expression link does not transfer between samples) and
clearly positive only for the fibroblast count (0.050 / 0.017) and the hypoxia (0.106 / 0.076) and IFNG (0.151 /
0.114) signatures, marginal for the stromal kernel (0.026 / -0.001), fibroblast distance, heat shock and the 8-TF
rule. Depth strata: the cytoplasmic advantage on the expression-free axes is largest in the shallowest
nuclear-depth decile (+0.05 to +0.15 R2 for the distance, kernel, crowding and core axes, +0.03 to +0.04 for
the neighbour counts, +0.48 for cell area) and shrinks to about 0.00-0.03 from decile 5 upward, where it persists for
the stromal, core and vascular axes and vanishes (or reverses in the deepest decile) for crowding, cell area
and the immune axes.

**Signal strength for question 2** (`comparison/signal_strength.csv/png`): the largest cytoplasmic-minus-
nuclear difference is approach 2 with the environmental archetype (AMI +0.015 [+0.014, +0.016], ratio 1.91);
approach 1 with the expression-free axes gives +0.016 R2 beyond sample [+0.014, +0.021] and AMI ratio 1.76. The
expression-based rows (signature axes 1.23, expression archetype 1.25) show the smallest differences.

Paired visualisations: `sankey_pair_stress_env.png` and `sankey_pair_stress_expr.png` (cytoplasmic substates |
stress labelling | nuclear substates on the cells kept in both labellings; ribbons coloured by the stress node,
the reproducibly enriched pairs of each compartment opaque: around the environmental archetype the cytoplasmic
side carries the Interferon / Complement -> Immune attack, Hypoxia -> Vascular deprivation, Stroma border ->
Stromal contact and ER stress -> Unstressed flows; the nuclear side carries Hypoxia -> Core isolation / Vascular
deprivation, Interferon -> Core isolation (not Immune attack) and EMT/stroma / Stress response -> Stromal
contact; around the expression archetype both sides carry the Interferon, Hypoxia and Inflammatory flows), `spatial_pair_<sample>.png` (cytoplasmic labels | nuclear labels | environmental
archetype | expression archetype; in CC the cytoplasmic Stress response and Hypoxia substates trace the nest
interiors that the Core isolation / Vascular deprivation archetypes mark, while the nuclear labelling paints the
same nests mostly Baseline / Interferon), `correspondence_summary.png`, `stress_subtype_comparison.png`,
`signal_strength.png`, `program_axis_difference.png`.

## 6. Fairness controls and caveats

- **Depth.** Nuclear counts of these genes are lower by construction (median gene-set transcripts per cell
  30-215 cytoplasmic vs 20-138 nuclear, 4.7 % of cells with no nuclear gene-set transcript). The depth strata of
  analysis C show that most of the cytoplasmic advantage on the environmental axes sits in the shallow nuclear
  deciles; the thinned cytoplasmic embedding (binomial downsampling to the nuclear count, re-embedded on HGCC)
  remains the definitive control and has not been run.
- **Circularity.** The expression signatures share genes with the 634 set (§3), so the signature axes and the
  expression archetype are mechanically closer to the cytoplasmic matrix, and their cytoplasmic lead for
  question 2 is not evidence on its own (the retired minus-634 arm showed that lead to be mostly the shared
  genes: AMI ratio 1.04-1.13 without them vs 1.20-1.23 with them in the review run). The expression-free axes and the
  environmental archetypes use no transcript of the tumor cell and are the fair evidence for question 2.
  Spatial-override labels are excluded from their own axes; the Hypoxia substate is validated by the
  expression-free vascular / core axes, not by the hypoxia signature alone.
- **Spatial autocorrelation.** Block bootstrap and block folds throughout; no cell-level tests.
- **Sample heterogeneity under pooling.** Many label-level associations are carried by one or two samples (OC and
  CC for Hypoxia-vessel, CC for the Stress response niche, Prostate for Complement-immune and for the nuclear
  Stress response); per-sample log2 ORs are in the correspondence tables and the forest plots, the
  Mantel-Haenszel estimates next to the pooled ones. LC is the recurrent exception (deep sample, its own Hypoxia
  cluster).
- **Annotation route.** BC, OC, CC annotated from 10x fine labels, LC and Prostate from hand-mapped graph clusters
  (§3.1); the label-free kernels agree with the label-based counts in all five samples (Spearman 0.48-0.86,
`stress_axes/labelfree_vs_label_agreement.csv`), and
  every immune finding is supported by the label-free cytotoxic / myeloid kernels.
- **Pooled z vs sample-centred PCs.** Analysis C is now reported both with pooled z (as in the first run) and with
  axes standardised within sample; the within-sample version is the one to quote for absolute R2.
- **Stress subtypes are composites of the same axes**, so approach 2 is not independent evidence; it is the
  label-against-label reading of the same stress measurements. The dominant-archetype rule assigns each cell one
  stress type, which under-represents cells under several stresses (§4.1, item 5); the two labellings are kept
  separate for the reasons given in §4.1.
- **Asymmetry of the question.** The comparison is "cytoplasmic vs nuclear expression of the target genes", which
  is the model's question; it is not a symmetric contest between compartments on a neutral gene set.

## 7. Answers to the two questions (2026-11-08, streamlined run)

**1. Are the cytoplasmic target-gene substates connected to tumor stress? Yes, coherently but weakly, and the
clearest reading comes from the stress groups across substates on the signature axes.** Every substate has a
reproducible environmental or expression signature that matches its name: Hypoxia cells are far from vessels
and in the tumor core and express the hypoxia programme (41 % of them fall in the hypoxia archetype); Interferon
cells sit closer to T cells, myeloid cells and cytotoxic neighbours and express the IFNG programme (53 % in the
IFNG archetype); the Stress response substate occupies a sparse, vessel-poor, immune-poor nest-interior niche
with high TNFA/NFKB (41 % of its cells in the Core isolation archetype); Complement / acute phase cells are
immune- and myeloid-adjacent and inflammatory; Baseline cells are vessel- and fibroblast-near and unstressed;
proliferating cells are the least hypoxic. Granule-level scores do not separate the substates. Effect sizes are
small (Cramer's V <= 0.40, AMI <= 0.06, a stress axis explains at most 0.115 of the substate variance beyond
sample) and the share of a substate's cells that sits in the matching stress subtype is 25-55 %, so stress
enriches the substates rather than defining them. **Approach to proceed with**: approach 1 (stress groups
across substates) on the canonical signature axes gives the strongest signal (max AMI 0.059, V 0.398); the
environmental archetype (approach 2, `stress_env`) is the companion when an expression-free stress
classification is needed, and the environment x response cross-tab of §5.1 links the two.

**2. Is cytoplasmic expression of the target genes more closely related to stress than nuclear expression?
Yes, on every measure, with the margin largest for expression-free stress and smallest for expression
signatures.** The cytoplasmic labelling has higher AMI with the stress groups on 19 of 20 expression-free axes
(ratio 1.76) and with both explicit stress labellings (environmental archetype ratio 1.91, expression archetype
1.25); its substates explain +0.016 more variance of the environmental axes beyond sample (CI [+0.014, +0.021]);
and the cytoplasmic PCs predict the expression-free axes better than the nuclear PCs on 19 of 20 (within-sample
block-CV R2 0.088 vs 0.060, pooled 0.072 vs 0.052), with higher graph smoothness on 18 of 20. Only the
cytoplasmic labelling recovers the immune neighbourhood of the Interferon substate and the sparse vessel-poor
niche of the Stress response substate. The qualifications stand: the advantage is concentrated in cells with few
nuclear transcripts (0.00-0.03 R2 from the fifth nuclear-depth decile upward, persisting for stromal, core and
vascular axes); on the signature axes it is mostly the shared genes (§3, §6); nuclear PCs are better for cell
area, SG marker expression and granule burden; both PC sets together beat either alone. **Approach to proceed
with**: for the cytoplasmic-vs-nuclear comparison the expression-free measures (approach 2 `stress_env`, the
largest difference, and approach 1 env axes) are the fair and the strongest evidence. The claim "more closely
related to stress" is supported as a counts-plus-localisation effect until the thinned-cytoplasm control is run.

## 8. Implementation (streamlined 2026-11-08)

- `code/1_exploration/stress_utils.py`: axes catalogue (`AXES`, `axis_set`: env / signature / other), reference
  sets, label-based axes (`label_axes`), Gaussian kernels (`kernel_score`, `kernel_axes`), expression axes
  (`expression_axes`, canonical lists), granule axes, pooled groups and blocks; composite stress scores and the
  two stress labellings (`COMPOSITES_ENV`, `COMPOSITES_EXPR`, `composite_scores`, `dominant_archetype`);
  analysis helpers (R2 beyond sample, epsilon2, Cramer's V, AMI incl. from a contingency table, log2 OR pooled /
  Mantel-Haenszel, label-vs-label enrichment with block-bootstrap CIs (`label_correspondence`), ridge
  out-of-fold predictions, within-sample standardisation, Moran's I, Spearman matrices).
- `code/1_exploration/stress_figures.py`: the figure helpers (spatial maps, boxplots in physical units, stacked
  group shares, dot plot without grid, enrichment heatmaps with cell counts, composition bars, forest plot,
  paired bars, Sankey-style alluvial flows between labellings (`alluvial_figure`: matplotlib node bars and
  Bezier ribbons, ribbon width = cells, highlighted ribbons opaque)); substate colours from `subtype_utils.LABEL_COLORS`, a fixed palette for the stress subtypes.
  The "reproducible" flag (|log2 OR or log2 enrichment| > 0.5, block-bootstrap CI excluding 0, same sign in >= 3
  of 5 samples) is a table column and the selection rule of the forest plot; it is no longer drawn as dots on
  the heatmaps.
- `code/1_exploration/6_stress_substates.ipynb` (local, AIVC-env;
  executed with `jupyter nbconvert --execute --inplace`, about 1 h). The output folder is emptied at the start
  of a run; the per-cell stress table is the one file carried over (`REBUILD = True` recomputes it, about 5 min;
  retired minus-634 columns are dropped from the cache on load). Parts: 1 stress table and per-sample summary;
  1b composites, stress labellings, environment x response cross-tab; 2 join with labellings, program scores,
  PCs, graphs; 3 per-compartment analysis and figures; 4 comparison; 5 summaries. An environment variable
  `STRESS_SMOKE=<folder>` runs a 40k-cell development subsample into that folder.
- Output layout `output/1_exploration/stress/`:
  - `stress_axes/`: `stress_scores.csv` (one row per pooled cell: 28 raw axes, their pooled z and groups, the
    spatial block, coordinates, and three supplementary columns that are not axes: `endothelial_n50`,
    `immune_dist`, `nontumor_n100`),
    `stress_axes_summary.csv`, `labelfree_vs_label_agreement.csv`, `composites.csv`, `stress_subtypes.csv`
    (two labellings and the ten composite scores per cell), `stress_subtype_sizes.csv`, `stress_env_vs_expr.csv`;
    figures `axis_distributions_by_sample`, `composite_correlation`, `stress_env_vs_expr`, `sankey_env_vs_expr`,
    `stress_subtype_composition_by_sample`, `spatial_axes_<sample>`, `spatial_stress_subtypes_<sample>`.
  - `cytoplasmic/` and `nuclear/` (identical file sets, question 1 for that labelling): `profile.csv`,
    `correspondence_axes.csv`, `correspondence_labels.csv` (approach 1), `contingency_stress_env.csv`,
    `contingency_stress_expr.csv`, `stress_subtype_correspondence.csv`, `stress_subtype_summary.csv`
    (approach 2), `signal_strength.csv`, `program_correlation.csv`, `summary.txt`; figures `profile_heatmap`,
    `odds_ratio_heatmap`, `dotplot_substate_axes`, `boxplot_axes_by_substate` (15 raw axes incl. the three
    granule-level ones), `stacked_groups_by_substate`, `stress_subtype_enrichment`, `stress_subtype_composition`,
    `sankey_substates_vs_stress`,
    `spatial_<substate>_vs_<axis>` x 3, `spatial_labels_vs_stress_<sample>`, `forest_headline`,
    `program_axis_heatmap`.
  - `comparison/` (question 2): `AB_comparison.csv`, `stress_subtype_comparison.csv`, `information.csv`,
    `information_within_sample.csv`, `information_by_depth.csv`, `morans_i.csv`, `C_summary.csv`,
    `signal_strength.csv`, `summary.txt`; figures `correspondence_summary`, `stress_subtype_comparison`,
    `information_cv`, `information_within_sample`, `information_by_depth`, `morans_i`, `signal_strength`,
    `program_axis_difference`, `sankey_pair_stress_env`, `sankey_pair_stress_expr`, `spatial_pair_<sample>`.
  - `summary.txt`: the answers and the approach with the strongest signal.
- Nothing upstream modified (634 set, labellings, PCs).

## 9. Open items

- Thinned cytoplasmic embedding (binomial downsampling to the nuclear count, HGCC) as the definitive depth
  control; the depth strata suggest the localisation-specific advantage is 0.01-0.03 R2 on the stromal / core /
  vascular axes.
- Why nuclear expression of the target genes predicts SG marker expression and granule burden better than
  cytoplasmic expression (nuclear retention of long SG-marker transcripts?).
- H&E morphology axis for OC (the only sample with H&E).
- Monotonic and oscillatory stress of `plans/AIVC_plan.md` §4.2 remain undefined: the data have no vessel-flow or
  directional-force proxy.
- Label transfer of the BC/OC/CC fine labels to LC and Prostate (one annotation route for all five samples) if
  the immune axes are used further; the label-free kernels currently cover the gap.
- Whether, later, the stress groups should be defined by the model's predicted stress instead of these proxies.

## Log

- 2026-11-07: document created with the analysis design (advisor's request). Scope fixed to the five subtyped
  samples and the pooled cells; expression signatures keep their canonical genes, with a minus-634 sensitivity
  version; SG marker expression added as a stress axis. No code yet.
- 2026-11-07: §3.1 added after checking the annotation route per sample (10x labels for BC/OC/CC, manual cluster
  mapping for LC/Prostate): coarse reference sets, Away/Close as the primary group rule, Unknown/Mixed excluded,
  label-free marker-kernel companions for every contact axis.
- 2026-11-08: analysis implemented (`stress_utils.py`, `6_stress_substates.ipynb`) and run on the pooled cells;
  §5-7 rewritten with the results. Substates correspond to stress coherently but weakly (max R2 0.11 beyond
  sample); cytoplasmic PCs read the environment better than nuclear PCs on 19/20 expression-free axes (median
  block-CV R2 0.072 vs 0.052) but the advantage sits mostly in shallow-nuclear-depth cells; nuclear adds
  independent information (joint 0.098). Partial support for the cytoplasmic outcome; thinned-cytoplasm control
  and within-sample standardised C listed as the next checks.
- 2026-11-08 (review): revised after the user's comments. (1) The canonical and the minus-634 gene lists are run
  and saved as peer variants, not as a sensitivity check; canonical gives the stronger signal. (2) Explicit
  stress subtypes added (environmental / expression archetypes by dominant composite, k-means clusters) and
  evaluated label against label next to the stress-groups-across-substates approach; approach 1 / canonical is
  the strongest signal for question 1, the expression-free measures (env axes, environmental archetypes) the
  fair and strongest evidence for question 2. (3) Outputs split into self-contained folders (`stress_axes/`,
  `cytoplasmic/`, `nuclear/`, `comparison/`), the old flat folder emptied. (4) Intuitive figures added (spatial
  maps, UMAPs, raw-unit boxplots, stacked shares, enrichment heatmaps with counts, composition bars, forest plots)
  and analysis C repeated with within-sample standardised axes (median R2 0.088 vs 0.060 on the expression-free
  axes). `stress_figures.py` added.
- 2026-11-08 (streamlining): the minus-634 arm retired on the user's instruction after the review run showed
  the canonical lists stronger on every metric with the same associations (cytoplasmic approach 1 max AMI 0.059
  vs 0.026, V 0.398 vs 0.270, R2 0.115 vs 0.052; expression archetype AMI 0.045 vs 0.021; question 2 AMI ratios
  1.20-1.23 vs 1.04-1.13); rationale in §3, the circularity it guarded against kept as a caveat in §6. The
  pooled k-means stress clustering dropped (it pooled the environmental and expression features; §4.1 gives the
  reasons for keeping the two labellings separate: exposure vs response, fairness of question 2, dominance of
  the environmental composites with silhouette 0.12, the environment x response cross-tab as an output). The
  Stress granule composite removed from the expression archetype (granule-level scores stay as stress axes
  across the substates, not as subtyping features). UMAP figures dropped, the "reproducible" dots removed from
  the heatmaps (flag kept in the tables), the dot-plot grid removed, `stress_env_vs_expr` added. Notebook re-run
  on the canonical arm (28 axes, two labellings); §5 and §7 rewritten with the new numbers.
- 2026-11-08 (Sankey figures): on the user's suggestion, Sankey-style flow figures added for the label-vs-label
  correspondence (`alluvial_figure` in `stress_figures.py`): environmental archetype | substates | expression
  archetype per compartment, environmental -> expression archetype in `stress_axes/`, and cytoplasmic substates
  | stress labelling | nuclear substates in `comparison/`; reproducibly enriched pairs drawn opaque, all other
  flows faint. Figures only; statistics and summaries unchanged. Notebook re-run.
- 2026-11-08 (audit): code, outputs and documents audited at the user's request. Outputs: the folder holds
  exactly the 103 files of the current notebook, all from the last run (the cached stress table carried over);
  no stale file. Every number quoted in §5 and §7 re-checked against the saved tables; all agree. Corrections to
  this document: the stress-table path (§2), the never-implemented CD8 distance column removed from §3 / §3.1,
  the agreement range (0.48-0.86, §6), the LOSO and depth-strata sentences of §5.4 made precise, the
  supplementary columns of `stress_scores.csv` listed in §8. Code: dead code removed from `stress_utils.py`
  (unused logistic branch of `cv_predictions`, unused imports), `axis_set` made explicit, docstrings updated,
  unused import dropped from the notebook; no numerical change, so the outputs were not regenerated (smoke run
  only). `CLAUDE.md` now lists `code/1_exploration/` and `output/1_exploration/`.
- 2026-11-14 (final audit of the exploratory stage): the three pooled notebooks were re-run after a fix of their batch-mixing / purity summaries; the cluster and
  label tables this analysis reads are byte-identical, so the stress analysis was not re-run.
  `stress_axes/stress_axes_summary.csv`,
  which had been re-saved by a spreadsheet application on 2026-11-10 (CRLF line ends, truncated precision), was
  regenerated from the notebook code; values unchanged. Document: zero-transcript share corrected (4.7 % of cells
  have no nuclear gene-set transcript), reference to a notebook builder script removed. No result changed.
