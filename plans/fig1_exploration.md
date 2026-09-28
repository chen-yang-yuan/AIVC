# Fig. 1 Exploration: Nuclear vs Cytoplasmic Expression in Tumor Cells

Last updated: 2026-09-27.

## 1. Purpose and scope

- Exploratory analysis of the six processed Xenium 5K FFPE samples (BC, OC, CC, LC, Prostate, Skin), tumor cells
  only, that motivates AIVC's compartment-level target and selects the candidate target gene set. Specification:
  `plans/AIVC Model Design Tracker.md`, "Next step: exploratory analysis for Fig. 1".
- Narrative: (A) the compartments are different; (B) the difference is biology rather than artifact; (C) the
  difference relates to the model's inputs (subtype, morphology, neighborhood).
- Outputs of the stage: candidate target gene set (or pathway set), the choice among target forms A/B/C, the
  identity source (total vs nuclear), the Fig. 1 panels.
- This document reflects the latest state only; locked decisions and notable updates go to the Log (section 7).

## 2. Implementation

Code: `code/2_exploration/` (execution order and file inventory in its `README.md`). Outputs:
`output/2_exploration/<dataset>/` (tables, `fig/`), `output/2_exploration/cross_sample/`; development runs on a
20k-cell subsample go to `<dataset>_quick/`.

- Light steps (0 QC, 1 compartment statistics, 2 concordance, 5 controls, 8-light pathways, 9 candidates,
  cross-sample) run in `fig1_exploration.ipynb` locally (`AIVC-env`) on all tumor cells.
- Heavy steps (3 clustering, 4 DE, 7 association, 8 pathway association) run on all kept tumor cells on HGCC via
  `1_heavy.py` / `1_heavy.sh` (SLURM array over the six datasets, `preprocessing-env`) and are read by the notebook
  after `make out-pull`. The same `pipeline.py` step functions serve both; `QUICK = True` runs everything locally on
  the subsample.

Inputs: `processed_data/{nuclear,cytoplasmic}_expression_matrix.npz` (tumor cells x 5001 genes; nuclear =
`in_nucleus == 1`, cytoplasmic = in cell and outside the nucleus; nuclear + cytoplasmic = in-cell total),
`cell_ids.npy`, `gene_ids.npy`, `intermediate_data/adata.h5ad` (obs: coordinates, areas, segmentation method,
cell types of all cells). Caveat: the matrices contain all decoded transcripts (no QV >= 20 filter; decision
2026-09-26 to keep them and record it).

## 3. Design

Definitions (cell i, gene g): k = cytoplasmic count, n = k + nuclear count (binomial trials), M_i and Q_i =
cytoplasmic and nuclear depth over the panel, d_i = min(M_i, Q_i).

Design rules applied everywhere:
1. Depth matching. Per-gene localisation statistics use all cells with a per-(cell, gene) leave-one-gene-out offset,
   `logit((M_i - k + 0.5) / (M_i - k + Q_i - (n - k) + 1))`, so a gene's log odds ratio (beta_g) is its departure from
   the cell's global nucleus/cytoplasm balance. Clustering, DE and cell-level concordance use matrices rarefied to d_i
   per cell (`nuc_m`, `cyto_m`, `total_m`), plus `total_full` at the true depth, on cells with one nucleus and
   d_i >= 30; complementary half-splits (`nuc_h1/h2`, `cyto_h1/h2`) by binomial thinning give disjoint reads at equal
   depth.
2. Reliability against counting noise. Per gene, `rel = 1 - Var(y_sim) / Var(y_obs)` with y the offset-corrected
   empirical logit on units with n >= 10 and y_sim from a parametric binomial null (B = 5); NA below 200 units. The
   analytic `1/(n p (1-p))` formula is kept as a cross-check column (it overestimates reliability at low n).
3. Positive controls. MALAT1, NEAT1 and XIST are not on the 5K panel. Nuclear-retained controls: MEG3, MIAT, PVT1,
   CRNDE, HOTAIR (all classified nuclear-retained in LC, beta -1.1 to -2.6). Cytoplasmic controls: NORAD plus the
   most abundant mRNAs (EPCAM, HSPA8, EEF1G ...). Strata checks: beta per segmentation method and per nucleus/cell
   area-ratio quintile, with the range across strata as an artifact indicator.
4. Tumor cells only, per sample first, then across samples.

Analyses:
- A1 Global description: per-cell nuclear fraction (overall, by segmentation method, by area ratio); per-gene beta_g
  with quasi-binomial CI vs abundance; classes nuclear-retained (beta < -log 2, CI excludes 0), cytoplasm-enriched
  (beta > log 2), balanced, low-coverage (< 200 transcripts).
- A2 Concordance: at cell level (genes with mean >= 0.5) and spatial-tile level (non-overlapping tiles holding ~15
  tumor cells; genes with mean >= 0.05): cross-compartment correlation r_nc on depth-normalised log values,
  within-compartment split-half reliabilities rel_n, rel_c (Spearman-Brown), disattenuated r_true = r_nc /
  sqrt(rel_n rel_c), reported where both reliabilities >= 0.1. Divergent gene: low r_true and in the lowest decile
  of its abundance bin.
- A3 Clustering: identical pipeline per matrix (shared gene set detected in >= 0.5% of cells, no HVG selection,
  normalize_total to the median depth, log1p, PCA 30, kNN 15, igraph Leiden at resolution 0.5 and 1.0) for the eight
  matrices. ARI/NMI for all pairs; the honesty ceiling is the within-compartment half-vs-half agreement against the
  cross-compartment half-vs-half agreement. Identity-source answer: ARI(nuc_m, total_full) and per-cluster recovery
  (best F1) of the total_full clusters by the nuclear partition. Cytoplasm-only cluster: max Jaccard < 0.3 with every
  nuclear and total cluster and reproduced in both cytoplasmic halves (Jaccard > 0.5). Spatial purity (k = 10 tumor
  neighbours) vs permutation; per-cluster enrichment of neighbourhood and morphology features (z vs permutation);
  per-cluster median depth.
- A3b Cytoplasm beyond nucleus (the "persuasive panel" test). Log-normalised `cyto_h1` is residualised on a nuclear
  design built from `nuc_h1`: nuclear PCs with split-half reliability >= 0.2 (errors-in-variables corrected with
  Sigma_e = (1 - rel) var), the `nuc_h1` cluster one-hot, a cubic spline of log depth, log cell area, nucleus/cell
  ratio and segmentation method. The `nuc_h2` half residualised on the same design is the empirical leakage floor
  (leakage + counting noise only); `nuc_h1` residualised on the analogous cytoplasmic design is the descriptive
  converse. Residual PCA (30); `cyto_h2` projected onto the `cyto_h1` residual loadings gives split-half reliability.
  Number of real cytoplasmic residual dimensions = singular values above the floor's top singular value.
  Primary test on residual PCs 1-3: reliability >= 0.3, Moran's I (k = 10 tumor kNN, scores centred within strata)
  >= 0.05, above the 99th percentile of a permutation WITHIN strata (nuclear cluster x depth quintile x
  segmentation method, min 20 cells) and above the floor, and a dose-response across tumor-fraction quintiles
  within strata of |effect| >= 0.3 SD with tile-bootstrap CI excluding 0 in both halves (criterion A). Residual
  clusters (Leiden on the residual scores; replicated on the projected `cyto_h2` scores) are eligible for the panel
  (criterion B) if: size >= max(200, 1%); repro F1 >= 0.6; median depth >= 0.5 x sample median; segmentation share
  <= 0.9 unless the sample's is >= 0.8; not morphology-driven (|d| < 0.5 SD on cell area and ratio); stratified
  purity excess >= 0.10 and above the floor's maximum; |d(tumor_frac)| >= 0.3 SD (strata-weighted standardised
  difference, tile-bootstrap CI excluding 0, same sign in `cyto_h2`); >= 5 cyto-only DE genes from a within-stratum
  paired DE on `nuc_h2` vs `cyto_h2`. Ranking by the lower CI bound of |d(tumor_frac)|, min over halves; immune
  density labels the cluster ("immune-rich") but does not rank. Verdicts: "panel" (A and B), "continuous axis only"
  (A), "cluster only (score fails)" (B), "no cytoplasmic structure beyond nucleus". The main-figure sample is the
  qualifying sample with the largest rank statistic; the zoom inset is the 600 um window holding the most cells of
  the rank-1 cluster. The same stratified statistics are reported for the plain `cyto_m` clusters (the nested
  version of A3) and for the floor and converse partitions.
- A4 DE: partition on `nuc_h1`, test on `nuc_h2` vs `cyto_h2` (no double dipping, equal depth): pseudobulk log2 FC
  per cluster in each compartment, delta = lfc_cyto - lfc_nuc with a tile block-bootstrap CI (200 reps); cyto-only DE
  = |lfc_cyto| > 1, |lfc_nuc| < 0.5, CI excludes 0 (nuc-only conversely). Wilcoxon top-50 marker overlap as the
  conventional panel. Converse with the `cyto_h1` partition.
- B5 Controls: abundance-bin (20 bins) percentile of every per-gene statistic; strata checks; SG marker set vs an
  expression-matched control set at set level.
- B6 Cross-sample: recurrence of candidate criteria per gene across samples; pairwise scatter of beta and r_true
  between samples; recurrence of cytoplasmic clusters by centroid correlation.
- C7 Association: categorical axes on the kept cells: subtype (nuclear-half partition; also cytoplasmic and total
  partitions), niche (k-means, k = 8, on Gaussian-kernel cell-type composition over all cells, sigma 30 um), crowding,
  tumor-fraction (boundary vs core), immune-fraction and distance-to-non-tumor quintiles, morphology (nucleus/cell
  ratio and cell-area quintiles, segmentation method). Three models per gene (mean >= 0.2) and axis: cytoplasmic
  fraction (offset-binomial), cytoplasmic count and nuclear count (Poisson with log depth). Deviance explained beyond
  two label permutations; count-weighted SD of level effects; neighbourhood beyond subtype from the subtype x axis
  cross-classification. Cytoplasm-specific sensitivity = dev(cyto) - dev(nuc) and dev(fraction).
- C8 Pathways: count aggregation over the 106 GMT sets plus the SG marker set, run through the same machinery
  (A1, A2 at tile level, C7), each against 50 abundance-matched random gene sets.
- C9 Candidates: per dataset, gene flags reliable (rel >= 0.1), divergent, cyto-only DE, axis-sensitive (fraction
  model or cyto-minus-nuc excess >= 0.005 on subtype, niche, tumor fraction, immune fraction, morphology ratio),
  control flags (reliability in the top decile of its abundance bin; segmentation-strata range not in the top
  decile), annotation (GMT membership, SG set, stress-related sets). candidate = reliable and segmentation-ok and
  (divergent or cyto-only DE or axis-sensitive). Cross-sample: recurrent = candidate in >= 3 samples.

Morphology in this stage is segmentation-derived only (H&E exists for OC only; UNI deferred).

### Annotation QC (step 0b, report only)

All six datasets share the 18-category `cell_type_merged` vocabulary but not the annotation process: BC, OC and CC
carry 10x supervised labels (tumor = the 10x "Tumor" class; a non-malignant epithelial class exists), while LC,
Prostate and Skin were annotated by hand from 10x graph clusters (5, 8 and 10 clusters called "Malignant cell");
Prostate and Skin have no non-malignant epithelial class, and 11% (Prostate, "Unknown") and 10% (Skin, "Mixed") of
cells are unresolved. Step 0b scores every cell on marker sets (epithelial, basal, cancer-type lineage, off-lineage,
proliferation, immune, endothelial, oncogenic; keratins and PTPRC are absent from the panel), profiles each class and
each manually called malignant cluster (benign-like and contaminated flags), computes purity flags on the tumor set,
and transfers the non-epithelial labels with a classifier trained on BC/OC/CC (leave-one-dataset-out accuracy as the
reference) to estimate agreement on LC/Prostate/Skin and non-epithelial contamination of their tumor sets. Decision
2026-09-26: report only; the tumor definition stays the label until the report has been read. Results: section 5.

## 4. Status

- 2026-09-27: full run complete on all six datasets (heavy steps on HGCC, all kept tumor cells; light steps and
  figures locally; notebook executed end to end without errors). All tables and figures are under
  `output/2_exploration/<ds>/` and `output/2_exploration/cross_sample/` (`decision_summary.csv` is the one-page view).
  Interpretation and recommendations in section 5.3; nothing is locked yet.

## 5. Results (full run, 2026-09-27)

### 5.1 Per sample

| | BC | OC | CC | LC | Prostate | Skin |
|---|---|---|---|---|---|---|
| Tumor cells / kept (one nucleus, matched depth >= 30) | 102k / 81% | 160k / 90% | 221k / 76% | 45k / 94% | 95k / 73% | 49k / 77% |
| Median matched depth; nuclear share of reads | 100; 0.54 | 170; 0.53 | 62; 0.45 | 447; 0.42 | 93; 0.46 | 208; 0.37 |
| Genes with mean in-cell count >= 0.5 | 88 | 164 | 21 | 859 | 86 | 311 |
| A1 nuclear-retained / cytoplasm-enriched genes | 2169 / 21 | 1918 / 0 | 1607 / 1 | 599 / 6 | 1379 / 0 | 372 / 20 |
| A1 genes with reliable localisation (rel >= 0.1, >= 200 units) | 15 | 26 | 21 | 209 | 26 | 92 |
| B5 SG set vs matched control, reliability | 0.14 vs 0.12 | 0.07 vs 0.06 | 0.17 vs 0.14 | 0.13 vs 0.01 | 0.06 vs 0.07 | 0.06 vs 0.03 |
| A2 tile-level r_true median (genes) | 0.46 (1419) | 0.66 (981) | 0.60 (1034) | 0.66 (2162) | 0.55 (999) | 0.64 (1648) |
| A3 ARI nuc/cyto vs ceilings nuc-halves, cyto-halves, cross-halves | 0.15 vs 0.20, 0.16, 0.13 | 0.22 vs 0.18, 0.25, 0.19 | 0.06 vs 0.14, 0.13, 0.12 | 0.40 vs 0.47, 0.44, 0.37 | 0.36 vs 0.16, 0.28, 0.25 | 0.12 vs 0.21, 0.22, 0.13 |
| A3 ARI nuc/total_full; cyto/total_full; cyto-only clusters | 0.14; 0.24; 0 | 0.14; 0.16; 0 | 0.12; 0.06; 0 | 0.44; 0.31; 0 | 0.33; 0.39; 0 | 0.15; 0.13; 0 |
| A3b verdict (residual PC1 gene; Moran's I vs floor) | none (XBP1; 0.07 vs 0.05) | continuous axis (H19; 0.26 vs 0.20) | none (EEF1G; 0.07 vs 0.06) | none (CTSH; 0.02 vs 0.03) | none (EEF1G; 0.08 vs 0.04) | none (S100A1; 0.01 vs 0.02) |
| A4 genes cyto-only / nuc-only / shared DE across nuclear subtypes | 190 / 198 / 342 | 468 / 321 / 850 | 268 / 269 / 243 | 24 / 328 / 342 | 142 / 127 / 414 | 73 / 110 / 1022 |
| C7 genes whose cytoplasmic fraction tracks subtype / niche / tumor fraction (dev. expl. beyond permutation >= 0.005) | 15 / 2 / 1 | 2 / 0 / 0 | 31 / 0 / 0 | 842 / 6 / 1 | 24 / 2 / 0 | 209 / 41 / 15 |
| C7 genes with niche signal beyond subtype | 1 | 0 | 0 | 6 | 0 | 1 |
| C8 pathways with reliability above the matched null (95th pct) | 21 | 8 | 27 | 25 | 12 | 18 |
| C9 candidate genes | 6 | 4 | 7 | 147 | 17 | 66 |

### 5.2 Across samples

- Per-gene localisation (log odds ratio) is reproducible between samples: r = 0.80-0.83 among BC, OC and CC, 0.72
  between LC and Skin, 0.35-0.65 otherwise. The two groups coincide with the two panel runs (BC/OC/CC: "5K with
  Cell Typing add-on", one run day; LC/Prostate/Skin: "5K Pan Tissue"), so part of the between-sample difference is
  batch. Positive controls hold in every sample (MEG3, MIAT, PVT1, CRNDE, HOTAIR nuclear-retained).
- Candidate genes: 211 in at least one sample, 28 in two, 6 in three or more (EEF1G, LDHA, NOTCH2NLA, H3F3B, YWHAZ,
  NDRG1). These are abundant genes whose localisation tracks the nuclear subtype; only NDRG1 and LDHA are
  stress-annotated. Genes cyto-only DE in >= 3 samples: 15 (NR4A1, HYOU1, SORD, PIK3R2, SRSF2, PIM1, DAXX, ...).
- Pathways whose cytoplasmic localisation is reliably variable in >= 4 samples: EMT, external encapsulating
  structure organisation, TNFa/NF-kB, RNA processing, RNA polymerase II transcription, interferon gamma response,
  E2F targets, RNA splicing, allograft rejection, mRNA metabolic process, and the SG marker set (4 of 6 samples above
  its matched null). Pathways whose cytoplasmic fraction tracks niche or tumor fraction beyond the matched null in
  several samples: KRAS signalling up (7 sample-axis hits), EMT (5), ER-stress response (4), TNFa/NF-kB (4), cellular
  response to chemical stress (4); effects are small (deviance explained 0.001-0.01).
- Step 3b: no sample yields an eligible cytoplasm-specific cluster. OC passes the continuous-axis criterion: residual
  PC1 (top gene H19) has split-half reliability 0.67, Moran's I 0.26 above the 99th permutation percentile and the
  leakage floor (0.20), and a core-vs-boundary dose-response of +0.6 SD in both halves; the map shows lobe-scale
  regions. Caveats: the floor is high (leakage), and in BC and CC no nuclear PC survived the reliability gate at half
  depth (median depth 100 and 62), so their residualisation rests on the cluster one-hot alone and the verdicts there
  are depth-limited rather than negative. Skin residual PC2 (EEF1G) shows a boundary-vs-core effect of +0.4 SD in
  both halves with I above the floor, failing only the singular-value gate (also under a scale-free version).
- Annotation QC (section 3, run 2026-09-26): the supervised tumor sets are clean; LC, Prostate and Skin carry
  benign-like or low-signal clusters (20-28% of their tumor cells).

### 5.3 Interpretation and recommendations (not locked)

1. The compartments are different, and the difference is reproducible: ~400-2200 genes per sample are nuclear-
   retained relative to the cell's own offset, with lncRNA controls at the expected end, and per-gene log odds
   ratios correlate 0.8 between samples of the same panel run. This is the Part A1 panel and it is solid.
2. Per-cell localisation is mostly counting noise. Reliable per-gene localisation exists for 15-26 genes in the four
   shallower samples and 92-209 in LC and Skin; tile pooling extends coverage to ~1000-2000 genes with median
   disattenuated nuclear-cytoplasmic correlation 0.5-0.66, i.e. substantial gene-level divergence between
   compartments once noise is removed. The SG set is more variable than its matched control only in LC.
3. Compartment-specific clustering does not exist in these data: nuclear and cytoplasmic partitions agree with each
   other as well as each agrees with its own half-split, and no cluster is cytoplasm-only. Nuclear-only clustering
   recovers the total-depth partition as well as cytoplasmic clustering does (identity-source question), but
   clusterability is limited by depth, not by compartment: at matched depth even total expression gives only 3-9
   clusters, while full-depth total finds up to 19 (OC). Recommendation: identity from nuclear reads (target form C)
   is defensible where nuclear depth is adequate (LC, Skin, OC); in BC, CC and Prostate a total-based identity with
   target genes excluded (form A) or a hybrid is safer.
4. Cytoplasm-specific differential expression across nuclear subtypes is real and large (24-468 genes per sample,
   symmetric with nucleus-only DE) and 15 genes recur in >= 3 samples. This is the strongest "post-transcriptional
   signal" panel (A4) and the best per-gene target source.
5. The environment adds little to per-gene localisation: beyond subtype, niche or boundary explain >= 0.5% deviance
   for at most 6 genes per sample; the naive-model conclusion generalises from the SG set to the panel. Where the
   environment shows, it shows at pathway level (KRAS up, EMT, ER stress, TNFa) and in the cytoplasmic COUNT rather
   than the fraction. Recommendation on target form: primary target = cytoplasmic count (A or C), with the
   compartment split (B) as a secondary readout for the subtype-sensitive genes, and count-aggregated pathway
   activity for the environment-sensitive programs.
6. The single persuasive panel (a cytoplasm-only, spatially coherent cluster at the boundary or in an immune niche)
   was not found under the pre-registered rule. The OC H19 residual axis is the closest and can be shown as a
   continuous map with its dose-response, labelled as such. Fig. 1 should therefore lead with A1 (compartment
   structure and controls), the reliability ceiling, A4 (cyto-only DE across subtypes), and C8 (pathway-level
   environment associations), with A3/A3b reported as the honest negative.
7. Target gene set: the recurrent candidates (6 genes in >= 3 samples, 28 in >= 2) plus the 15 recurrent cyto-only DE
   genes form a first list of ~40 genes; it is dominated by abundant housekeeping-like mRNAs (EEF1G, YWHAZ, H3F3B,
   LDHA), so a pathway-level target remains the recommended primary form for the model, with these genes as the
   per-gene track.

## 6. Open items

- Harmonise the tumor definition? Options after reading the step 0b report: (a) keep the label; (b) label AND marker
  QC pass (drop lineage-negative, immune-like and endothelial-like cells, and benign-like clusters); (c) add an
  inferCNV-style copy-number score (gene coordinates from Ensembl, immune/stromal reference) as a label-independent
  malignant call, mainly for Prostate and Skin.

- Run the heavy steps on HGCC for all six datasets and the notebook on the results; fill sections 4-5 with the full
  numbers and write the interpretation (target form, identity source, gene vs pathway).
- Decide the Fig. 1 panel selection from the produced figures; the A3b verdict per sample (decision summary) says
  whether the "cytoplasm-only spatially coherent cluster" panel exists and in which sample.
- Threshold review after the full run: D_MIN, MEAN_MIN_*, ASSOC_MIN_EXCESS, R_TRUE_DIVERGENT, Jaccard thresholds.
- Step 3b uses unscaled log values (as A3 does), so a residual PC can be dominated by one abundant, regionally
  expressed gene (LC quick: residual PC1 is GRP). If that recurs in the full run, add a gene-scaled residual PCA as a
  sensitivity analysis.
- Step 3b quick-mode results (LC, 20k cells): 13 nuclear PCs kept (reliability 0.94/0.88/0.81 for PCs 1-3), EIV
  correction applied with mean shrink 0.6, residual scores nearly orthogonal to nuclear half 2 (max canonical
  correlation 0.14, ARI 0.03), noise-control Moran's I 0.00; verdict "no cytoplasmic structure beyond nucleus"
  (no residual dimension above the leakage floor at half depth on 20k cells). To be re-evaluated on all cells.

## 7. Log

- **2026-09-27** — Full run on all six datasets completed (HGCC heavy steps, local light steps and figures); results
  in section 5. Headline: compartments differ reproducibly per gene; no cytoplasm-only cluster in any sample under
  the pre-registered rule (OC has a continuous residual axis, H19); cytoplasm-specific DE across nuclear subtypes is
  substantial (24-468 genes per sample, 15 recurrent); the environment explains little per-gene localisation and
  shows mainly at pathway level. Recommendations in 5.3; decisions pending.

- **2026-09-26 (evening)** — Step 0b annotation QC added and run on all six datasets (report only). Finding: the
  supervised tumor sets (BC, OC, CC) are clean; the manually annotated sets contain benign-like or low-signal
  clusters inside "Malignant cell" (LC ~20%, Prostate ~28%, Skin ~22% of tumor cells). Harmonisation deferred.
- **2026-09-26 (later)** — Step 3b added after the question whether cytoplasm-defined clusters had been checked for
  spatial pattern: the original A3 tests used global permutation nulls and a whole-cluster Jaccard rule, which restate
  the nuclear subtype and cannot produce the panel when the cytoplasm splits a subtype spatially. The residual design,
  leakage floor, stratified nulls, tile-bootstrap effects and the eligibility / panel rule above are locked before the
  HGCC run. The cell-level "discordance" idea was rejected (with ARI ~0.3 most cells are discordant by depth).

- **2026-09-26** — Document created. Design locked as in section 3; implementation in `code/2_exploration/`
  verified on LC and CC in quick mode. Decisions: morphology = segmentation-derived features only; existing
  compartment matrices kept (no QV filter) with the caveat recorded; heavy steps on HGCC, light steps in the
  notebook; positive controls = panel lncRNAs (MEG3, MIAT, PVT1, CRNDE, HOTAIR nuclear; NORAD cytoplasmic).
