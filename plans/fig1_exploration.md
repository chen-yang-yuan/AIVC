# Fig. 1 Exploration: Nuclear vs Cytoplasmic Expression in Tumor Cells

Last updated: 2026-09-26.

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

- 2026-09-26: pipeline implemented and verified end to end on LC in quick mode (20k cells; light steps ~2 min,
  heavy steps ~5 min locally) and on CC in quick mode (the shallow sample). Full runs on all six datasets pending
  (HGCC submission of `1_heavy.sh`, then the notebook).

## 5. Measurements (quick-mode, LC, 20k-cell subsample; to be replaced by the full run)

- Depth: median nuclear/cytoplasmic 568/704 per tumor cell; nuclear share of in-cell reads 0.42; 93.5% of cells pass
  the matched-depth filter; tile edge 130 um.
- Positive controls: MEG3 -2.1, MIAT -1.7, PVT1 -1.9, CRNDE -2.6, HOTAIR -1.1 (all nuclear-retained); NORAD +0.32,
  EPCAM +0.37, HSPA8 +0.46, EEF1G +0.62 (cytoplasm side of the offset). 599 genes nuclear-retained, 6 cytoplasm-
  enriched, 3927 balanced. Most reliable per-gene localisation: interferon-response genes (CXCL10, IRF1, CXCL9),
  rel 0.5-0.6.
- SG marker set vs matched control (set level): beta -0.23 vs +0.01, reliability 0.13 vs 0.01.
- Clustering (res 0.5): ARI nuc/cyto 0.28 against ceilings nuc-halves 0.51, cyto-halves 0.64, cross-halves 0.37;
  ARI nuc/total_full 0.41; no cluster met the cytoplasm-only criteria on the subsample.
- DE (nuclear partition): cyto-only 37, nuc-only 406, shared 426 gene-cluster pairs; top-50 marker Jaccard 0.45.
- Association: 1795 genes; largest fraction-model effects come from segmentation method (EEF1G 0.21) and are
  excluded from the candidate axes; pathway fraction model above the matched null in 23 sets for subtype, 7 for
  niche, 10 for tumor fraction.

CC (quick mode, the shallow sample): 75.8% of cells pass the matched-depth filter (median matched depth 42); per-gene
reliability is reportable for only 22 genes at cell level, so CC will rest on tile-level and pathway-level statistics;
1607 genes classify as nuclear-retained at the offset (to be checked against segmentation strata in the full run);
positive controls remain nuclear-retained (MEG3 -1.3, CRNDE -2.1, HOTAIR -1.5); half-split clustering ceilings are at
the noise floor (ARI 0.03-0.04), which is the expected honest result at this depth. Light steps take ~3.5 min
(quick) and the heavy quick steps ~2 min.

### Annotation QC results (step 0b, all cells, 2026-09-26)

Provenance and tumor-set purity (`output/2_exploration/cross_sample/annotation_provenance.csv`):

| Sample | Source | Tumor frac | Unresolved | Non-malignant epithelium class | Malignant cells confidently epithelial (transfer) | Immune-like | Benign-like malignant clusters (cells) |
|---|---|---|---|---|---|---|---|
| BC | 10x supervised | 0.25 | 0 | yes | 0.97 | 0.018 | n/a |
| OC | 10x supervised | 0.49 | 0 | yes | 0.99 | 0.022 | n/a |
| CC | 10x supervised | 0.38 | 0 | yes | 0.96 | 0.041 | n/a |
| LC | manual clusters | 0.16 | 0.03 | yes | 0.92 | 0.030 | 2 of 5 (20%): cluster 31 basal-high (656 cells), cluster 14 lineage-low, immune-high (8.4k) |
| Prostate | manual clusters | 0.49 | 0.11 | no | 0.76 | 0.018 | 2 of 8 (28%): cluster 9 TP63-high, AMACR-low = benign basal/glands (10.5k); cluster 4 low signal on every set (15.8k) |
| Skin | manual clusters | 0.44 | 0.13 | no | 0.64 | 0.011 | 5 of 10 (22%): clusters 20/26/27/28 keratinocyte-like, TP63-high, melanocytic markers absent (3.1k); cluster 5 near-zero on every set (7.8k) |

- No sample has lineage-negative tumor cells by the per-cell rule (the sets are broad), and immune or endothelial
  contamination of the tumor sets is at most 4% by markers and at most 3.6% by the transfer classifier. The
  supervised sets are clean by every check.
- The manually annotated sets carry benign or low-signal clusters inside "Malignant cell": ~20% (LC), ~28% (Prostate),
  ~22% (Skin) of their tumor cells. Prostate cluster 9 and the Skin keratinocyte clusters are the clearest cases;
  Prostate cluster 4 and Skin cluster 5 look like low-count cells rather than a cell type.
- Non-epithelial label transfer (classifier trained on BC/OC/CC; leave-one-dataset-out accuracy 0.6-0.9 for most
  classes, ~0.75 median): endothelial, B, mast, pericyte and fibroblast labels agree at 0.7-0.97 in LC/Prostate/Skin;
  the T-cell subtypes and dendritic vs myeloid distinctions do not transfer well (0.2-0.5), as expected across
  tissues. Prostate "Unknown" (21k cells) is predicted with low confidence for any class, i.e. mostly low-quality
  cells; LC and Skin "Mixed" are ~60% confidently immune.
- Implication for the exploration: results on LC, Prostate and Skin should be read with the benign-like clusters in
  mind; the nuclear "subtype" partitions there may separate benign from malignant cells. The fix, if adopted, is a
  harmonised tumor set (label AND not benign-like cluster AND marker/transfer pass), which would drop 20-28% of the
  tumor cells in those three samples. Not applied (decision: report only).
- Caveats: the off-lineage flag is noisy (several off-lineage genes, e.g. MSLN, TP63, SOX2, CDKN2A, are expressed in
  other carcinomas); melanoma cells map to the classifier's "epithelial" class because that class is the only
  non-stromal, non-immune class it knows.

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
