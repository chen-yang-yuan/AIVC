# AIVC Project Plan: Background, Data, and Model Design

Last updated: 2026-11-14.

**Status note.** Sections 2 and 4 describe the original design from the ideas document (stress-score and
pathway heads for nucleus, cytoplasm and RNA granules). The design has since been reframed as a predictive,
interpretable model of cytoplasmic activity of a data-selected gene set; the current state of every design
component, the exploration results and the decision log live in `plans/AIVC Model Design Tracker.md`, and the
original design below is kept until those decisions are locked and folded back in.

## 1. Purpose and scope

- This document covers the project's background and rationale, the data, and the model design and architecture.
  It is refined through brainstorming and always reflects the latest thinking; stale content is removed on update
  and locked decisions are recorded in the **Log** (section 6).
- Pipeline and code-level status live in `plans/data_processing.md`; the exploratory stage is documented in
  `plans/fig1_exploration.md`, `plans/subtype_exploration.md` and `plans/stress_substate_exploration.md`. Later
  stages (model implementation, simulation studies, real data analysis) will get their own documents.
- Source material: `plans/AIVC_ideas_12-08-2026.docx`, preliminary ideas from the advisor and Chenyang (the
  filename date is likely a typo for Dec 2025, which matches the git history).

## 2. Background and rationale

**AIVC** (AI-based Virtual Cell) is framed as a universal multi-modal biological state representation plus a set of
"virtual instruments" (VIs): neural networks that decode or manipulate that representation, acting as computational
twins of lab experiments.

This project builds a **virtual tumor cell model of subcellular responses to environmental stress**, trained on
subcellular-resolution Xenium 5K data. Three goals:

1. Predict per-cell **stress scores** (hypoxia, mechanical, monotonic, immune attack, optionally oscillatory).
2. Predict the **subcellular response**: pathway activity in nucleus, cytoplasm, and RNA granules (originally
   framed as the 50 Hallmark pathways).
3. Act as a **virtual stress simulator**: manipulate the spatial-context input and observe how predicted stress
   and response change, without altering the cell's own expression.

Downstream: apply the model to large scRNA-seq or H&E cohorts (e.g. TCGA) to build a virtual population with
subcellular resolution for population-scale clinical discovery.

Biological motivation: cancer cells face hypoxia, nutrient deprivation, mechanical and osmotic stress, oxidative
stress, and drug stress; they adapt via integrated stress response (ISR) activation, EMT/MET plasticity, senescence,
and aneuploidy. Mechanical stress relates to morphology (actomyosin genes such as ACTG1). Stress responses are
expected to show up as compartment-specific transcript localization (nuclear retention, cytoplasmic translation,
RNA aggregates such as stress granules), which is what subcellular-resolution spatial transcriptomics can resolve.

## 3. Data

Inventory from the ideas doc ("25 samples across 8 cancer types"):

| Source | Samples | In repo? |
|---|---|---|
| Public Xenium 5K | OC (FF), OC (FFPE), LC, BC, CC, Skin melanoma, Prostate; lymph node dropped | 6 processed, all FFPE: `Xenium_5K_{BC,OC,CC,LC,Prostate,Skin}` (OC = `Xenium_Prime_Ovarian_Cancer_FFPE`) |
| Linghua lung Xenium 5K | 3 AAH, 3 AIS, 6 LUAD; 6 patients; clinical info | 12 empty folders `data/Linghua_P*_{AAH,AIS,LUAD}` |
| NC paper Xenium 5K | OV, HCC, COAD (FFPE); 3 patients; clinical info | Not present |
| NC paper CosMx 6K | OV, HCC, COAD (FFPE) | Not present |

Processed public datasets (cells after cleaning; tumor = `cell_type_merged == "Malignant cell"`):

| Dataset | Cells | Tumor cells |
|---|---|---|
| BC | 402,871 | 102,180 |
| OC | 327,607 | 160,250 |
| CC | 588,646 | 221,355 |
| LC | 278,328 | 44,624 |
| Prostate | 193,000 | 95,429 |
| Skin | 112,551 | 49,126 |
| Total | 1,903,003 | 672,964 |

Cell typing:
- BC, OC, CC: 10x cell-type annotations mapped to a merged vocabulary (`data/_utils/cell_type_dict.pkl`).
- LC, Prostate, Skin: 10x graph clusters assigned by hand (`data/_utils/cell_type_dict_manual.pkl`).
- Merged vocabulary: Adipocyte, B cell, CD4+/CD8+ T cell, T cell, Dendritic, Endothelial, Lymphatic endothelial,
  Epithelial (non-malignant), Fibroblast (CAF), Malignant cell, Mast, Mesothelial, Myeloid, Pericyte, Smooth
  muscle, Mixed, Unknown.

Gene panel: the intersection of the six public panels (5,001 genes, `data/_utils/shared_genes.npy`). Adding the
Linghua and NC-paper datasets will require recomputing the intersection or handling missing genes (CosMx 6K in
particular has a different panel).

Modalities available per cell: 5K expression, transcript coordinates with nucleus overlap, H&E image (aligned via
`HE_alignment.csv`; alignment verified visually on OC only), spatial neighbors and their cell types.

## 4. Model design

### 4.1 Inputs and targets

- Inputs: own expression (5K), H&E patch centered on the cell, and optionally neighbor expression as the
  spatial-context signal. Any subset may be missing at inference (e.g. scRNA-seq without images, H&E without
  expression).
- Targets: (a) stress scores, one per stress type; (b) subcellular response, one vector of pathway scores per
  compartment. Only tumor cells are training samples. Cross-sample comparability requires per-cell normalization;
  the ideas doc suggests checking key genes (KRT family, EPCAM) across samples.

### 4.2 Stress-type definitions

| Stress | Definition in ideas doc | Current operational definition | Status |
|---|---|---|---|
| Hypoxia | 8 marker genes; score = k-th highest marker expression (k = 2/3/4), scaled 0–1 | k-th highest of HIF1A, EPAS1, NFE2L2, CREB1, RELA, RELB, NFKB1, NFKB2 on log-normalized expression; median split → Low/High | Implemented (not 0–1 scaled) |
| Heat shock | not in doc | same scheme, 11 HSP genes | Implemented (extra) |
| Mechanical | overall force in neighbor graph; high at tumor core | inverse mean distance to 10 nearest tumor cells, z-scored and clipped; median split | Implemented as crowding proxy |
| Monotonic | directional force in neighbor graph; high at tumor boundary | none | Not implemented |
| Immune attack | proximity to activated CD8 T / NK / inflammatory macrophages | T-cell count within 100 µm, optional Gaussian kernel weighted by GZMA/GZMB/GZMK; also broad immune-cell count within 50 µm | Implemented (T cell + immune proximity) |
| Oscillatory | vessel-based (optional, from Jiahui's work) | none | Not implemented |

The table is the state of `code/2_outcome/5_score.ipynb`. The exploratory stress analysis
(`plans/stress_substate_exploration.md` section 3) redefined the axes for the five carcinoma samples: 28 axes (20
expression-free environmental axes, including vascular, core / boundary, stromal and myeloid ones; 5 whole-cell
signature scores; 3 others), crowding computed on all cells instead of tumor cells only, and hallmark signature
scores instead of the 8-TF rule as the hypoxia measure (the rule is kept as a secondary axis). In `5_score.ipynb`
the T-cell count uses a 50 um radius; the 100 um radius applies to the GZM-weighted kernel only.

### 4.3 Outcome variable (OPEN)

Current candidate (not locked): the sample-relative log cytoplasmic activity of the 634 target genes and of the
recurrent substate programs, with a count link (cytoplasmic depth offset, per-gene per-sample intercept); see
`plans/subtype_exploration.md` section 7 and the Tracker. It is a refinement of direction G2 below. The rest of
this section is the September analysis of the granule question that led there.

The ideas doc specifies pathway activity in nucleus, cytoplasm, and RNA granules. The granule component is the
weak link.

Evidence (measured 2026-09-13; details in `data_processing.md`):
- With mcDETECT at the current setting, only 22% (CC) to 60% (LC) of tumor cells contain any stress granule, and
  13% to 58% have ≥10 genes detected in granules; the rest get zero-filled granule scores. A relaxed setting raises
  coverage to 52–70% but with smaller granules.
- Per-cell counts are sparse in every compartment: only 3% (CC) to 36% (LC) of panel genes are detected in ≥10% of
  tumor cells in nucleus or cytoplasm.
- Gene-wise, the cytoplasmic profile correlates with the granule profile (median Pearson 0.14–0.36; 0.26–0.40 on
  granule-positive cells) while nucleus vs granule is near zero. Paired DE found only about one granule-enriched
  gene per dataset.

Candidate directions (named G1-G3 here; the target forms A / B / C of the Tracker are a different, later
classification):
- **G1. Keep mcDETECT granule pathway scores** (current). Cons: low coverage, zeros dominate, few granule-enriched
  genes, sensitive to detection parameters.
- **G2. Drop granule construction; use cytoplasmic localization directly.** Outcome = cytoplasmic (vs nuclear)
  expression or enrichment of selected genes or gene sets: per-cell in-cytoplasm ratio of SG-marker genes,
  cytoplasmic pathway scores, or cytoplasm-minus-nucleus pathway contrasts. Pros: defined for every tumor cell,
  carries most of the granule-associated signal, no detection hyperparameters. Cons: loses the explicit aggregate
  phenotype; cytoplasmic counts are still sparse per cell.
- **G3. Hybrid.** A continuous per-cell "granule propensity" (cytoplasmic SG-marker aggregate, or a relaxed
  mcDETECT setting) as one head, plus nucleus and cytoplasm pathway heads.

Criteria: fraction of tumor cells with a non-degenerate value; signal-to-noise (variance not driven by zeros);
biological interpretability; correlation with stress scores; consistency across datasets; fit with the multi-head
regression design.

Analyses that addressed it: coverage and noise of a per-cell cytoplasmic outcome (Log 2026-09-22 and
`plans/fig1_exploration.md`: per-cell localization of single genes is mostly counting noise, so the outcome must
pool transcripts over genes or programs); cross-sample consistency of the cytoplasmic substates
(`plans/subtype_exploration.md`); correlation with stress (`plans/stress_substate_exploration.md`: granule-level
scores do not separate the substates, cytoplasmic expression is related to the environmental stress axes). Not
done: whether relaxed mcDETECT settings change the picture enough to keep G1 or G3.

### 4.4 Architecture, losses, evaluation

- **Encoders**: H&E patch → UNI → Embed1; own expression → scGPT → Embed2; neighbor expression → scGPT → Embed3
  (optional; this is the spatial-context embedding that gets swapped for virtual cells).
- **Fusion** into a unified representation (UR): attention-based mixing or simple concatenation; sparse attention
  for dimension reduction; must tolerate missing modalities (learned gating, modality dropout during training).
- **Heads**: multi-task regression heads for stress scores and for each compartment's pathway scores. Possibly
  per-task attention over modalities so each stress type weighs modalities differently.
- **Losses**: MSE on stress scores; MSE on pathway scores; latent alignment across modalities; single-modality
  auxiliary prediction losses; a term encouraging pan-cancer rather than cancer-type-specific responses (cancer type
  as confounder).
- **Evaluation**: hold out whole samples; stress and pathway prediction accuracy (pathways possibly binarized →
  AUROC); single vs multi modality; modality-attention visualization per stress type; correlation of each stress
  type with pathways per compartment to identify stress-specific responses.
- **Counterfactuals**: swap spatial-context embeddings and observe predicted stress and response; order cells by
  predicted stress and inspect response gradients.

Assets available: scGPT weights (`data/_utils/scGPT/`), UNI2-h weights (`data/_utils/UNI/`), and trial notebooks
in `code/3_model/` that extract scGPT cell embeddings on a merged-data subsample and UNI embeddings from
cell-centered H&E patches on OC.

## 5. Open questions

- Outcome variable (section 4.3; candidate recorded, decision pending in the Tracker).
- Uncertainty estimation for predictions: desired, method unspecified.
- Train/test split: by sample, but stratification across cancer types and handling of paired samples (Linghua
  patients) is undecided.
- Pathway-scoring metric per compartment: ssGSEA vs GSVA vs AUCell; may differ by compartment.
- Stress definitions still missing: monotonic (directional force at boundary), oscillatory (vessel-based); 0–1
  scaling of stress scores.
- Whether p53 expression tracks stress level.
- Gene-panel strategy once datasets with different panels are added.

## 6. Log

- **2026-09-13** — Documentation split: `CLAUDE.md` reduced to essentials; this plan document created. Outcome
  variable placed under review; no pipeline changes.
- **2026-09-15** — Document rescoped to background, data, and model design. Pipeline detail, measurements, and
  known code issues moved to `plans/data_processing.md`. Maintenance rule adopted: latest state only, stale content
  removed on update, decisions logged here.
- **2026-09-22** — Naive-model stage started (`plans/naive_model.md`, code in `naive_model/`). It probes outcome
  candidate **B** of §4.3 (cytoplasmic localization) with a refinement: the outcome is an offset binomial log odds
  ratio, `k` = cytoplasmic and `n` = cytoplasmic + nuclear counts over the 129 SG marker genes, with a per-cell
  offset equal to the logit of the cell's global non-SG cytoplasmic fraction. This cancels cell depth and the
  cell's global nucleus/cytoplasm balance algebraically, which matters because both naive alternatives are
  dominated by them: abundance is 51-88% explained by depth alone, and the raw cytoplasmic fraction is still 51%
  explained by depth. The genuine incremental signal from expression is ΔR² ≈ +0.065 (Prostate, held-out spatial
  half), not the ~0.58 a naive setup would report.
- **2026-09-22** — **Two findings that bear on §4.3 directly.**
  **(a) The per-cell SG-localization outcome is severely noise-limited.** With a median of 7-29 SG transcripts per
  tumor cell, 65-100% of the target's variance is binomial counting noise (reliability 0.00 in CC to 0.34 in
  Prostate), so the ceiling on any model's R² is ~0.05-0.34, not 1.0. Any future per-cell compartment outcome
  should be reported against that ceiling.
  **(b) Beyond the noise, the signal is genuinely absent in five of six samples.** A synthetic-signal test
  (inject a known linear log-OR of sd 0.5, keep each sample's real trial counts) recovers +0.027 in CC and +0.018
  in Prostate at their real depth -- so a moderate effect *would* be detected there. The real data gives +0.0003
  and -0.0016. Only LC shows real signal (+0.041, equivalent to a true log-OR sd of ~0.5), and LC is also the
  sample with the highest mcDETECT granule coverage (60.3% vs 21.9-48.2%, `data_processing.md` §4). So the one
  sample with the most stress-granule biology is the one where it is predictable; with n = 1 sample per cancer
  type we cannot say whether that is lung-specific or specimen-specific.
  **Implication for the outcome decision:** candidate B (cytoplasmic localization) is well-defined and
  confound-free in the offset-binomial form, but at Xenium 5K depth it is only informative in specimens with
  strong SG biology. Options are to aggregate (niche-level over cells, or pathway-level over genes -- which is
  what the existing ssGSEA scores already do), to weight samples by reliability, or to treat SG localization as a
  sample-level rather than cell-level phenotype. See `plans/naive_model.md` §4.
- **2026-11-14** — Status update after the exploratory stage (Fig. 1, subtype and stress explorations): status note
  added at the top; §3 OC sample identified; §4.2 points to the redefined stress axes; §4.3 records the current
  outcome candidate and renames the granule directions G1-G3 to avoid a clash with the Tracker's target forms; the
  outcome variable stays open. The naive model was dropped from the project on 2026-09-26: `plans/naive_model.md`
  and `naive_model/` cited in the two 2026-09-22 entries no longer exist (retained lessons: Tracker, "Evidence from
  the naive model"), and "candidate B" in those entries is direction G2 of §4.3.
