# Data Processing and Outcome Definition

Last updated: 2026-09-15.

## 1. Purpose and scope

- This document tracks the data processing and outcome-definition stage: what the pipeline does today, what it has
  produced, measurements that inform decisions, known issues, and planned changes.
- It reflects the latest state only. On update, rewrite affected sections and delete stale content; record locked
  decisions and notable updates in the **Log** (section 7).
- Background, data inventory, and model design live in `plans/AIVC_plan.md`.

## 2. Current pipeline

`data/_utils/` (run once):
- `shared_genes.ipynb` → `shared_genes.npy/.csv`: intersection of the six public panels (5,001 genes).
- `cell_type_dict.ipynb` → `cell_type_dict.pkl` (10x label → merged label) and `cell_type_dict_manual.pkl`
  (graph-cluster id → merged label for LC, Prostate, Skin).
- `hallmark_pathways.R` → `all_pathways_filtered.gmt`: Hallmark + keyword-filtered Reactome and GO-BP sets
  (stress, translation, RNA processing, proteostasis keywords), restricted to the panel, ≥50 overlapping genes,
  Jaccard-deduplicated at 0.5 (the message text says 0.7). Also writes `filtered_pathways_{hallmark,reactome,gobp}.csv`.
- `hallmark_pathways_test.R` → `test_filtered_pathways_*_<dataset>.csv`: same idea with panel = top-100
  granule-detected genes per dataset and lower overlap thresholds. Exploratory; outputs are untracked.
- `SG_markers.xlsx`: stress-granule marker genes with "Fraction of RNA molecules in SGs"; threshold 0.4 downstream.

`code/1_preprocessing/`:
1. `1_clean.ipynb` — per dataset: filter transcripts to the panel; derive `in_cell`, `in_nucleus`,
   `overlaps_nucleus` (= `in_nucleus == in_cell`, so 1 for nuclear and for unassigned transcripts); write
   `processed_data/transcripts.parquet`; gene-wise in-cytoplasm ratio → `output/<ds>/in_cytoplasm_ratio.csv`; build
   AnnData from `cells.parquet` + `cell_feature_matrix.h5`; H&E pixel coordinates via the inverse of
   `HE_alignment.csv` and the pixel size in `experiment.xenium`; cell typing; write `intermediate_data/adata.h5ad`.
2. `2_merge.ipynb` — concatenate all datasets with `label="batch"`, shifting `global_x/global_y` per dataset so
   slides tile without overlap; tumor-only in-cytoplasm ratio → `output/merged_data/in_cytoplasm_ratio_tumor.csv`;
   write `merged_data/adata_all_raw.h5ad`.
3. `3_embedding.py` (SLURM) — normalize, log1p, PCA(100), t-SNE, UMAP → `merged_data/adata_embedded.h5ad`.
   **Not present on disk** (verified 2026-09-22); either never run to completion or the output was removed.
4. `4_plot.py` (SLURM) — UMAP and t-SNE colored by batch, merged cell type, EPCAM, KRT20.

`code/2_outcome/`:
1. `1_split_expression.ipynb` — malignant cells only: nuclear (`in_nucleus == 1`) and cytoplasmic
   (`overlaps_nucleus == 0`) count matrices → `processed_data/{nuclear,cytoplasmic}_expression_matrix.npz`, with
   `cell_ids.npy` / `gene_ids.npy` defining row and column order. Downstream steps assert these match the current
   tumor AnnData; re-run after any change to `1_clean`.
2. `2_ssGSEA.ipynb` — gseapy ssGSEA (CPM + log1p, rank normalization, 2,000-cell chunks, NES) on both matrices
   against `all_pathways_filtered.gmt` → `processed_data/ssgsea_hallmark_{nuclear,cytoplasmic}.parquet`;
   per-pathway Wilcoxon + BH; boxplots and spatial plots of top pathways.
3. `3_detection.py` (SLURM) — mcDETECT on tumor-cell transcripts. Detection collapses all SG markers (fraction
   > 0.4, in panel) into one pseudo-gene "Merged" (eps 1.5, minspl 3, size_thr 4, in_nucleus_thr 0.1–0.9);
   profiling then counts all panel genes per granule (buffer 0.05). Outputs `granules_{thr}_{eps}.csv`,
   `granule_adata_{thr}_{eps}.h5ad`, `output/<ds>/in_SG_ratio.csv`, granule maps.
4. `4_ssGSEA_SG.ipynb` — aggregate granule counts per cell → `SG_expression_matrix.npz`; ssGSEA on cells passing
   `min_total_counts` / `min_nnz_genes`, zeros elsewhere → `ssgsea_hallmark_sg.parquet`.
5. `5_score.ipynb` — on merged tumor cells: granule counts and the three score tables into `obs` (prefixes
   `nuclear_`, `cytoplasmic_`, `sg_`); stress scores and binary subtypes (definitions in `AIVC_plan.md` §4.2);
   correlations; write `merged_data/adata_tumor_scored.h5ad`.
- `plot_SG_markers.R` — in-cytoplasm ratio vs SG fraction; SG-marker vs non-marker boxplots with t-tests.

Conventions: figure size derived from slide extent (`scale = 5 / short_edge`, 10 for merged); axes and spines
stripped; JPEG at 300–500 dpi into `output/<dataset>/`.

## 3. Outputs produced so far

Per public dataset (`data/<dataset>/processed_data/`): `transcripts.parquet`; nuclear, cytoplasmic, and SG
expression matrices; `cell_ids.npy`, `gene_ids.npy`; granule files at two settings (`_0.25_1` and `_0.4_1.5`);
`ssgsea_hallmark_{nuclear,cytoplasmic,sg}.parquet`.

Merged (`data/merged_data/`): `adata_all_raw.h5ad` and `adata_tumor_scored.h5ad` (`adata_embedded.h5ad` is
**absent**)
(672,964 tumor cells; 31 pathways × 3 compartments; hypoxia, heat-shock, mechanical, T-cell-attack, and
immune-proximity scores with binary subtypes).

## 4. Measurements

Granule coverage (measured 2026-09-13 from `processed_data/` files):

| Dataset | Tumor cells | Granules (thr 0.4, eps 1.5) | Tumor cells with ≥1 granule | Granules (thr 0.25, eps 1.0) | Coverage |
|---|---|---|---|---|---|
| BC | 102,180 | 39,241 | 29.2% | 177,704 | 69.6% |
| OC | 160,250 | 68,591 | 32.6% | 149,230 | 52.4% |
| CC | 221,355 | 58,097 | 21.9% | 202,203 | 52.9% |
| LC | 44,624 | 62,787 | 60.3% | 89,816 | 66.0% |
| Prostate | 95,429 | 38,998 | 29.5% | 120,630 | 56.8% |
| Skin | 49,126 | 47,658 | 48.2% | 84,517 | 57.3% |

From `code/test/SG_vs_nuclei.ipynb` (gitignored; thr 0.4 granules):
- Tumor cells with ≥10 genes detected in SG: 13% (CC) to 58% (LC).
- Panel genes detected in ≥10% of tumor cells in nucleus or cytoplasm: 3% (CC) to 36% (LC). Median per-gene
  non-zero fraction is about 1–4% in both compartments.
- Gene-wise correlation across cells: cytoplasm vs SG median Pearson 0.14–0.36 (0.26–0.40 on SG-positive cells);
  nucleus vs SG about 0.01–0.03; nucleus vs cytoplasm about 0.02–0.04.
- Paired DE (nucleus vs SG on SG-positive cells) found about one SG-enriched gene per dataset.

Interpretation: the mcDETECT granule signal is sparse and zero-dominated, while cytoplasmic expression exists for
every cell and carries most of the granule-associated signal. This feeds the outcome-variable decision in
`AIVC_plan.md` §4.3.

## 5. Known issues

- Granule filenames have drifted: `granule_adata.h5ad` (read by `4_ssGSEA_SG`), `granule_adata_0.25_1.h5ad`,
  `granule_adata_0.4_1.5.h5ad` (written by `3_detection`, read by `5_score`).
- `ssgsea_hallmark_*` filenames are legacy; scores come from the combined Hallmark + Reactome + GO-BP GMT
  (31 pathways after filtering and dedup).
- The per-dataset coordinate shift table is duplicated in `2_merge.ipynb` and `5_score.ipynb`.
- `4_ssGSEA_SG.ipynb` is mostly commented out; its outputs exist from an earlier run.
- `hallmark_pathways_test.R` outputs (`test_filtered_pathways_*.csv`) are untracked in git.
- H&E-to-Xenium alignment has been checked visually on OC only.
- Notebooks often loop over a single dataset while iterating; outputs may not be current for all six.
- `adata_embedded.h5ad` does not exist on disk, though §2 lists it as a product of `3_embedding.py`.
- Only `ssgsea_hallmark_*` parquets exist; no GO-BP or Reactome ssGSEA outputs were ever written, despite scores
  coming from the combined GMT.
- `cell_type_merged` carries only 9-16 of the 18 categories per dataset (BC 16, OC 15, CC 15, LC 15, Prostate 9,
  Skin 10). Any per-dataset code that derives a cell-type vocabulary from `.cat.categories` will misalign columns
  across datasets; use a hard-coded union (see `naive_model/code/config.py:CELL_TYPES_18`).
- `nucleus_area` is NaN for cells segmented without a nucleus (1.5-8% of tumor cells); such cells have a
  degenerate nuclear/cytoplasmic split.
- The compartment matrices are built from all decoded transcripts: `1_clean` drops the raw `qv` column, so the
  QV >= 20 filter that 10x applies to `cell_feature_matrix.h5` is not applied to `nuclear/cytoplasmic_expression_matrix`.
  Nuclear + cytoplasmic therefore does not equal the `X` counts. Decision 2026-09-26: keep the matrices as they are
  for the Fig. 1 exploration and record the caveat; a Q20 re-split from `raw_data/transcripts.parquet` is the fix
  if it is ever needed.
- All six public samples are FFPE (`raw_data/experiment.xenium` -> `preservation_method`); the OC sample is
  `Xenium_Prime_Ovarian_Cancer_FFPE`.

## 6. Planned changes

Pending the outcome-variable decision (`AIVC_plan.md` §4.3):
1. Run the coverage and stress-correlation analyses for each candidate outcome (per dataset and pooled).
2. Revise `code/2_outcome/` accordingly; if granule construction is dropped, demote `3_detection` and
   `4_ssGSEA_SG` and make compartment-level outcomes the primary product.
3. Clean up the known issues above once the outcome pipeline is settled.
4. Implement the missing stress definitions (monotonic, oscillatory) and 0–1 scaling if adopted.
5. Ingest the Linghua and NC-paper datasets; settle the gene-panel strategy across panels.

## 7. Log

- **2026-09-15** — Document created from pipeline content split out of `plans/AIVC_plan.md`. No pipeline changes.
- **2026-09-26** — Known issues extended: no QV filter in the compartment matrices (kept, caveat recorded); FFPE
  status verified for all six samples. `code/2_exploration/` (Fig. 1 exploration) reads the compartment matrices and
  `intermediate_data/adata.h5ad` directly; no pipeline changes.
- **2026-09-22** — Corrections from the naive-model stage: `adata_embedded.h5ad` is absent; no GO-BP/Reactome
  ssGSEA outputs exist; per-dataset `cell_type_merged` vocabularies are incomplete; `nucleus_area` is NaN for
  nucleus-free cells. No pipeline changes.
