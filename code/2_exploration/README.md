# Fig. 1 exploration: nuclear vs cytoplasmic expression in tumor cells

Exploratory analysis of the six processed Xenium 5K FFPE samples that motivates AIVC's compartment-level target and
selects the candidate target gene set. Design, status and decisions: `plans/fig1_exploration.md`. Specification:
`plans/AIVC Model Design Tracker.md`, section "Next step: exploratory analysis for Fig. 1".

## Execution order

Run everything from this folder (`code/2_exploration`). Local work uses the `AIVC-env` conda env; HGCC jobs use
`preprocessing-env` (same package versions for this code).

| Order | What | Where | Command | Writes |
|---|---|---|---|---|
| 1 | Development check (optional but recommended): light + heavy steps on a 20k-cell subsample of one dataset | local | `python 1_heavy.py --dataset Xenium_5K_LC --quick`, then run the notebook with `QUICK = True`, `RUN_DATASETS = ["Xenium_5K_LC"]` | `output/2_exploration/Xenium_5K_LC_quick/` |
| 2 | Heavy steps on all kept tumor cells, all six datasets (clustering, cytoplasm-beyond-nucleus residual analysis, DE, association, pathway association) | HGCC | `make push` (local) → `git pull` (HGCC) → `cd ~/hulab/projects/AIVC/code/2_exploration && mkdir -p logs && sbatch 1_heavy.sh` | `output/2_exploration/<ds>/` on HGCC |
| 3 | Retrieve the heavy results | local | `make out-pull` | `output/2_exploration/<ds>/` locally |
| 4 | Main notebook: light steps (0, 1, 2, 5, 8-light), loads the heavy results, all figures, per-dataset candidates, cross-sample section, decision summary | local (Jupyter) | open `fig1_exploration.ipynb`, `QUICK = False`, run all | `output/2_exploration/<ds>/`, `.../<ds>/fig/`, `output/2_exploration/cross_sample/` |
| 5 | Record the results and decisions | — | update `plans/fig1_exploration.md` (status + Log) and the tracker | — |

The processed compartment matrices (`data/<ds>/processed_data/{nuclear,cytoplasmic}_expression_matrix.npz`,
`cell_ids.npy`, `gene_ids.npy`) and `intermediate_data/adata.h5ad` must exist on HGCC (`make data-push-Xenium`).
A single dataset on HGCC: `sbatch --array=3 1_heavy.sh` (index into `config.DATASETS`: 0 BC, 1 OC, 2 CC, 3 LC,
4 Prostate, 5 Skin). A subset of steps: `python 1_heavy.py --dataset Xenium_5K_CC --steps 7,8` (earlier steps are
loaded from cache); step names are 3, 3b, 4, 7, 8.

Every step caches its tables and is skipped when its files exist; use the `FORCE` dict in the notebook or
`--force` on the driver to recompute. Expected runtimes: light steps 1–4 min per dataset locally (step 8-light is the
slowest, ~1–5 min); heavy steps on HGCC roughly 20 min (LC, Skin) to 1.5 h (CC) per dataset; the quick local check
takes ~5 min.

## Files

| File | Role |
|---|---|
| `config.py` | Datasets, paths, thresholds, cell-type union, positive controls, marker sets, quick-mode overrides |
| `annotation_qc.py` | Step 0b (report only): marker-set scores for all cells, class profiles, tumor-set purity flags, profiles of the manually annotated malignant clusters, cross-dataset transfer of the non-epithelial labels |
| `io_utils.py` | Loaders with order assertions, Dropbox hydration, GMT reader, pathway indicator matrix, cache helpers |
| `compartment.py` | Depths, leave-one-gene-out offsets, per-gene log odds ratio and reliability, gene classes, rarefaction, half-splits, tiles, split-half concordance, tile bootstrap |
| `clustering.py` | Identical clustering pipeline per matrix (igraph Leiden; `cluster_from_pcs` for any embedding), ARI/NMI/Jaccard, cyto-only clusters, spatial purity, cluster feature enrichment |
| `residual.py` | Step 3b: nuclear design (reliability-weighted PCs, cluster one-hot, depth spline, covariates), errors-in-variables residualisation, residual PCA, leakage floor, strata, stratified purity/effects, Moran's I, dose-response, the pre-registered eligibility rule |
| `de.py` | Pseudobulk paired log-FC with tile bootstrap (global and within-stratum), Wilcoxon top-k overlap, cyto-only / nuc-only calls, centroids |
| `neighborhood.py` | Kernel geometry over all cells, composition, crowding, boundary/core, niches, quantile bins |
| `association.py` | Vectorised categorical GLMs (Poisson closed form, offset-binomial Newton), deviance explained, permutation baseline, cross-classification increments |
| `controls.py` | Abundance-bin percentiles, matched control genes, matched random gene sets |
| `plotting.py` | Figure helpers |
| `pipeline.py` | The step functions (`prepare`, `step0_qc` … `step8_pathways_assoc`, `run_heavy`, `load_heavy`) shared by the notebook and the driver |
| `cross_sample.py` | Per-dataset candidate calls (step 9), recurrence, between-sample agreement, cluster recurrence, decision summary |
| `1_heavy.py`, `1_heavy.sh` | HGCC driver for steps 3, 4, 7, 8 (SLURM array over the six datasets) |
| `fig1_exploration.ipynb` | Main notebook (built by `_build_notebook.py`; edit the notebook directly afterwards) |

## Output inventory (`output/2_exploration/<ds>/`)

Light steps: `qc_summary.json`, `cells.parquet` (per-cell depths, filters, morphology bins, tile id),
`annotation_provenance.json`, `annotation_qc_flags.json`, `annotation_qc_class_scores.parquet`, `annotation_qc_cells.parquet`,
`annotation_qc_tumor_clusters.parquet` (manual datasets), `annotation_transfer.parquet` (step 0b; the classifier and its
leave-one-dataset-out reference are cached in `cross_sample/`),
`genes_compartment.parquet` (A1 + B5 percentiles + strata), `genes_concordance_{cell,tile}.parquet` (A2),
`controls_sets.parquet` (SG vs matched control), `pathways_compartment.parquet`, `pathways_concordance.parquet` (C8),
`candidates.parquet` (C9).

Heavy steps: `nbhd_features.parquet`, `clusters.parquet` (labels for the 8 matrices × 2 resolutions), `pcs.npz`,
`clustering_pairwise.parquet`, `jaccard_*.csv`, `cyto_only_clusters.parquet`, `cluster_recovery_total_by_nuc.parquet`,
`cluster_spatial.parquet`, `cluster_features.parquet` (A3); `resid_diagnostics.json` (verdict, panel rule), `resid_moran.parquet`,
`resid_cluster_summary.parquet`, `resid_clusters.parquet`, `resid_scores.parquet`, `resid_loadings.npz`,
`cluster_spatial_stratified.parquet`, `de_paired_residpart.parquet` (A3b); `de_paired_{nucpart,cytopart}.parquet`, `de_topk_*.json`,
`de_topk_overlap.parquet`, `centroids_cyto_m.npz` (A4); `assoc_genes_{long,wide}.parquet` (C7);
`pathways_assoc_long.parquet` (C8); `cells_heavy.parquet` (per-cell labels and neighbourhood features).

Figures in `fig/` (`A0_*`, `A1_*`, `A2_*`, `A3_*`, `A3b_*`, `A4_*`, `C7_*`, `C8_*`); cross-sample tables and figures in
`output/2_exploration/cross_sample/` (`gene_recurrence.parquet`, `candidates_final.csv`, `beta_by_sample.csv`,
`decision_summary.csv`, `fig/B6_*`).

## Conventions

- Tumor cells only (`cell_type_merged == "Malignant cell"`). The label comes from 10x supervised labels in BC/OC/CC and
  from manually annotated graph clusters in LC/Prostate/Skin; step 0b reports how consistent the tumor sets are but
  does not change the definition. Nuclear = `in_nucleus == 1`, cytoplasmic = in cell and
  not in nucleus; depths are sums over the 5001 shared panel genes (never obs `total_counts`).
- Matched-depth analyses (clustering, DE, cell-level concordance) use cells with a single nucleus and
  `min(nuclear, cytoplasmic depth) >= D_MIN`; per-gene localisation statistics use all tumor cells.
- The compartment matrices contain all decoded transcripts (no QV >= 20 filter); recorded as a caveat.
