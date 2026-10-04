# code/1_exploration

Embedding and clustering of tumor cells on nuclear vs cytoplasmic expression, per Xenium 5K sample.

## Execution order

| Step | Where | Command | Output |
|---|---|---|---|
| 1. Embedding | HGCC | `mkdir -p logs && sbatch embedding.sh` from this folder | `output/1_exploration/<sample>/adata_<compartment>_embedded_<n_neighbors>_neighbors_<n_pcs>_pcs.h5ad` |
| 2. Pull results | local | `make out-pull` from the repo root | same files, local |
| 3. Clustering | local, `AIVC-env` | `clustering.ipynb` | clusters, ARI, UMAP and spatial plots |

## Files

- `embedding.py`: one work unit = (sample, compartment, n_pcs); normalize, log1p, PCA(100), t-SNE on `n_pcs`
  components, then kNN + UMAP for each `n_neighbors`. Grid: 6 samples x {nuclear, cytoplasmic} x n_pcs {20, 50}
  x n_neighbors {50, 100}. `python3 embedding.py --list` prints the 24 task ids; `--sample/--compartment/--n-pcs`
  runs one unit; without a selector it runs all units sequentially. Existing outputs are skipped unless `--force`.
- `embedding.sh`: SLURM array (`--array=0-23`) over the 24 units, 16 CPUs and 64 GB each; logs in `logs/`
  (gitignored). Partial submissions: `sbatch --array=0-3 embedding.sh` (BC only), `--array=0-23%8` (throttle).
- `clustering.ipynb`: Louvain clustering on the embedded objects, nuclear vs cytoplasmic ARI, plots. It currently
  reads the un-suffixed names `adata_<compartment>_embedded.h5ad` (an earlier single-parameter run); point it at one
  of the parameterised files above.
