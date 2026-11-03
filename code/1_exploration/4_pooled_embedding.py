"""Embed the tumor cells of the five Xenium 5K samples POOLED on the cytoplasmic expression of the gene set.

Purpose: test whether the subtypes found per sample (3_clustering.ipynb) are genuine cross-sample substates or
per-sample artifacts. The pooled object is processed with exactly the per-sample procedure of 2_embedding.py, only
on all cells at once; 4_pooled.ipynb then clusters and annotates it with the per-sample code (subtype_utils.py) and
compares the pooled labels with the per-sample labels cell by cell.

Procedure (cytoplasmic compartment only; functions and parameters imported from 2_embedding.py)
1. Per sample: load the tumor cells and the cytoplasmic compartment matrix, normalize_total(1e4) on the FULL panel
   (size factor = total cytoplasmic transcripts, obs["compartment_total_counts"]), subset to the gene set
   (overlap_genes_3plus_samples.npy), log1p -> layers["lognorm"]. obs keeps sample, cell_id, coordinates, areas, cell
   type and depth; obs_names = "<sample>:<cell_id>" so that cells stay unique after pooling.
2. Concatenate the five objects (same gene order asserted), then on the POOLED matrix: regress out
   log10(compartment_total_counts + 1), scale (z-score per gene, clipped at 10), PCA(50), t-SNE and kNN / UMAP on the
   first N_PCS components. The regression intercept is the pooled gene mean, so differences between the sample means
   are deliberately left in: nothing is regressed on the sample, because the batch effect between samples is what
   the pooled analysis examines (3rd question of 4_pooled.ipynb).
3. Write ../../output/1_exploration/pooled/adata_pooled_cytoplasmic_embedded.h5ad: X = regressed + scaled,
   layers["lognorm"], obsm X_pca / X_tsne / X_umap, the kNN graph in obsp, uns["embedding_params"].

Runs as one SLURM job (4_pooled_embedding.sh; ~624k cells x 634 genes, dense float32 = 1.6 GB, 96 GB requested).
Spearman(PC1-5, log depth) is printed as the depth-correction check (|rho| > 0.3 means the correction failed), and the
per-sample cell numbers must add up to the pooled number.

Usage (from this folder):
    python3 4_pooled_embedding.py                # all five samples, with t-SNE
    python3 4_pooled_embedding.py --no-tsne      # skip t-SNE (the slowest step)
    python3 4_pooled_embedding.py --force        # recompute even if the output exists
"""
import argparse
import importlib
import os
import time

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

emb = importlib.import_module("2_embedding")   # per-sample embedding: same functions and parameters

SAMPLES = emb.SAMPLES
OUTPUT_DIR = "../../output/1_exploration/pooled/"
OUTPUT_FILE = OUTPUT_DIR + "adata_pooled_cytoplasmic_embedded.h5ad"
OBS_COLS = ["cell_id", "global_x", "global_y", "cell_area", "nucleus_area", "cell_type_merged"]

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
sc.settings.n_jobs = N_JOBS


def load_sample_lognorm(sample, genes):
    """Tumor cells of one sample on the cytoplasmic compartment: normalize_total on the full panel, gene subset, log1p.
    Returns an AnnData with X = log-normalized (sparse) and the metadata columns needed downstream."""
    tag = f"[{sample}]"
    t0 = time.time()
    adata_tumor = emb.load_tumor(sample)
    processed_path = f"../../data/{sample}/processed_data/"
    compartment_expression = sparse.load_npz(processed_path + "cytoplasmic_expression_matrix.npz")
    assert compartment_expression.shape[0] == adata_tumor.shape[0], \
        "compartment_expression should have the same number of cells as adata_tumor"
    adata = sc.AnnData(X=compartment_expression, obs=adata_tumor.obs[OBS_COLS].copy(), var=adata_tumor.var[[]].copy())
    del adata_tumor

    # size factor = total cytoplasmic transcripts over the full panel, computed BEFORE the gene subset (2_embedding.py)
    sc.pp.normalize_total(adata, target_sum=1e4, key_added="compartment_total_counts")
    mask = adata.var_names.isin(genes)
    assert mask.sum() == len(genes), f"{len(genes) - mask.sum()} genes of the set are missing from the panel"
    adata = adata[:, mask].copy()
    n_zero = int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum())
    sc.pp.log1p(adata)

    adata.obs["sample"] = sample
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs_names = sample + ":" + adata.obs["cell_id"].values
    assert adata.obs_names.is_unique, f"{tag} duplicated cell ids"
    print(f"{tag} {adata.shape[0]} tumor cells, {adata.shape[1]} genes; {n_zero} cells with zero target-gene counts (kept); "
          f"median cytoplasmic depth {np.median(adata.obs['compartment_total_counts']):.0f}  {time.time() - t0:.0f}s", flush=True)
    return adata


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--no-tsne", action="store_true", help="skip t-SNE")
    p.add_argument("--force", action="store_true", help="recompute even if the output exists")
    args = p.parse_args()

    if os.path.exists(OUTPUT_FILE) and not args.force:
        print(f"{OUTPUT_FILE} exists, skipping (use --force to recompute)", flush=True)
        return
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"n_jobs={N_JOBS}; samples={SAMPLES}; n_pcs={emb.N_PCS}, n_neighbors={emb.N_NEIGHBORS}", flush=True)
    genes = emb.load_gene_set()

    # 1. per-sample log-normalized objects, pooled into one
    parts = [load_sample_lognorm(sample, genes) for sample in SAMPLES]
    for part in parts[1:]:
        assert (part.var_names == parts[0].var_names).all(), "gene order differs between samples"
    adata = ad.concat(parts, join="inner", merge="same", index_unique=None)
    adata.var = parts[0].var.copy()
    n_per_sample = {sample: int(part.shape[0]) for sample, part in zip(SAMPLES, parts)}
    del parts
    assert adata.shape[0] == sum(n_per_sample.values())
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    adata.layers["lognorm"] = adata.X.copy()
    print(f"[pooled] {adata.shape[0]} cells x {adata.shape[1]} genes: {n_per_sample}", flush=True)

    # 2. depth correction and embedding on the pooled matrix (2_embedding.py procedure)
    t0 = time.time()
    adata.obs["log_compartment_total_counts"] = np.log10(adata.obs["compartment_total_counts"].values + 1)
    adata.X = emb.regress_out_covariate(adata.X, adata.obs["log_compartment_total_counts"].values)
    sc.pp.scale(adata, max_value=emb.SCALE_MAX_VALUE)
    sc.tl.pca(adata, n_comps=emb.N_COMPS_PCA, svd_solver="arpack")
    for k in range(5):
        rho = spearmanr(adata.obsm["X_pca"][:, k], adata.obs["log_compartment_total_counts"].values).correlation
        print(f"[pooled] Spearman(PC{k + 1}, log10 cytoplasmic depth) = {rho:+.2f} (|rho| > 0.3 means depth still dominates)", flush=True)
    var_cum = float(np.cumsum(adata.uns["pca"]["variance_ratio"])[emb.N_PCS - 1])
    print(f"[pooled] regress_out/scale/PCA({emb.N_COMPS_PCA})  {time.time() - t0:.0f}s; first {emb.N_PCS} PCs explain {var_cum:.1%} of variance", flush=True)

    if not args.no_tsne:
        t0 = time.time()
        sc.tl.tsne(adata, n_pcs=emb.N_PCS, n_jobs=N_JOBS)
        print(f"[pooled] t-SNE  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    sc.pp.neighbors(adata, n_neighbors=emb.N_NEIGHBORS, n_pcs=emb.N_PCS)
    sc.tl.umap(adata)
    print(f"[pooled] neighbors/UMAP  {time.time() - t0:.0f}s", flush=True)

    # 3. write
    adata.uns["embedding_params"] = {
        "gene_set_file": os.path.basename(emb.GENE_SET_FILE),
        "n_genes": int(adata.shape[1]),
        "samples": list(SAMPLES),
        "n_cells_per_sample": n_per_sample,
        "normalization": "per sample: normalize_total(1e4) on full panel before gene subset, then log1p (layers['lognorm'])",
        "depth_correction": f"pooled OLS residuals on log10(compartment_total_counts + 1) + scale(max_value={emb.SCALE_MAX_VALUE}); "
                            "no correction for sample",
        "n_comps_pca": emb.N_COMPS_PCA,
        "n_pcs": emb.N_PCS,
        "n_neighbors": emb.N_NEIGHBORS,
        "tsne": not args.no_tsne,
    }
    adata.write(OUTPUT_FILE)
    print(f"[pooled] written {OUTPUT_FILE}", flush=True)


if __name__ == "__main__":
    t_start = time.time()
    main()
    print(f"done in {(time.time() - t_start) / 60:.1f} min", flush=True)
