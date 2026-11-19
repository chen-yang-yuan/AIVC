"""Embed the tumor cells of the five Xenium 5K samples POOLED on the cytoplasmic (or nuclear) expression of the gene set.

Purpose: test whether the subtypes found per sample (3_clustering.ipynb) are genuine cross-sample substates or
per-sample artifacts. The pooled object is processed with exactly the per-sample procedure of 2_embedding.py, only
on all cells at once; 4_pooled.ipynb then clusters and annotates it with the per-sample code (subtype_utils.py) and
compares the pooled labels with the per-sample labels cell by cell.

Compartment (--compartment, default cytoplasmic): the same procedure runs on the NUCLEAR compartment matrix of the
same 634 genes (--compartment nuclear, output folder with the suffix _nuclear, object
adata_pooled_nuclear_embedded.h5ad) to test whether the nuclear expression of the cytoplasm-enriched genes carries
the same substates; 4_pooled_sample_corrected_nuclear.ipynb compares the nuclear pooled labels with the cytoplasmic
sample-corrected pooled labels cell by cell. Only the matrix file and the size factor (total transcripts of that
compartment) change.

Procedure (one compartment; functions and parameters imported from 2_embedding.py)
1. Per sample: load the tumor cells and the compartment matrix, normalize_total(1e4) on the FULL panel
   (size factor = total transcripts of the compartment, obs["compartment_total_counts"]), subset to the gene set
   (overlap_genes_3plus_samples.npy), log1p -> layers["lognorm"]. obs keeps sample, cell_id, coordinates, areas, cell
   type and depth; obs_names = "<sample>:<cell_id>" so that cells stay unique after pooling.
2. Depth correction, in one of two modes (--correction):
   - none (default; output folder pooled/): concatenate the five objects (same gene order asserted), then on the
     POOLED matrix regress out log10(compartment_total_counts + 1) and scale (z-score per gene, clipped at 10). The
     regression intercept is the pooled gene mean and the scaling centres on it, so the differences between the
     sample means are deliberately left in: nothing is corrected for sample, because the batch effect between
     samples is what this run examines (3rd question of 4_pooled.ipynb). Result 2026-11-03: the samples separate
     completely, driven by 10-20 lineage genes (H19, XBP1, NKX3-1, SORD, EPDR1, ...) whose sample means differ by
     1-2 z-units.
   - sample (output folder pooled_sample_corrected/): the depth regression and the scaling run WITHIN each sample,
     exactly the two calls of build_compartment() in 2_embedding.py, and the five transformed matrices are stacked.
     Every gene then has mean 0 and unit variance within every sample, so sample-level mean differences and the
     sample-specific depth slopes are removed before pooling (the pooled X equals the X of the five per-sample
     adata_cytoplasmic_embedded.h5ad objects). Differences between samples in average program activity are
     removed with them, so pooled label frequencies are relative within each sample; differences in the
     within-sample covariance structure are not corrected. The uncorrected run stays the reference for the raw
     sample effect.
   Then PCA(50; arpack solver here, scanpy's default solver in 2_embedding.py), t-SNE and kNN / UMAP on the first N_PCS components,
   identical in both modes.
3. Write <output folder>/adata_pooled_<compartment>_embedded.h5ad: X = regressed + scaled, layers["lognorm"],
   obsm X_pca / X_tsne / X_umap, the kNN graph in obsp, uns["embedding_params"] (records mode and compartment).

Runs as one SLURM job (4_pooled_embedding.sh; ~624k cells x 634 genes, dense float32 = 1.6 GB, 96 GB requested).
Spearman(PC1-5, log depth) is printed as the depth-correction check (|rho| > 0.3 means the correction failed), and the
per-sample cell numbers must add up to the pooled number.

Usage (from this folder):
    python3 4_pooled_embedding.py                      # uncorrected pool, with t-SNE -> pooled/
    python3 4_pooled_embedding.py --correction sample  # per-sample transform, pooled PCA -> pooled_sample_corrected/
    python3 4_pooled_embedding.py --correction sample --compartment nuclear   # nuclear counterpart -> pooled_sample_corrected_nuclear/
    python3 4_pooled_embedding.py --no-tsne            # skip t-SNE (the slowest step)
    python3 4_pooled_embedding.py --force              # recompute even if the output exists
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
COMPARTMENTS = emb.COMPARTMENTS   # ["nuclear", "cytoplasmic"]
OUTPUT_DIRS = {"none": "../../../output/1_exploration/pooled/",
               "sample": "../../../output/1_exploration/pooled_sample_corrected/"}


def output_dir(correction, compartment):
    """pooled/ or pooled_sample_corrected/ for the cytoplasmic compartment, with the suffix _nuclear for the nuclear one."""
    return OUTPUT_DIRS[correction].rstrip("/") + ("" if compartment == "cytoplasmic" else f"_{compartment}") + "/"


def output_name(compartment):
    return f"adata_pooled_{compartment}_embedded.h5ad"

OBS_COLS = ["cell_id", "global_x", "global_y", "cell_area", "nucleus_area", "cell_type_merged"]

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
sc.settings.n_jobs = N_JOBS


def load_sample_lognorm(sample, genes, correction="none", compartment="cytoplasmic"):
    """Tumor cells of one sample on one compartment: normalize_total on the full panel, gene subset, log1p
    (kept in layers["lognorm"]). With correction="sample", X is additionally depth-regressed and scaled within the
    sample (build_compartment() of 2_embedding.py); otherwise X stays the sparse log-normalized matrix.
    Returns an AnnData with the metadata columns needed downstream."""
    tag = f"[{sample} {compartment}]"
    t0 = time.time()
    adata_tumor = emb.load_tumor(sample)
    processed_path = f"../../../data/{sample}/processed_data/"
    compartment_expression = sparse.load_npz(processed_path + f"{compartment}_expression_matrix.npz")
    assert compartment_expression.shape[0] == adata_tumor.shape[0], \
        "compartment_expression should have the same number of cells as adata_tumor"
    adata = sc.AnnData(X=compartment_expression, obs=adata_tumor.obs[OBS_COLS].copy(), var=adata_tumor.var[[]].copy())
    del adata_tumor

    # size factor = total transcripts of the compartment over the full panel, computed BEFORE the gene subset (2_embedding.py)
    sc.pp.normalize_total(adata, target_sum=1e4, key_added="compartment_total_counts")
    mask = adata.var_names.isin(genes)
    assert mask.sum() == len(genes), f"{len(genes) - mask.sum()} genes of the set are missing from the panel"
    adata = adata[:, mask].copy()
    n_zero = int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum())
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    adata.obs["log_compartment_total_counts"] = np.log10(adata.obs["compartment_total_counts"].values + 1)

    if correction == "sample":
        # per-sample depth correction and scaling, the same two calls as build_compartment() in 2_embedding.py:
        # removes this sample's gene means and depth slopes, so nothing sample-specific at the mean level is pooled
        adata.X = emb.regress_out_covariate(adata.X, adata.obs["log_compartment_total_counts"].values)
        sc.pp.scale(adata, max_value=emb.SCALE_MAX_VALUE)

    adata.obs["sample"] = sample
    adata.obs["cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs_names = sample + ":" + adata.obs["cell_id"].values
    assert adata.obs_names.is_unique, f"{tag} duplicated cell ids"
    print(f"{tag} {adata.shape[0]} tumor cells, {adata.shape[1]} genes; {n_zero} cells with zero target-gene counts (kept); "
          f"median {compartment} depth {np.median(adata.obs['compartment_total_counts']):.0f}"
          f"{'; depth regressed and scaled within the sample' if correction == 'sample' else ''}  {time.time() - t0:.0f}s",
          flush=True)
    return adata


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--correction", choices=["none", "sample"], default="none",
                   help="none: depth regression and scaling on the pooled matrix (no sample correction); "
                        "sample: depth regression and scaling within each sample, then pooled PCA")
    p.add_argument("--compartment", choices=COMPARTMENTS, default="cytoplasmic",
                   help="compartment matrix to embed (default cytoplasmic); nuclear writes to <folder>_nuclear/")
    p.add_argument("--no-tsne", action="store_true", help="skip t-SNE")
    p.add_argument("--force", action="store_true", help="recompute even if the output exists")
    args = p.parse_args()

    out_dir = output_dir(args.correction, args.compartment)
    output_file = out_dir + output_name(args.compartment)
    if os.path.exists(output_file) and not args.force:
        print(f"{output_file} exists, skipping (use --force to recompute)", flush=True)
        return
    os.makedirs(out_dir, exist_ok=True)
    print(f"correction={args.correction}; compartment={args.compartment}; n_jobs={N_JOBS}; samples={SAMPLES}; n_pcs={emb.N_PCS}, "
          f"n_neighbors={emb.N_NEIGHBORS}", flush=True)
    genes = emb.load_gene_set()

    # 1. per-sample objects (log-normalized; with --correction sample also depth-regressed and scaled), pooled into one
    parts = [load_sample_lognorm(sample, genes, args.correction, args.compartment) for sample in SAMPLES]
    for part in parts[1:]:
        assert (part.var_names == parts[0].var_names).all(), "gene order differs between samples"
    adata = ad.concat(parts, join="inner", merge="same", index_unique=None)
    adata.var = parts[0].var.copy()
    n_per_sample = {sample: int(part.shape[0]) for sample, part in zip(SAMPLES, parts)}
    del parts
    assert adata.shape[0] == sum(n_per_sample.values())
    assert "lognorm" in adata.layers
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    print(f"[pooled] {adata.shape[0]} cells x {adata.shape[1]} genes: {n_per_sample}", flush=True)

    # 2. depth correction (pooled mode only; already done per sample otherwise) and embedding (2_embedding.py procedure)
    t0 = time.time()
    if args.correction == "none":
        adata.X = emb.regress_out_covariate(adata.X, adata.obs["log_compartment_total_counts"].values)
        sc.pp.scale(adata, max_value=emb.SCALE_MAX_VALUE)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)
    sc.tl.pca(adata, n_comps=emb.N_COMPS_PCA, svd_solver="arpack")
    for k in range(5):
        rho = spearmanr(adata.obsm["X_pca"][:, k], adata.obs["log_compartment_total_counts"].values).correlation
        print(f"[pooled] Spearman(PC{k + 1}, log10 {args.compartment} depth) = {rho:+.2f} (|rho| > 0.3 means depth still dominates)", flush=True)
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
        "compartment": args.compartment,
        "gene_set_file": os.path.basename(emb.GENE_SET_FILE),
        "n_genes": int(adata.shape[1]),
        "samples": list(SAMPLES),
        "n_cells_per_sample": n_per_sample,
        "normalization": "per sample: normalize_total(1e4) on full panel before gene subset, then log1p (layers['lognorm'])",
        "depth_correction": (f"{'pooled' if args.correction == 'none' else 'per-sample'} OLS residuals on "
                             f"log10(compartment_total_counts + 1) + scale(max_value={emb.SCALE_MAX_VALUE})"),
        "sample_correction": ("none: sample means left in" if args.correction == "none" else
                              "per-sample depth regression and scaling before pooling (gene means and depth slopes "
                              "removed within each sample), then pooled PCA / kNN / UMAP / t-SNE"),
        "n_comps_pca": emb.N_COMPS_PCA,
        "n_pcs": emb.N_PCS,
        "n_neighbors": emb.N_NEIGHBORS,
        "tsne": not args.no_tsne,
    }
    adata.write(output_file)
    print(f"[pooled] written {output_file}", flush=True)


if __name__ == "__main__":
    t_start = time.time()
    main()
    print(f"done in {(time.time() - t_start) / 60:.1f} min", flush=True)
