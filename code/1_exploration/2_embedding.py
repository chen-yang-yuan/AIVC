"""Embed tumor cells of five Xenium 5K samples on nuclear and cytoplasmic expression of cytoplasm-enriched genes.

Gene set: ../../output/1_exploration/enrichment/overlap_genes_3plus_samples.npy from 1_enrichment.ipynb, i.e. genes
that are highly variable in >= 2 samples and cytoplasm-enriched in >= 3 of the five samples (~600 genes). Skin
melanoma is not part of the enrichment analysis and is excluded here as well.

One work unit = (sample, compartment): load the tumor cells, normalize_total(1e4) on the FULL panel (size factor =
total transcripts of the compartment, stored in obs["compartment_total_counts"]), subset to the gene set, log1p
(kept in layers["lognorm"] for DE and plots), regress out log10 compartment depth, scale, PCA(50), t-SNE and
kNN/UMAP on the first N_PCS components, write one .h5ad. Normalization comes before the subset on purpose:
normalizing the subset would use the target-gene count as the size factor, inflating cells with few target-gene
transcripts and letting the abundant genes of the set dominate the denominator.

Depth correction: the compartment vectors are sparse (BC cytoplasm: median 130 transcripts per cell on the panel,
39 on the gene set), so normalize_total + log1p alone leaves PC1 = depth (Spearman 0.91 with log depth in BC) and
Louvain then returns depth strata instead of substates. Regressing out log10(compartment_total_counts + 1) and
z-scoring the genes (clipped at 10) removes this (Spearman -0.2) and yields interpretable clusters (proliferation,
interferon/complement, hypoxia/lysosomal, ...). Spearman(PC1, log depth) is printed to the log as a sanity check;
|rho| > 0.3 means the correction failed. All cells are kept (the notebooks assert the same cell set as adata_tumor);
cells with (almost) no target-gene transcripts form a low-depth cluster that the clustering notebook flags.

The 10 units run in parallel as a SLURM array (2_embedding.sh). Per sample the outputs are
    adata_{nuclear,cytoplasmic}_embedded.h5ad      X = regressed + scaled, layers["lognorm"] = log-normalized

Parameters (fixed, no sweep):
    N_PCS = 30         ~600 genes; 20 PCs / 15 neighbors give the same clusters in BC, so the graph is robust to
                       these choices. After scaling the first N_PCS components explain only ~7% of the variance,
                       which is expected for z-scored sparse data; the value is printed to the log.
    N_NEIGHBORS = 30   depends on the cell number (45k-220k tumor cells per sample): 15 (scanpy default) gives a
                       noisy graph at this scale, 50-100 over-smooths small subtypes.

Usage (from this folder):
    python3 2_embedding.py --list                      # task id -> (sample, compartment)
    python3 2_embedding.py --task 7                    # one array task (default: $SLURM_ARRAY_TASK_ID)
    python3 2_embedding.py --sample Xenium_5K_LC --compartment cytoplasmic
    python3 2_embedding.py                             # no selector, no SLURM: all 10 tasks sequentially

Existing outputs are skipped unless --force is given, so a resubmission only redoes missing units.
"""
import argparse
import itertools
import os
import time

import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

SAMPLES = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate"]
COMPARTMENTS = ["nuclear", "cytoplasmic"]
GENE_SET_FILE = "../../output/1_exploration/enrichment/overlap_genes_3plus_samples.npy"
N_COMPS_PCA = 50
N_PCS = 30
N_NEIGHBORS = 30
SCALE_MAX_VALUE = 10   # clip z-scores after scaling (scanpy convention)

# sample-major: ids 0-1 BC, 2-3 OC, 4-5 CC, 6-7 LC, 8-9 Prostate; within a sample: nuclear, cytoplasmic
TASKS = list(itertools.product(SAMPLES, COMPARTMENTS))

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
sc.settings.n_jobs = N_JOBS


def output_file(sample, label):
    return f"../../output/1_exploration/{sample}/adata_{label}_embedded.h5ad"


def load_gene_set():
    genes = [str(g) for g in np.load(GENE_SET_FILE, allow_pickle=True)]
    assert len(genes) == len(set(genes)), "gene set contains duplicates"
    return genes


def load_tumor(sample):
    intermediate_path = f"../../data/{sample}/intermediate_data/"
    processed_path = f"../../data/{sample}/processed_data/"
    adata_all = sc.read_h5ad(intermediate_path + "adata.h5ad")
    cell_ids = np.load(processed_path + "cell_ids.npy", allow_pickle=True)
    assert adata_all.shape[0] > cell_ids.shape[0], "cell_ids should be a subset of adata.obs_names"
    return adata_all[adata_all.obs["cell_id"].isin(cell_ids)].copy()


def build_compartment(adata_tumor, sample, label, genes, tag):
    processed_path = f"../../data/{sample}/processed_data/"
    compartment_expression = sparse.load_npz(processed_path + f"{label}_expression_matrix.npz")
    assert compartment_expression.shape[0] == adata_tumor.shape[0], \
        "compartment_expression should have the same number of cells as adata_tumor"
    adata = sc.AnnData(X=compartment_expression, obs=adata_tumor.obs, var=adata_tumor.var)

    # size factor = total transcripts of the compartment over the full panel, computed BEFORE the gene subset;
    # normalizing after subsetting would let the target-gene count act as the size factor (cells with few
    # target-gene transcripts get inflated, abundant set genes dominate the denominator)
    sc.pp.normalize_total(adata, target_sum=1e4, key_added="compartment_total_counts")

    # subset to the gene set, keeping the panel's gene order; everything downstream (log1p, PCA, graph) uses it
    mask = adata.var_names.isin(genes)
    assert mask.sum() == len(genes), f"{len(genes) - mask.sum()} genes of the set are missing from the panel"
    adata = adata[:, mask].copy()
    n_zero = int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum())  # zero iff raw target-gene counts are zero
    print(f"{tag} subset to {adata.shape[1]} genes; {n_zero} cells with zero target-gene counts (kept)", flush=True)

    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()  # log-normalized values for DE and expression plots

    # depth correction: regress out log10 compartment depth and z-score the genes, otherwise PC1 = depth and the
    # clusters are depth strata (see module docstring)
    adata.obs["log_compartment_total_counts"] = np.log10(adata.obs["compartment_total_counts"].values + 1)
    sc.pp.regress_out(adata, "log_compartment_total_counts", n_jobs=N_JOBS)
    sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)

    sc.tl.pca(adata, n_comps=N_COMPS_PCA, svd_solver="auto")
    rho = spearmanr(adata.obsm["X_pca"][:, 0], adata.obs["log_compartment_total_counts"].values).correlation
    print(f"{tag} Spearman(PC1, log10 compartment depth) = {rho:+.2f} (|rho| > 0.3 means depth still dominates)",
          flush=True)
    return adata


def run_task(sample, label, force=False):
    tag = f"[{sample} {label}]"
    out = output_file(sample, label)
    if os.path.exists(out) and not force:
        print(f"{tag} output exists, skipping", flush=True)
        return
    os.makedirs(f"../../output/1_exploration/{sample}/", exist_ok=True)
    genes = load_gene_set()

    t0 = time.time()
    adata_tumor = load_tumor(sample)
    print(f"{tag} loaded {adata_tumor.shape[0]} tumor cells  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    adata = build_compartment(adata_tumor, sample, label, genes, tag)
    var_cum = float(np.cumsum(adata.uns["pca"]["variance_ratio"])[N_PCS - 1])
    print(f"{tag} normalize/log1p/regress_out/scale/PCA({N_COMPS_PCA})  {time.time() - t0:.0f}s; "
          f"first {N_PCS} PCs explain {var_cum:.1%} of variance", flush=True)

    t0 = time.time()
    sc.tl.tsne(adata, n_pcs=N_PCS, n_jobs=N_JOBS)
    print(f"{tag} t-SNE  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    sc.tl.umap(adata)
    adata.uns["embedding_params"] = {
        "gene_set_file": os.path.basename(GENE_SET_FILE),
        "n_genes": int(adata.shape[1]),
        "normalization": "normalize_total(1e4) on full panel before gene subset, then log1p (layers['lognorm'])",
        "depth_correction": f"regress_out(log10(compartment_total_counts + 1)) + scale(max_value={SCALE_MAX_VALUE})",
        "n_comps_pca": N_COMPS_PCA,
        "n_pcs": N_PCS,
        "n_neighbors": N_NEIGHBORS,
    }
    adata.write(out)
    print(f"{tag} neighbors/UMAP written  {time.time() - t0:.0f}s", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task", type=int, default=None,
                   help=f"task id 0-{len(TASKS) - 1} (default: $SLURM_ARRAY_TASK_ID if set)")
    p.add_argument("--sample", choices=SAMPLES)
    p.add_argument("--compartment", choices=COMPARTMENTS)
    p.add_argument("--force", action="store_true", help="recompute even if the output exists")
    p.add_argument("--list", action="store_true", help="print the task table and exit")
    return p.parse_args()


def main():
    args = parse_args()
    if args.list:
        for i, (s, c) in enumerate(TASKS):
            print(f"{i:3d}  {s:20s} {c}")
        return

    explicit = [args.sample, args.compartment]
    if any(v is not None for v in explicit):
        assert all(v is not None for v in explicit), "--sample and --compartment must be given together"
        tasks = [(args.sample, args.compartment)]
    else:
        task = args.task if args.task is not None else os.environ.get("SLURM_ARRAY_TASK_ID")
        if task is None:
            tasks = TASKS
        else:
            task = int(task)
            assert 0 <= task < len(TASKS), f"task id must be in 0-{len(TASKS) - 1}"
            tasks = [TASKS[task]]

    print(f"n_jobs={N_JOBS}; {len(tasks)} task(s); n_pcs={N_PCS}, n_neighbors={N_NEIGHBORS}", flush=True)
    for sample, label in tasks:
        t0 = time.time()
        run_task(sample, label, force=args.force)
        print(f"[{sample} {label}] done in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
