"""Embed tumor cells of five Xenium 5K samples on nuclear and cytoplasmic expression of cytoplasm-enriched genes.

Gene sets come from 1_enrichment.ipynb (../../output/1_exploration/enrichment/):
    shared5  overlap_genes_all.npy            genes enriched in cytoplasm in all five samples (~500)
    shared3  overlap_genes_3plus_samples.npy  genes enriched in cytoplasm in at least three samples (~1100)
Skin melanoma is not part of the enrichment analysis and is excluded here as well.

One work unit = (sample, compartment, gene set): load the tumor cells, subset the compartment matrix to the gene
set, normalize, log1p, PCA(50), t-SNE and kNN/UMAP on the first n_pcs components, write one .h5ad. The 20 units run
in parallel as a SLURM array (2_embedding.sh). Per sample the outputs are
    adata_{nuclear,cytoplasmic}_{shared5,shared3}_embedded.h5ad

Parameters (fixed per gene set, no sweep):
    n_pcs        30 for shared5 (~500 genes), 50 for shared3 (~1100 genes). The number of informative PCs scales
                 with the gene space; for sparse Xenium counts the variance curve is flat beyond these values. The
                 cumulative variance explained by the chosen n_pcs is printed so the choice can be checked in the log.
    n_neighbors  30 for both. It depends on the cell number (45k-220k tumor cells per sample), not on the gene
                 number: 15 (scanpy default) gives a noisy graph at this scale, 50-100 over-smooths small subtypes.

Usage (from this folder):
    python3 2_embedding.py --list                      # task id -> (sample, compartment, gene set)
    python3 2_embedding.py --task 7                    # one array task (default: $SLURM_ARRAY_TASK_ID)
    python3 2_embedding.py --sample Xenium_5K_LC --compartment cytoplasmic --gene-set shared5
    python3 2_embedding.py                             # no selector, no SLURM: all 20 tasks sequentially

Existing outputs are skipped unless --force is given, so a resubmission only redoes missing units.
"""
import argparse
import itertools
import os
import time

import numpy as np
import scanpy as sc
from scipy import sparse

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

SAMPLES = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate"]
COMPARTMENTS = ["nuclear", "cytoplasmic"]
ENRICHMENT_PATH = "../../output/1_exploration/enrichment/"
GENE_SETS = {
    "shared5": {"file": "overlap_genes_all.npy", "n_pcs": 30, "n_neighbors": 30},
    "shared3": {"file": "overlap_genes_3plus_samples.npy", "n_pcs": 50, "n_neighbors": 30},
}
N_COMPS_PCA = 50

# sample-major: ids 0-3 BC, 4-7 OC, 8-11 CC, 12-15 LC, 16-19 Prostate
# within a sample: nuclear/shared5, nuclear/shared3, cytoplasmic/shared5, cytoplasmic/shared3
TASKS = list(itertools.product(SAMPLES, COMPARTMENTS, list(GENE_SETS)))

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
sc.settings.n_jobs = N_JOBS


def output_file(sample, label, gene_set):
    return f"../../output/1_exploration/{sample}/adata_{label}_{gene_set}_embedded.h5ad"


def load_gene_set(gene_set):
    genes = np.load(ENRICHMENT_PATH + GENE_SETS[gene_set]["file"], allow_pickle=True)
    genes = [str(g) for g in genes]
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

    # subset to the gene set, keeping the panel's gene order
    mask = adata.var_names.isin(genes)
    assert mask.sum() == len(genes), f"{len(genes) - mask.sum()} genes of the set are missing from the panel"
    adata = adata[:, mask].copy()
    n_zero = int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum())
    print(f"{tag} subset to {adata.shape[1]} genes; {n_zero} cells with zero counts (kept)", flush=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=N_COMPS_PCA, svd_solver="auto")
    return adata


def run_task(sample, label, gene_set, force=False):
    tag = f"[{sample} {label} {gene_set}]"
    out = output_file(sample, label, gene_set)
    if os.path.exists(out) and not force:
        print(f"{tag} output exists, skipping", flush=True)
        return
    os.makedirs(f"../../output/1_exploration/{sample}/", exist_ok=True)
    n_pcs = GENE_SETS[gene_set]["n_pcs"]
    n_neighbors = GENE_SETS[gene_set]["n_neighbors"]
    genes = load_gene_set(gene_set)

    t0 = time.time()
    adata_tumor = load_tumor(sample)
    print(f"{tag} loaded {adata_tumor.shape[0]} tumor cells  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    adata = build_compartment(adata_tumor, sample, label, genes, tag)
    var_cum = float(np.cumsum(adata.uns["pca"]["variance_ratio"])[n_pcs - 1])
    print(f"{tag} normalize/log1p/PCA({N_COMPS_PCA})  {time.time() - t0:.0f}s; "
          f"first {n_pcs} PCs explain {var_cum:.1%} of variance", flush=True)

    t0 = time.time()
    sc.tl.tsne(adata, n_pcs=n_pcs, n_jobs=N_JOBS)
    print(f"{tag} t-SNE  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    adata.uns["embedding_params"] = {
        "gene_set": gene_set,
        "gene_set_file": GENE_SETS[gene_set]["file"],
        "n_genes": int(adata.shape[1]),
        "n_comps_pca": N_COMPS_PCA,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
    }
    adata.write(out)
    print(f"{tag} neighbors/UMAP written  {time.time() - t0:.0f}s", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task", type=int, default=None,
                   help=f"task id 0-{len(TASKS) - 1} (default: $SLURM_ARRAY_TASK_ID if set)")
    p.add_argument("--sample", choices=SAMPLES)
    p.add_argument("--compartment", choices=COMPARTMENTS)
    p.add_argument("--gene-set", choices=list(GENE_SETS))
    p.add_argument("--force", action="store_true", help="recompute even if the output exists")
    p.add_argument("--list", action="store_true", help="print the task table and exit")
    return p.parse_args()


def main():
    args = parse_args()
    if args.list:
        for i, (s, c, g) in enumerate(TASKS):
            print(f"{i:3d}  {s:20s} {c:12s} {g}  (n_pcs={GENE_SETS[g]['n_pcs']}, "
                  f"n_neighbors={GENE_SETS[g]['n_neighbors']})")
        return

    explicit = [args.sample, args.compartment, args.gene_set]
    if any(v is not None for v in explicit):
        assert all(v is not None for v in explicit), "--sample, --compartment and --gene-set must be given together"
        tasks = [(args.sample, args.compartment, args.gene_set)]
    else:
        task = args.task if args.task is not None else os.environ.get("SLURM_ARRAY_TASK_ID")
        if task is None:
            tasks = TASKS
        else:
            task = int(task)
            assert 0 <= task < len(TASKS), f"task id must be in 0-{len(TASKS) - 1}"
            tasks = [TASKS[task]]

    print(f"n_jobs={N_JOBS}; {len(tasks)} task(s)", flush=True)
    for sample, label, gene_set in tasks:
        t0 = time.time()
        run_task(sample, label, gene_set, force=args.force)
        print(f"[{sample} {label} {gene_set}] done in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
