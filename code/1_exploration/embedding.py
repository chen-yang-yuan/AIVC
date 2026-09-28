"""Embed tumor cells of the six Xenium 5K samples on nuclear and cytoplasmic expression.

One work unit = (sample, compartment, n_pcs): load the tumor cells, normalize, log1p, PCA(100), t-SNE on the first
n_pcs components, then for each n_neighbors build the kNN graph and UMAP and write one .h5ad. The units are run in
parallel as a SLURM array (embedding.sh); t-SNE is the slow step and depends only on n_pcs, so it is computed once
per unit and shared by the n_neighbors outputs.

Usage (from this folder):
    python3 embedding.py --list                       # task id -> (sample, compartment, n_pcs)
    python3 embedding.py --task 7                     # one array task (default: $SLURM_ARRAY_TASK_ID)
    python3 embedding.py --sample Xenium_5K_LC --compartment nuclear --n-pcs 20 --n-neighbors 50
    python3 embedding.py                              # no selector, no SLURM: all 24 tasks sequentially

Outputs: ../../output/1_exploration/<sample>/adata_<compartment>_embedded_<n_neighbors>_neighbors_<n_pcs>_pcs.h5ad
Existing outputs are skipped unless --force is given, so a resubmission only redoes missing combinations.
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

SAMPLES = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate", "Xenium_5K_Skin"]
COMPARTMENTS = ["nuclear", "cytoplasmic"]
N_PCS = [20, 50]
N_NEIGHBORS = [50, 100]
N_COMPS_PCA = 100

# sample-major: ids 0-3 are BC, 4-7 OC, 8-11 CC, 12-15 LC, 16-19 Prostate, 20-23 Skin
TASKS = list(itertools.product(SAMPLES, COMPARTMENTS, N_PCS))

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
sc.settings.n_jobs = N_JOBS


def output_file(sample, label, n_neighbors, n_pcs):
    return f"../../output/1_exploration/{sample}/adata_{label}_embedded_{n_neighbors}_neighbors_{n_pcs}_pcs.h5ad"


def load_tumor(sample):
    intermediate_path = f"../../data/{sample}/intermediate_data/"
    processed_path = f"../../data/{sample}/processed_data/"
    adata_all = sc.read_h5ad(intermediate_path + "adata.h5ad")
    cell_ids = np.load(processed_path + "cell_ids.npy", allow_pickle=True)
    assert adata_all.shape[0] > cell_ids.shape[0], "cell_ids should be a subset of adata.obs_names"
    return adata_all[adata_all.obs["cell_id"].isin(cell_ids)].copy()


def build_compartment(adata_tumor, sample, label):
    processed_path = f"../../data/{sample}/processed_data/"
    compartment_expression = sparse.load_npz(processed_path + f"{label}_expression_matrix.npz")
    assert compartment_expression.shape[0] == adata_tumor.shape[0], \
        "compartment_expression should have the same number of cells as adata_tumor"
    adata = sc.AnnData(X=compartment_expression, obs=adata_tumor.obs, var=adata_tumor.var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=N_COMPS_PCA, svd_solver="auto")
    return adata


def run_task(sample, label, n_pcs, n_neighbors_list, force=False):
    tag = f"[{sample} {label} n_pcs={n_pcs}]"
    todo = [k for k in n_neighbors_list if force or not os.path.exists(output_file(sample, label, k, n_pcs))]
    for k in n_neighbors_list:
        if k not in todo:
            print(f"{tag} n_neighbors={k}: output exists, skipping", flush=True)
    if not todo:
        return
    os.makedirs(f"../../output/1_exploration/{sample}/", exist_ok=True)

    t0 = time.time()
    adata_tumor = load_tumor(sample)
    print(f"{tag} loaded {adata_tumor.shape[0]} tumor cells  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    adata = build_compartment(adata_tumor, sample, label)
    print(f"{tag} normalize/log1p/PCA({N_COMPS_PCA})  {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    sc.tl.tsne(adata, n_pcs=n_pcs, n_jobs=N_JOBS)
    print(f"{tag} t-SNE  {time.time() - t0:.0f}s", flush=True)

    for n_neighbors in todo:
        t0 = time.time()
        adata_tmp = adata.copy()
        sc.pp.neighbors(adata_tmp, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata_tmp)
        adata_tmp.write(output_file(sample, label, n_neighbors, n_pcs))
        print(f"{tag} n_neighbors={n_neighbors}: neighbors/UMAP written  {time.time() - t0:.0f}s", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task", type=int, default=None,
                   help=f"task id 0-{len(TASKS) - 1} (default: $SLURM_ARRAY_TASK_ID if set)")
    p.add_argument("--sample", choices=SAMPLES)
    p.add_argument("--compartment", choices=COMPARTMENTS)
    p.add_argument("--n-pcs", type=int, choices=N_PCS)
    p.add_argument("--n-neighbors", type=int, nargs="+", default=N_NEIGHBORS, choices=N_NEIGHBORS)
    p.add_argument("--force", action="store_true", help="recompute even if the output exists")
    p.add_argument("--list", action="store_true", help="print the task table and exit")
    return p.parse_args()


def main():
    args = parse_args()
    if args.list:
        for i, (s, c, n) in enumerate(TASKS):
            print(f"{i:3d}  {s:20s} {c:12s} n_pcs={n}")
        return

    explicit = [args.sample, args.compartment, args.n_pcs]
    if any(v is not None for v in explicit):
        assert all(v is not None for v in explicit), "--sample, --compartment and --n-pcs must be given together"
        tasks = [(args.sample, args.compartment, args.n_pcs)]
    else:
        task = args.task if args.task is not None else os.environ.get("SLURM_ARRAY_TASK_ID")
        if task is None:
            tasks = TASKS
        else:
            task = int(task)
            assert 0 <= task < len(TASKS), f"task id must be in 0-{len(TASKS) - 1}"
            tasks = [TASKS[task]]

    print(f"n_jobs={N_JOBS}; {len(tasks)} task(s); n_neighbors={args.n_neighbors}", flush=True)
    for sample, label, n_pcs in tasks:
        t0 = time.time()
        run_task(sample, label, n_pcs, args.n_neighbors, force=args.force)
        print(f"[{sample} {label} n_pcs={n_pcs}] done in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
