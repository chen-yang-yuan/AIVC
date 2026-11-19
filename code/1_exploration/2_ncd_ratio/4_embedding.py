"""NCD ratio analysis, step 4: pooled embedding and clustering (plans/ncd_ratio_exploration.md, section 5).

Embeds and clusters the tumor cells of the five samples on the ratio score r = log2((c + 1) / (n + 1)) of a target
gene list, pooled after a correction within each sample. Heavy: runs on HGCC as a SLURM array (4_embedding.sh), one
task per combination of

  list        primary (V-abs first 500), larger (V-abs first 1000), union (V-abs 500 + V-null 500, 630 genes)
  correction  zscore, depth                          (ncd_ratio_utils.correct_within_sample)
  run         data, null, half_a, half_b

= 3 x 2 x 4 = 24 tasks (task id = index in TASKS, list slowest; `python 4_embedding.py --status` prints them).

Procedure of one task
1. Per sample: cells of the analysis and protein-coding genes (ncd_ratio_utils.load_counts); counts of the genes of
   the list. Cell covariates: depth = all transcripts of the cell over the protein-coding genes (c + n), nuclear
   share of these transcripts, coordinates and areas from intermediate_data/adata.h5ad. The covariates are always
   those of the observed cell, also in the null and half runs.
2. Counts of the run:
     data    observed
     null    one draw of the count-preserving null: every (cell, gene) pair keeps its total, the transcripts are
             reassigned with the cytoplasmic fraction of the gene in the sample (ncd_ratio_utils.null_counts)
     half_a, half_b   every transcript of a pair goes to half A with probability 0.5, otherwise to half B (the two
             tasks use the same random split, so the halves sum to the data)
3. r for every pair; empty pairs are 0. All cells of the analysis are kept.
4. Correction within the sample, per gene (ncd_ratio_utils.correct_within_sample):
     zscore  z-score of r over the cells of the sample
     depth   residual of the OLS of r on log10(depth + 1), then z-score
   Clipped at +- SCALE_MAX. Nothing else is removed at cell level; the cell-wide lean to one compartment stays in.
5. Pooled: the five corrected matrices stacked; PCA (N_COMPS_PCA, arpack); neighbours (N_NEIGHBORS on N_PCS PCs);
   UMAP (not for the halves); clustering on the neighbour graph at every resolution of RES_GRID (Louvain as in the
   cytoplasmic analysis; Leiden if the louvain package is not installed, recorded in summary.txt).
6. Checks: Spearman of PC1 to PC10 with the cell covariates (pooled and per sample); sample effect (share of the
   variance of each PC between samples; share of a cell's neighbours from its own sample against the sample's share
   of the cells); the correction itself. Before clipping every gene has mean 0 and SD 1 in every sample (and slope 0
   on depth in the depth version) by construction; the clipping changes this for rarely detected genes, whose few
   non-zero scores lie far out, so the check reports the share of clipped scores and the mean, SD and slope after
   clipping.

Output: output/1_exploration/ncd_ratio/4_embedding/<list>/<correction>/<run>/   (--quick: 4_embedding_quick/)
  cells.csv.gz              per cell: sample, cell_id, coordinates, areas, depth, nuclear share, detected genes of
                            the list, UMAP, cluster at each resolution (cluster_r<resolution>)
  pca.npy, pca_variance.csv cells x N_COMPS_PCA (float32, order of cells.csv.gz); variance ratio per PC
  gene_scaling.csv          per sample and gene: mean and SD of r, slope on depth, detection rate
  cluster_gene_means.csv.gz per resolution, cluster and gene: mean corrected score, mean r, detection rate
  cluster_sizes.csv         per resolution and cluster: cells per sample
  checks.csv                the correlations and the sample effect
  connectivities.npz        the neighbour graph (data runs only)
  summary.txt               written last; its presence marks the task as complete
No count matrix is written: the cells are those of load_counts in the same order (cells.csv.gz has the ids).

Usage (from this folder):
    python 4_embedding.py --task 0                       # one task by id (the SLURM array uses this)
    python 4_embedding.py --list primary --correction depth --run data
    python 4_embedding.py --all --quick [--quick-cells 4000]   # development: all tasks on a subsample, local
    python 4_embedding.py --status                       # which tasks are complete
    add --force to recompute a complete task
"""
import argparse
import os
import shutil
import time
import warnings

warnings.filterwarnings("ignore")

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

import ncd_ratio_utils as nr

sc.settings.verbosity = 0
N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
sc.settings.n_jobs = N_JOBS

TASKS = [(l, c, r) for l in nr.LISTS for c in nr.CORRECTIONS for r in nr.RUNS]
STEP = "4_embedding"
COVARIATES = ["log_depth", "nuclear_share", "nucleus_area_share", "cell_area", "detected_genes"]
N_PC_CHECK = 10
OBS_COLS = ["global_x", "global_y", "cell_area", "nucleus_area"]


def cluster_method():
    try:
        import louvain  # noqa: F401
        return "louvain"
    except Exception:
        import leidenalg  # noqa: F401
        return "leiden"


def task_dir(task, quick):
    return nr.step_dir(STEP, quick) + "/".join(task) + "/"


def is_done(task, quick):
    return os.path.exists(task_dir(task, quick) + "summary.txt")


def sample_part(j, sample, genes_list, correction, run, quick_cells):
    """Corrected score matrix of one sample (dense, cells x genes of the list), the sparse r, the cell table and
    the per-gene scaling table."""
    Xc, Xn, cells, genes = nr.load_counts(sample, quick_cells, np.random.default_rng([nr.SEED, j, 0]))
    tc, tn = np.asarray(Xc.sum(axis=1)).ravel(), np.asarray(Xn.sum(axis=1)).ravel()
    pos = pd.Series(np.arange(len(genes)), index=genes).reindex(genes_list).values
    assert not np.isnan(pos).any(), "a gene of the list is not among the genes of the analysis"
    pos = pos.astype(int)
    rows, cols, c, n = nr.pairs(Xc[:, pos].tocsr(), Xn[:, pos].tocsr())
    del Xc, Xn
    G = len(genes_list)
    if run == "null":
        p = np.bincount(cols, weights=c, minlength=G) / np.maximum(np.bincount(cols, weights=c + n, minlength=G), 1)
        c0, n0 = nr.null_counts(c, n, p[cols], np.random.default_rng([nr.SEED, j, 1]))
        assert (c0 + n0 == c + n).all()
        c, n = c0, n0
    elif run in ("half_a", "half_b"):
        rng = np.random.default_rng([nr.SEED, j, 2])   # the same split in both half tasks
        ca, na = rng.binomial(c, 0.5), rng.binomial(n, 0.5)
        c, n = (ca, na) if run == "half_a" else (c - ca, n - na)
        keep = (c + n) > 0
        rows, cols, c, n = rows[keep], cols[keep], c[keep], n[keep]
    r = nr.ratio_score(c, n).astype(np.float32)
    R = sparse.csr_matrix((r, (rows, cols)), shape=(len(cells), G))
    D = sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(cells), G))   # detected pairs

    ad = anndata.read_h5ad(f"../../../data/{sample}/intermediate_data/adata.h5ad", backed="r")
    obs = ad.obs[["cell_id"] + OBS_COLS].copy()
    obs["cell_id"] = obs["cell_id"].astype(str)
    obs = obs.set_index("cell_id").reindex(cells)
    assert not obs["cell_area"].isna().any(), "a cell of the analysis is missing from the AnnData"
    table = pd.DataFrame({"sample": nr.short(sample), "cell_id": cells, **{k: obs[k].values for k in OBS_COLS},
                          "depth": (tc + tn).astype(int), "log_depth": np.log10(tc + tn + 1.0),
                          "nuclear_share": tn / np.maximum(tc + tn, 1),
                          "nucleus_area_share": (obs["nucleus_area"] / obs["cell_area"]).values,
                          "detected_genes": np.asarray(D.sum(axis=1)).ravel().astype(int)})
    X, sc_ = nr.correct_within_sample(R.toarray(), table["log_depth"].values, correction)
    scaling = pd.DataFrame({"sample": nr.short(sample), "gene": genes_list, **sc_,
                            "detection_rate": np.asarray(D.sum(axis=0)).ravel() / len(cells),
                            "mean_corrected": X.mean(axis=0), "sd_corrected_final": X.std(axis=0, ddof=1)})
    d = table["log_depth"].values - table["log_depth"].values.mean()
    scaling["slope_after"] = (d @ X.astype(np.float64)) / (d @ d)
    return X, R, D, table, scaling


def spearman_table(pcs, cov):
    """Spearman of every column of pcs with every column of cov (rows with a missing covariate dropped)."""
    out = {}
    rp = pd.DataFrame(pcs).rank().values
    for k in cov.columns:
        v = cov[k].values.astype(float)
        ok = np.isfinite(v)
        rv = pd.Series(v[ok]).rank().values
        a = rp[ok] if ok.all() else pd.DataFrame(pcs[ok]).rank().values
        a = a - a.mean(axis=0); rv = rv - rv.mean()
        out[k] = (a * rv[:, None]).sum(axis=0) / np.sqrt((a ** 2).sum(axis=0) * (rv ** 2).sum())
    return pd.DataFrame(out)


def run_task(task, quick_cells, force):
    lst, correction, run = task
    quick = quick_cells is not None
    out = task_dir(task, quick)
    if is_done(task, quick) and not force:
        print(f"[{'/'.join(task)}] complete, skipped (--force to recompute)", flush=True)
        return
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    t0 = time.time()
    method = cluster_method()
    genes_list = nr.load_target_genes(lst)
    tag = f"[{'/'.join(task)}]"
    print(f"{tag} {len(genes_list)} genes; n_jobs {N_JOBS}; clustering {method}{'; QUICK ' + str(quick_cells) if quick else ''}", flush=True)

    Xs, Rs, Ds, tables, scalings = [], [], [], [], []
    for j, sample in enumerate(nr.SAMPLES):
        X, R, D, table, scaling = sample_part(j, sample, genes_list, correction, run, quick_cells)
        Xs.append(X); Rs.append(R); Ds.append(D); tables.append(table); scalings.append(scaling)
        print(f"{tag} {nr.short(sample)}: {len(table)} cells, {R.nnz / R.shape[0] / R.shape[1]:.3f} of scores non-zero ({time.time() - t0:.0f} s)", flush=True)
    cells = pd.concat(tables, ignore_index=True)
    scaling = pd.concat(scalings, ignore_index=True)
    scaling.to_csv(out + "gene_scaling.csv", index=False)
    X = np.vstack(Xs); R = sparse.vstack(Rs).tocsr(); D = sparse.vstack(Ds).tocsr()
    del Xs, Rs, Ds
    n_cells, G = X.shape
    assert n_cells == len(cells) and (cells["sample"] + ":" + cells["cell_id"]).is_unique

    # the correction did what it should (genes without variance in a sample are 0 and excluded)
    live = scaling["sd_corrected"] > 1e-12
    clipped = float((np.abs(X) >= nr.SCALE_MAX).mean())
    check_corr = {"genes_without_variance_in_a_sample": int((~live).sum()), "share_of_scores_clipped": clipped,
                  "max_abs_mean_after_clipping": float(scaling.loc[live, "mean_corrected"].abs().max()),
                  "sd_after_clipping_min": float(scaling.loc[live, "sd_corrected_final"].min()),
                  "sd_after_clipping_max": float(scaling.loc[live, "sd_corrected_final"].max()),
                  "max_abs_slope_on_depth_after_clipping": float(scaling.loc[live, "slope_after"].abs().max()),
                  "median_abs_slope_on_depth_before": float(scaling["depth_slope"].abs().median())}

    adata = anndata.AnnData(X=X, obs=pd.DataFrame(index=(cells["sample"] + ":" + cells["cell_id"]).values), var=pd.DataFrame(index=genes_list))
    sc.tl.pca(adata, n_comps=min(nr.N_COMPS_PCA, G - 1), svd_solver="arpack")
    pcs = adata.obsm["X_pca"]
    np.save(out + "pca.npy", pcs.astype(np.float32))
    vr = adata.uns["pca"]["variance_ratio"]
    pd.DataFrame({"pc": np.arange(1, len(vr) + 1), "variance_ratio": vr, "cumulative": np.cumsum(vr)}).to_csv(out + "pca_variance.csv", index=False)
    print(f"{tag} PCA: first {nr.N_PCS} PCs explain {vr[:nr.N_PCS].sum():.1%} ({time.time() - t0:.0f} s)", flush=True)

    sc.pp.neighbors(adata, n_neighbors=nr.N_NEIGHBORS, n_pcs=nr.N_PCS, random_state=nr.SEED)
    print(f"{tag} neighbours ({time.time() - t0:.0f} s)", flush=True)
    if run in ("data", "null"):
        sc.tl.umap(adata, random_state=nr.SEED)
        cells["umap_1"], cells["umap_2"] = adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1]
        print(f"{tag} UMAP ({time.time() - t0:.0f} s)", flush=True)
    for res in nr.RES_GRID:
        key = f"cluster_r{res}"
        (sc.tl.louvain if method == "louvain" else sc.tl.leiden)(adata, resolution=res, key_added=key, random_state=nr.SEED)
        cells[key] = adata.obs[key].astype(int).values
        print(f"{tag} resolution {res}: {cells[key].nunique()} clusters ({time.time() - t0:.0f} s)", flush=True)
    cells.to_csv(out + "cells.csv.gz", index=False)
    if run == "data":
        sparse.save_npz(out + "connectivities.npz", adata.obsp["connectivities"].tocsr().astype(np.float32))

    # ---- checks
    npc = min(N_PC_CHECK, pcs.shape[1])
    rows = []
    rho = spearman_table(pcs[:, :npc], cells[COVARIATES])
    for k in COVARIATES:
        for i in range(npc):
            rows.append({"check": "spearman_pc_covariate", "sample": "pooled", "covariate": k, "pc": i + 1, "value": rho[k].iloc[i]})
    codes = cells["sample"].astype("category").cat.codes.values
    shorts = list(cells["sample"].astype("category").cat.categories)
    for s_i, s in enumerate(shorts):
        m = codes == s_i
        rs = spearman_table(pcs[m, :npc], cells.loc[m, COVARIATES])
        for k in COVARIATES:
            for i in range(npc):
                rows.append({"check": "spearman_pc_covariate", "sample": s, "covariate": k, "pc": i + 1, "value": rs[k].iloc[i]})
    onehot = sparse.csr_matrix((np.ones(n_cells), (np.arange(n_cells), codes)), shape=(n_cells, len(shorts)))
    share = np.asarray(onehot.sum(axis=0)).ravel() / n_cells
    grand = pcs.mean(axis=0)
    means = (onehot.T @ pcs) / (share[:, None] * n_cells)
    eta = (share[:, None] * n_cells * (means - grand) ** 2).sum(axis=0) / ((pcs - grand) ** 2).sum(axis=0)
    for i in range(pcs.shape[1]):
        rows.append({"check": "variance_between_samples", "sample": "pooled", "covariate": "sample", "pc": i + 1, "value": eta[i]})
    eta_w = float((eta[:nr.N_PCS] * vr[:nr.N_PCS]).sum() / vr[:nr.N_PCS].sum())
    knn = (adata.obsp["distances"] > 0).astype(np.float32)      # the neighbours of every cell
    same = np.asarray((knn @ onehot).multiply(onehot).sum(axis=1)).ravel() / np.maximum(np.asarray(knn.sum(axis=1)).ravel(), 1)
    for s_i, s in enumerate(shorts):
        rows.append({"check": "same_sample_neighbour_share", "sample": s, "covariate": "sample", "pc": 0, "value": same[codes == s_i].mean()})
        rows.append({"check": "sample_share_of_cells", "sample": s, "covariate": "sample", "pc": 0, "value": share[s_i]})
    for k, v in check_corr.items():
        rows.append({"check": k, "sample": "pooled", "covariate": "", "pc": 0, "value": v})
    checks = pd.DataFrame(rows)
    checks.to_csv(out + "checks.csv", index=False)

    # ---- clusters: sizes by sample, mean scores by gene
    size_rows, mean_parts = [], []
    for res in nr.RES_GRID:
        lab = cells[f"cluster_r{res}"].values
        K = lab.max() + 1
        M = sparse.csr_matrix((np.ones(n_cells), (lab, np.arange(n_cells))), shape=(K, n_cells))
        sizes = np.asarray(M.sum(axis=1)).ravel()
        by_sample = np.asarray((M @ onehot).todense())
        for k in range(K):
            size_rows.append({"resolution": res, "cluster": k, "cells": int(sizes[k]), **{s: int(by_sample[k, i]) for i, s in enumerate(shorts)},
                              **{f"median_{c_}": float(np.nanmedian(cells[c_].values[lab == k])) for c_ in COVARIATES}})
        mz, mr, md = (M @ X) / sizes[:, None], np.asarray((M @ R).todense()) / sizes[:, None], np.asarray((M @ D).todense()) / sizes[:, None]
        mean_parts.append(pd.DataFrame({"resolution": res, "cluster": np.repeat(np.arange(K), G), "gene": np.tile(genes_list, K),
                                        "mean_corrected": mz.ravel().astype(np.float32), "mean_r": mr.ravel().astype(np.float32),
                                        "detection_rate": md.ravel().astype(np.float32)}))
    sizes_table = pd.DataFrame(size_rows)
    sizes_table.to_csv(out + "cluster_sizes.csv", index=False)
    pd.concat(mean_parts, ignore_index=True).to_csv(out + "cluster_gene_means.csv.gz", index=False, float_format="%.5g")

    # ---- summary (written last: marks the task as complete)
    pd.set_option("display.width", 250)
    pooled = checks[(checks["check"] == "spearman_pc_covariate") & (checks["sample"] == "pooled")].pivot(index="covariate", columns="pc", values="value").loc[COVARIATES]
    nn = checks[checks["check"].isin(["same_sample_neighbour_share", "sample_share_of_cells"])].pivot(index="sample", columns="check", values="value").loc[shorts]
    lines = [f"NCD ratio analysis, step 4: pooled embedding and clustering{' (QUICK RUN)' if quick else ''}", "",
             f"Task: list {lst} ({G} genes), correction {correction}, run {run}",
             f"Cells: {n_cells} (" + ", ".join(f"{s} {int(round(share[i] * n_cells))}" for i, s in enumerate(shorts)) + ")",
             f"Parameters: PCA {pcs.shape[1]} components (arpack), neighbours {nr.N_NEIGHBORS} on {nr.N_PCS} PCs, clustering {method} at {nr.RES_GRID}, clip +-{nr.SCALE_MAX:g}, seed {nr.SEED}", "",
             "Correction check", "\n".join(f"  {k}: {v:.4g}" for k, v in check_corr.items()), "",
             f"PCA: variance explained by PC1 to PC10: " + " ".join(f"{v:.3f}" for v in vr[:10]) + f"; first {nr.N_PCS} PCs: {vr[:nr.N_PCS].sum():.3f}", "",
             "Spearman of the PCs with the cell covariates (pooled)", pooled.to_string(float_format=lambda v: f"{v:+.2f}"), "",
             "Sample effect",
             "  share of the variance between samples, PC1 to PC10: " + " ".join(f"{v:.3f}" for v in eta[:10]),
             f"  the same, weighted over the first {nr.N_PCS} PCs: {eta_w:.4f}",
             nn.to_string(float_format=lambda v: f"{v:.3f}"), "",
             "Clusters by resolution (cells per sample)", sizes_table[["resolution", "cluster", "cells"] + shorts].to_string(index=False), "",
             f"done in {(time.time() - t0) / 60:.1f} min"]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", type=int, help=f"task id, 0 to {len(TASKS) - 1}")
    ap.add_argument("--list", choices=nr.LISTS)
    ap.add_argument("--correction", choices=nr.CORRECTIONS)
    ap.add_argument("--run", choices=nr.RUNS)
    ap.add_argument("--all", action="store_true", help="all tasks one after the other (development, with --quick)")
    ap.add_argument("--quick", action="store_true", help="subsample of the cells, output in 4_embedding_quick/")
    ap.add_argument("--quick-cells", type=int, default=nr.QUICK_CELLS, help="cells per sample of a quick run")
    ap.add_argument("--force", action="store_true", help="recompute a complete task")
    ap.add_argument("--status", action="store_true", help="list the tasks and whether they are complete")
    args = ap.parse_args()
    if args.status:
        done = [is_done(t, args.quick) for t in TASKS]
        for i, (t, d) in enumerate(zip(TASKS, done)):
            print(f"{i:2d}  {'/'.join(t):28s} {'complete' if d else 'MISSING'}")
        print(f"{sum(done)} of {len(TASKS)} tasks complete" + ("" if all(done) else "; rerun the missing ones with: sbatch --array=" + ",".join(str(i) for i, d in enumerate(done) if not d) + " 4_embedding.sh"))
        return
    if args.all:
        tasks = TASKS
    elif args.task is not None:
        tasks = [TASKS[args.task]]
    else:
        assert args.list and args.correction and args.run, "give --task, --all, or --list, --correction and --run"
        tasks = [(args.list, args.correction, args.run)]
    for t in tasks:
        run_task(t, args.quick_cells if args.quick else None, args.force)


if __name__ == "__main__":
    main()
