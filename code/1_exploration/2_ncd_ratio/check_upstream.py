"""NCD ratio analysis: audit of the compartment count matrices the analysis reads (plans/ncd_ratio_exploration.md,
section 10.8). Not a step of the analysis; it changes no result of steps 0 to 3.

The matrices data/<sample>/processed_data/{cytoplasmic,nuclear}_expression_matrix.npz are built by
code/1_preprocessing/1_clean.ipynb and code/2_outcome/1_split_expression.ipynb. This script

  A  rebuilds them independently from raw_data/transcripts.parquet (10x columns cell_id, overlaps_nucleus,
     feature_name, qv, is_gene; nuclear = overlaps_nucleus == 1, cytoplasmic = overlaps_nucleus == 0, tumor cells of
     cell_ids.npy, genes of gene_ids.npy) and compares them entry by entry with the stored matrices; with Q20 applied
     the sum of the two must equal the 10x count matrix (adata.X);
  B  measures the transcripts below Q20 (the stored matrices have no quality filter), by compartment and gene;
  C  measures the effect on the gene selection of two corrections, alone and together, against the reference "all
     tumor cells with all transcripts": removing the tumor cells with zero nuclear transcripts (the cells without a
     segmented nucleus plus a few with an empty small nucleus; all their transcripts are cytoplasmic), and keeping
     only transcripts with Q20 or better. All four variants use the protein-coding genes. Decided 2026-11-18: the
     cells are removed (ncd_ratio_utils.load_counts), Q20 is not applied; so the removal variant is the analysis,
     and the script checks that it reproduces steps 1 and 3. For each variant: the gene filter (step 1 rule) and
     the V-abs ranking (step 2 rule, variance of r on detected cells; no null needed);
  D  measures how far the score follows cell geometry, on the cells of the analysis: correlation of the cell's
     nuclear share of counts with the nucleus share of the cell area, and the share of the variance of r of the
     primary target genes explained by each.

Run from this folder (AIVC-env, about 4 min):  python check_upstream.py [--quick] [--keep-cache]
Needs: raw_data/transcripts.parquet, intermediate_data/adata.h5ad and the processed_data matrices of the five
samples; the outputs of steps 1 and 3 (filtered_genes.npy, target_genes.npy).
Output: output/1_exploration/ncd_ratio/upstream_check/ (emptied first; --quick: first sample, part of the raw file,
written to upstream_check_quick/)
  rebuild_check.csv      A, per sample
  qv_by_compartment.csv, qv_by_gene.csv      B
  variant_summary.csv, variant_by_sample.csv, variant_genes.csv      C
  geometry.csv, geometry_by_gene.csv      D
  summary.txt, upstream_overview.png
The Q20 matrices are cached in <output>/cache/ during the run and deleted at the end unless --keep-cache.
"""
import argparse
import os
import shutil
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy import sparse
from scipy.stats import spearmanr

import ncd_ratio_utils as nr

QV_MIN = 20
MIN_TOTAL = 50   # cells with at least this many transcripts enter the geometry correlations
ALL, REMOVED = "all tumor cells", "cells without nuclear transcripts removed"
VARIANTS = {ALL: (False, False), REMOVED: (True, False), "Q20": (False, True), "both": (True, True)}   # (remove cells, Q20)
QUICK_GROUPS = 12


def rebuild(sample, cells, genes, quick):
    """Nuclear and cytoplasmic count matrices of the tumor cells from the raw 10x transcripts, split by quality:
    {(compartment, quality): csr}, compartment in 'n', 'c', quality in 'hi' (qv >= QV_MIN), 'lo'."""
    pf = pq.ParquetFile(f"../../../data/{sample}/raw_data/transcripts.parquet")
    cell_set, gene_set = pa.array(list(cells), type=pa.large_string()), pa.array(list(genes), type=pa.large_string())
    G = len(genes)
    keys = {k: [] for k in [("n", "hi"), ("n", "lo"), ("c", "hi"), ("c", "lo")]}
    groups = range(min(QUICK_GROUPS, pf.metadata.num_row_groups) if quick else pf.metadata.num_row_groups)
    n_gene_tx = 0
    for g in groups:
        t = pf.read_row_group(g, columns=["cell_id", "overlaps_nucleus", "feature_name", "qv", "is_gene"])
        t = t.filter(t["is_gene"])
        n_gene_tx += t.num_rows
        row = pc.index_in(t["cell_id"], value_set=cell_set).to_numpy(zero_copy_only=False)
        col = pc.index_in(t["feature_name"], value_set=gene_set).to_numpy(zero_copy_only=False)
        ok = ~np.isnan(row) & ~np.isnan(col)
        key = row[ok].astype(np.int64) * G + col[ok].astype(np.int64)
        nuc = t["overlaps_nucleus"].to_numpy()[ok] == 1
        hi = t["qv"].to_numpy()[ok] >= QV_MIN
        for (comp, q) in keys:
            m = (nuc if comp == "n" else ~nuc) & (hi if q == "hi" else ~hi)
            u, cnt = np.unique(key[m], return_counts=True)
            keys[comp, q].append((u, cnt))
    out = {}
    for k, parts in keys.items():
        u = np.concatenate([p[0] for p in parts]); cnt = np.concatenate([p[1] for p in parts])
        out[k] = sparse.coo_matrix((cnt, (u // G, u % G)), shape=(len(cells), G)).tocsr()   # duplicates are summed
    return out, n_gene_tx


def moments(Xc, Xn, keep_cells, cols_keep):
    """Per gene of cols_keep, on the detected cells of keep_cells: cells detected, mean and variance of r."""
    Xc, Xn = Xc[keep_cells][:, cols_keep].tocsr(), Xn[keep_cells][:, cols_keep].tocsr()
    _, cols, c, n = nr.pairs(Xc, Xn)
    G = int(cols_keep.sum())
    r = nr.ratio_score(c, n)
    det = np.bincount(cols, minlength=G)
    mean = np.bincount(cols, weights=r, minlength=G) / np.maximum(det, 1)
    var = np.bincount(cols, weights=r * r, minlength=G) / np.maximum(det, 1) - mean ** 2
    return det, mean, var


def gene_corr(cols, r, x, G):
    """Per gene: squared Pearson correlation of r with the cell-level number x over the gene's pairs."""
    cnt = np.maximum(np.bincount(cols, minlength=G), 1).astype(float)
    s = lambda w: np.bincount(cols, weights=w, minlength=G) / cnt
    cov = s(r * x) - s(r) * s(x)
    vr, vx = s(r * r) - s(r) ** 2, s(x * x) - s(x) ** 2
    return np.where((vr > 1e-12) & (vx > 1e-12), cov ** 2 / np.maximum(vr * vx, 1e-24), np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="first sample, part of the raw file; no exact comparison")
    ap.add_argument("--keep-cache", action="store_true", help="keep the Q20 matrices in <output>/cache/")
    args = ap.parse_args()
    t0 = time.time()
    out = nr.step_dir("upstream_check", args.quick)
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out + "cache/")
    samples = nr.SAMPLES[:1] if args.quick else nr.SAMPLES
    shorts = [nr.short(s) for s in samples]
    filtered = nr.load_filtered_genes()
    target = nr.load_target_genes()                        # the primary list
    pool = nr.load_target_genes(nr.LIST_SIZE_POOL)
    coding_names = nr.coding_genes()

    # ---- pass 1: rebuild, checks A and B, expression of every gene under the four variants, geometry
    rows_a, rows_b, rows_d, qv_gene, geo_gene = [], [], [], [], []
    mean_count = {v: [] for v in VARIANTS}
    has_nuc, genes = {}, None
    hex_data = None
    for sample in samples:
        s = nr.short(sample)
        Xc, Xn, cells, genes = nr.load_counts(sample, all_cells=True, all_genes=True)
        ad = anndata.read_h5ad(f"../../../data/{sample}/intermediate_data/adata.h5ad", backed="r")
        tumor = np.where((ad.obs["cell_type_merged"] == "Malignant cell").values)[0]
        obs = ad.obs.iloc[tumor]
        assert (obs["cell_id"].astype(str).values == cells).all() and len(set(cells)) == len(cells), "cell ids do not match the AnnData"
        assert (ad.var_names.astype(str).values == genes).all(), "genes do not match the AnnData"
        nuc_ok = ~obs["nucleus_area"].isna().values
        has_nuc[s] = np.asarray(Xn.sum(axis=1)).ravel() > 0   # the cell rule of the analysis: at least one nuclear transcript
        assert not (has_nuc[s] & ~nuc_ok).any(), "a cell without a nucleus area has nuclear transcripts"
        R, n_gene_tx = rebuild(sample, cells, genes, args.quick)
        Rn, Rc = R["n", "hi"] + R["n", "lo"], R["c", "hi"] + R["c", "lo"]
        X10 = ad.X[tumor]
        X10 = sparse.csr_matrix(X10)
        hi_sum = (R["n", "hi"] + R["c", "hi"]).astype(np.float64)
        tn, tc = np.asarray(Xn.sum(axis=1)).ravel(), np.asarray(Xc.sum(axis=1)).ravel()
        rows_a.append({
            "sample": s, "tumor_cells": len(cells), "gene_transcripts_in_raw_file": n_gene_tx,
            "transcripts_of_tumor_cells": int(Rn.sum() + Rc.sum()),
            "nuclear_entries_differing": int((Rn != Xn).nnz), "cytoplasmic_entries_differing": int((Rc != Xc).nnz),
            "nuclear_sum_stored": int(Xn.sum()), "nuclear_sum_rebuilt": int(Rn.sum()),
            "cytoplasmic_sum_stored": int(Xc.sum()), "cytoplasmic_sum_rebuilt": int(Rc.sum()),
            "q20_entries_differing_from_10x_matrix": int((hi_sum != X10.astype(np.float64)).nnz),
            "q20_sum": int(hi_sum.sum()), "sum_10x_matrix": int(X10.sum()),
            "cells_without_nucleus": int((~nuc_ok).sum()), "share_without_nucleus": (~nuc_ok).mean(),
            "cells_with_zero_nuclear_count": int((tn == 0).sum()), "no_nucleus_cells_with_nuclear_transcripts": int(((tn > 0) & ~nuc_ok).sum()),
            "cells_with_zero_cytoplasmic_count": int((tc == 0).sum()),
            "transcripts_in_no_nucleus_cells": (tc + tn)[tn == 0].sum() / (tc + tn).sum()})
        lo_n, lo_c, all_n, all_c = R["n", "lo"].sum(), R["c", "lo"].sum(), Rn.sum(), Rc.sum()
        rows_b.append({"sample": s, "below_q20_nuclear": lo_n / all_n, "below_q20_cytoplasmic": lo_c / all_c,
                       "below_q20_all": (lo_n + lo_c) / (all_n + all_c),
                       "cytoplasmic_share_all_transcripts": all_c / (all_n + all_c),
                       "cytoplasmic_share_q20": R["c", "hi"].sum() / (R["c", "hi"].sum() + R["n", "hi"].sum())})
        gn = lambda M: np.asarray(M.sum(axis=0)).ravel()
        qv_gene.append(pd.DataFrame({"sample": s, "gene": genes, "transcripts": gn(Rn) + gn(Rc),
                                     "below_q20": (gn(R["n", "lo"]) + gn(R["c", "lo"])) / np.maximum(gn(Rn) + gn(Rc), 1),
                                     "below_q20_nuclear": gn(R["n", "lo"]) / np.maximum(gn(Rn), 1),
                                     "below_q20_cytoplasmic": gn(R["c", "lo"]) / np.maximum(gn(Rc), 1)}))
        sparse.save_npz(out + f"cache/{s}_c.npz", R["c", "hi"]); sparse.save_npz(out + f"cache/{s}_n.npz", R["n", "hi"])
        for v, (drop, q20) in VARIANTS.items():
            keep = has_nuc[s] if drop else np.ones(len(cells), bool)
            A = (R["c", "hi"] + R["n", "hi"]) if q20 else (Xc + Xn)
            mc = np.asarray(A[keep].sum(axis=0)).ravel() / keep.sum()
            mean_count[v].append(np.where(np.isin(genes, coding_names), mc, -1.0))   # non-coding genes rank last and are never kept

        # D: geometry, on the cells of the analysis (all of them have a nucleus area, asserted above)
        ana = has_nuc[s]
        share_n = tn / np.maximum(tc + tn, 1)
        area = (obs["nucleus_area"] / obs["cell_area"]).values
        use = ana & (tc + tn >= MIN_TOTAL)
        cols_t = np.isin(genes, target)
        rows_idx, cols, c, n = nr.pairs(Xc[ana][:, cols_t].tocsr(), Xn[ana][:, cols_t].tocsr())
        rr = nr.ratio_score(c, n)
        r2_area = gene_corr(cols, rr, area[ana][rows_idx], int(cols_t.sum()))
        r2_share = gene_corr(cols, rr, share_n[ana][rows_idx], int(cols_t.sum()))
        geo_gene.append(pd.DataFrame({"sample": s, "gene": genes[cols_t], "r2_area_ratio": r2_area, "r2_nuclear_share_of_counts": r2_share}))
        rows_d.append({"sample": s, "cells_used": int(use.sum()),
                       "nucleus_share_of_cell_area_median": np.nanmedian(area[use]),
                       "nuclear_share_of_counts_median": np.median(share_n[use]),
                       "pearson_nuclear_share_vs_area_ratio": np.corrcoef(share_n[use], area[use])[0, 1],
                       "spearman": spearmanr(share_n[use], area[use])[0],
                       "target_genes_median_r2_with_area_ratio": np.nanmedian(r2_area),
                       "target_genes_median_r2_with_nuclear_share_of_counts": np.nanmedian(r2_share)})
        if hex_data is None:
            hex_data = (s, area[use], share_n[use])
        print(f"[{s}] rebuilt: {rows_a[-1]['nuclear_entries_differing']} nuclear and {rows_a[-1]['cytoplasmic_entries_differing']} cytoplasmic "
              f"entries differ; Q20 vs 10x matrix: {rows_a[-1]['q20_entries_differing_from_10x_matrix']} ({time.time() - t0:.0f} s)", flush=True)
        del R, Rn, Rc, X10, hi_sum, Xc, Xn
    check_a, check_b, geometry = pd.DataFrame(rows_a), pd.DataFrame(rows_b), pd.DataFrame(rows_d)
    check_a.to_csv(out + "rebuild_check.csv", index=False)
    check_b.to_csv(out + "qv_by_compartment.csv", index=False)
    pd.concat(qv_gene).to_csv(out + "qv_by_gene.csv", index=False)
    geometry.to_csv(out + "geometry.csv", index=False)
    pd.concat(geo_gene).to_csv(out + "geometry_by_gene.csv", index=False)

    # ---- pass 2: gene filter and V-abs ranking per variant
    G = len(genes)
    need = min(nr.MIN_SAMPLES, len(samples))
    kept, stats = {}, {}
    for v in VARIANTS:
        ranks = np.column_stack([pd.Series(-m).rank(method="first").values for m in mean_count[v]])
        kept[v] = (np.sort(ranks, axis=1)[:, need - 1] <= nr.N_FILTER) & np.isin(genes, coding_names)
    union = np.any(list(kept.values()), axis=0)
    for v, (drop, q20) in VARIANTS.items():
        det, mean, var = [], [], []
        for sample in samples:
            s = nr.short(sample)
            if q20:
                Xc, Xn = sparse.load_npz(out + f"cache/{s}_c.npz"), sparse.load_npz(out + f"cache/{s}_n.npz")
            else:
                Xc, Xn, _, _ = nr.load_counts(sample, all_cells=True, all_genes=True)
            keep = has_nuc[s] if drop else np.ones(Xc.shape[0], bool)
            d, m, w = moments(Xc, Xn, keep, union)
            det.append(d / keep.sum()); mean.append(m); var.append(w)
        stats[v] = {"det": np.column_stack(det), "mean": np.column_stack(mean), "var": np.column_stack(var)}
        print(f"[{v}] statistics done ({time.time() - t0:.0f} s)", flush=True)
    ug = genes[union]
    lists, rank_v = {}, {}
    for v in VARIANTS:
        k = kept[v][union]
        var = np.where(k[:, None], stats[v]["var"], np.nan)
        if len(samples) < nr.MIN_SAMPLES:
            entry = np.sort(nr.ranks_desc(var, np.ones_like(var, bool)), axis=1)[:, need - 1]
        else:
            entry = nr.third_best(nr.ranks_desc(var, np.ones_like(var, bool)))
        tie = -np.nanmedian(var, axis=1)
        lists[v] = {n: nr.top_genes(entry, tie, n) for n in (nr.LIST_SIZE, nr.LIST_SIZE_CHECK)}
        order = np.lexsort((tie, entry))
        rk = np.full(len(ug), np.inf)
        fin = np.isfinite(entry[order])
        rk[order[fin]] = np.arange(1, fin.sum() + 1)
        rank_v[v] = rk
    if not args.quick:   # the variant with the cell rule of the analysis must reproduce steps 1 and 3
        assert (genes[kept[REMOVED]] == filtered).all(), "the removal variant does not reproduce the gene filter"
        assert (ug[np.argsort(rank_v[REMOVED])][:nr.LIST_SIZE_POOL] == pool).all(), "the removal variant does not reproduce the target list"

    ref = ALL
    rows, rows_s = [], []
    for v in VARIANTS:
        both_kept = kept[ref][union] & kept[v][union]
        row = {"variant": v, "genes_kept_by_filter": int(kept[v].sum()), "kept_also_in_reference": int((kept[v] & kept[ref]).sum()),
               "filter_genes_lost": int((kept[ref] & ~kept[v]).sum()), "filter_genes_gained": int((kept[v] & ~kept[ref]).sum())}
        for n in (nr.LIST_SIZE, nr.LIST_SIZE_CHECK):
            row[f"v_abs_{n}_shared_with_reference"] = int((lists[v][n] & lists[ref][n]).sum())
        row["spearman_rank_v_abs_with_reference"] = spearmanr(rank_v[ref][both_kept], rank_v[v][both_kept])[0]
        t = lists[v][nr.LIST_SIZE]
        row["list_median_detection_rate"] = float(np.median(np.median(stats[v]["det"], axis=1)[t]))
        row["list_cytoplasmic_leaning"] = int((np.median(stats[v]["mean"] - np.nanmedian(np.where(kept[v][union][:, None], stats[v]["mean"], np.nan), axis=0), axis=1)[t] > 0).sum())
        rows.append(row)
        for j, s in enumerate(shorts):
            rows_s.append({"variant": v, "sample": s,
                           "median_mean_r": np.median(stats[v]["mean"][kept[v][union], j]), "median_var_r": np.median(stats[v]["var"][kept[v][union], j]),
                           "spearman_mean_r_with_reference": spearmanr(stats[ref]["mean"][both_kept, j], stats[v]["mean"][both_kept, j])[0],
                           "spearman_var_r_with_reference": spearmanr(stats[ref]["var"][both_kept, j], stats[v]["var"][both_kept, j])[0],
                           "median_change_of_mean_r": np.median((stats[v]["mean"] - stats[ref]["mean"])[both_kept, j]),
                           "median_ratio_of_var_r": np.median((stats[v]["var"] / stats[ref]["var"])[both_kept, j])})
    variant_summary, variant_by_sample = pd.DataFrame(rows), pd.DataFrame(rows_s)
    variant_summary.to_csv(out + "variant_summary.csv", index=False)
    variant_by_sample.to_csv(out + "variant_by_sample.csv", index=False)
    vg = pd.DataFrame({"gene": ug})
    for v in VARIANTS:
        vg[f"kept [{v}]"] = kept[v][union]
        vg[f"rank_v_abs [{v}]"] = rank_v[v]
        vg[f"median_var_r [{v}]"] = np.median(stats[v]["var"], axis=1)
        vg[f"median_mean_r [{v}]"] = np.median(stats[v]["mean"], axis=1)
    vg.sort_values(f"rank_v_abs [{ref}]").to_csv(out + "variant_genes.csv", index=False)

    # ---- figure
    fig, ax = plt.subplots(2, 3, figsize=(18, 10.5))
    a = ax[0, 0]
    x = np.arange(len(shorts)); w = 0.38
    a.bar(x - w / 2, check_a["cells_with_zero_nuclear_count"] / check_a["tumor_cells"], w, color="#7570b3", label="tumor cells without nuclear transcripts")
    a.bar(x + w / 2, check_a["transcripts_in_no_nucleus_cells"], w, color="#bbbbbb", label="their share of the transcripts")
    a.set_xticks(x); a.set_xticklabels(shorts)
    a.set(ylabel="share", title="Tumor cells without nuclear transcripts (removed from the analysis)")
    a.legend(frameon=False)
    a = ax[0, 1]
    a.bar(x - w / 2, check_b["below_q20_nuclear"], w, color="#1b9e77", label="nuclear")
    a.bar(x + w / 2, check_b["below_q20_cytoplasmic"], w, color="#d95f02", label="cytoplasmic")
    a.set_xticks(x); a.set_xticklabels(shorts)
    a.set(ylabel="share of the transcripts below Q20", title="Transcripts below Q20 by compartment")
    a.legend(frameon=False)
    a = ax[0, 2]
    s0, ar, sh = hex_data
    hb = a.hexbin(ar, sh, gridsize=60, cmap="viridis", mincnt=1, bins="log")
    a.plot([0, 1], [0, 1], color="w", lw=0.8, ls="--")
    a.set(xlabel="nucleus share of the cell area", ylabel="nuclear share of the cell's transcripts",
          title=f"Cell geometry and the cell-wide lean ({s0}; r = {geometry['pearson_nuclear_share_vs_area_ratio'].iloc[0]:.2f})")
    fig.colorbar(hb, ax=a, label="cells")
    for a, v in zip(ax[1], [REMOVED, "Q20", "both"]):
        both_kept = kept[ref][union] & kept[v][union]
        a.scatter(rank_v[ref][both_kept], rank_v[v][both_kept], s=3, color="#999999")
        a.axvline(nr.LIST_SIZE, color="k", lw=0.6, ls=":"); a.axhline(nr.LIST_SIZE, color="k", lw=0.6, ls=":")
        r_ = variant_summary.set_index("variant").loc[v]
        a.set(xlabel="rank in V-abs, all tumor cells", ylabel=f"rank in V-abs, {v}",
              title=f"{'cells removed' if v == REMOVED else v}: {int(r_[f'v_abs_{nr.LIST_SIZE}_shared_with_reference'])} of {nr.LIST_SIZE} genes shared")
    fig.tight_layout()
    fig.savefig(out + "upstream_overview.png", dpi=150)
    plt.close(fig)

    # ---- summary
    pd.set_option("display.width", 320); pd.set_option("display.max_columns", 50)
    num = lambda v: f"{v:.0f}" if float(v).is_integer() else f"{v:.4f}"
    fmt = lambda d: d.set_index(d.columns[0]).T.to_string(float_format=num)
    size = nr.LIST_SIZE
    moved = lambda v, a_, b_: " ".join(ug[np.argsort(rank_v[a_])][np.isin(ug[np.argsort(rank_v[a_])], ug[lists[a_][size] & ~lists[b_][size]])][:60])
    lines = [f"NCD ratio analysis: audit of the compartment count matrices{' (QUICK RUN)' if args.quick else ''}", "",
             "A. Independent rebuild from raw_data/transcripts.parquet", fmt(check_a), "",
             f"B. Transcripts below Q{QV_MIN}", fmt(check_b), "",
             "C. Effect on the gene selection", fmt(variant_summary), "",
             "   per sample", variant_by_sample.to_string(index=False, float_format=lambda v: f"{v:.3f}"), ""]
    for v in list(VARIANTS)[1:]:
        lines += [f"   [{v}] genes leaving the V-abs list of {size} (first 60 by the rank with all tumor cells): " + moved(v, ref, v),
                  f"   [{v}] genes entering (first 60 by new rank): " + moved(v, v, ref), ""]
    lines += [f"D. Cell geometry (correlations: cells of the analysis with at least {MIN_TOTAL} transcripts; r2 per primary target gene over its detected cells of the analysis)", fmt(geometry), ""]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    if not args.keep_cache:
        shutil.rmtree(out + "cache/")
    print(f"done in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
