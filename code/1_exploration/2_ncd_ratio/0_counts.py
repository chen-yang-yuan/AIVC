"""NCD ratio analysis, step 0: counts check and choice of N (plans/ncd_ratio_exploration.md, sections 3 and 7).

How sparse is the ratio score r = log2((c + 1) / (n + 1)) on the tumor cells of the five samples, and how many of
the most expressed genes (N) should the gene filter keep?

Cells: the tumor cells with at least one nuclear transcript; genes: the protein-coding panel genes
(ncd_ratio_utils.load_counts).
Per sample and gene: expression (mean count per cell, cytoplasmic + nuclear), detection, cells with two or more
transcripts, pair types (cytoplasm-only, nucleus-only, two-sided), mean and variance of r over the detected cells.
Gene set of a given N (section 3): genes in the top N by mean count per tumor cell in >= 3 of 5 samples. Hence a
gene enters at entry_N = its third-best expression rank over the samples.
Measurable gene: >= MIN_CELLS_GE2 cells with two or more transcripts in >= 3 of 5 samples (what the variation
statistic of section 4.2 needs; single-transcript cells carry no information on variation).

Rule for N (fixed before the run):
  upper bound  largest N at which at least half of the last WINDOW genes that entered are measurable
  lower bound  smallest N that keeps KEEP_SHARE of the measurable genes of the known gene groups
  N by the rule = upper bound rounded down to ROUND if it is >= the lower bound; otherwise both are reported.
Decision after the review of the run (2026-11-18): 2000, 2500 and 3000 are all reasonable; the analysis proceeds
with N = 3000 (ncd_ratio_utils.N_FILTER), which the summary and the figure mark next to the value of the rule.

Run from this folder (AIVC-env, about 40 s):  python 0_counts.py [--quick]
Output: output/1_exploration/ncd_ratio/0_counts/ (emptied first; --quick writes to 0_counts_quick/).
"""
import argparse
import os
import shutil
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # shared subtype_utils.py
import ncd_ratio_utils as nr
import subtype_utils as su

N_GRID = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, None]   # None: all genes, set at run time
MIN_SAMPLES = nr.MIN_SAMPLES
MIN_CELLS_GE2 = nr.MIN_CELLS_GE2
WINDOW = 250
KEEP_SHARE = 0.90
ROUND = 500
CAP = 5   # counts are capped at CAP in the (c, n) table
TARGET_634 = "../../../output/1_exploration/enrichment/overlap_genes_3plus_samples.npy"
TOTAL_BINS = ["1", "2", "3-4", "5+"]


def sample_counts(sample, quick, rng):
    """Per-gene table, the (gene, c, n) count table with c, n capped at CAP, and the pairs for per-cell counts."""
    Xc, Xn, cells, genes = nr.load_counts(sample, nr.QUICK_CELLS if quick else None, rng)
    rows, cols, c, n = nr.pairs(Xc, Xn)
    assert c.sum() == Xc.sum() and n.sum() == Xn.sum()
    G, n_cells = len(genes), len(cells)
    key = np.minimum(c, CAP) * (CAP + 1) + np.minimum(n, CAP)
    K = np.bincount(cols * (CAP + 1) ** 2 + key, minlength=G * (CAP + 1) ** 2).reshape(G, CAP + 1, CAP + 1)
    r = nr.ratio_score(c, n)
    detected = np.bincount(cols, minlength=G)
    sum_r = np.bincount(cols, weights=r, minlength=G)
    sum_r2 = np.bincount(cols, weights=r * r, minlength=G)
    mean_r = sum_r / np.maximum(detected, 1)
    total = np.bincount(cols, weights=c + n, minlength=G)
    cyto = np.bincount(cols, weights=c, minlength=G)
    table = pd.DataFrame({
        "sample": nr.short(sample), "gene": genes, "n_cells": n_cells, "mean_count": total / n_cells,
        "detection_rate": detected / n_cells, "cells_detected": detected,
        "cells_ge2": detected - K[:, 1, 0] - K[:, 0, 1],
        "pairs_cyto_only": K[:, 1:, 0].sum(axis=1), "pairs_nuc_only": K[:, 0, 1:].sum(axis=1),
        "pairs_two_sided": K[:, 1:, 1:].sum(axis=(1, 2)), "pairs_nonzero_score": detected - np.bincount(cols[c == n], minlength=G),
        "cytoplasmic_fraction": cyto / np.maximum(total, 1), "mean_r_detected": mean_r,
        "var_r_detected": sum_r2 / np.maximum(detected, 1) - mean_r ** 2, "mean_r_all_cells": sum_r / n_cells})
    table["rank"] = table["mean_count"].rank(ascending=False, method="first").astype(int)
    return table, K, (rows, cols, c != n), n_cells


def total_bin_masks():
    """Boolean (CAP+1, CAP+1) masks of the capped (c, n) table for the total-count bins (exact: a capped value
    means a total of at least CAP)."""
    cc, nn = np.meshgrid(np.arange(CAP + 1), np.arange(CAP + 1), indexing="ij")
    t = cc + nn
    capped = (cc == CAP) | (nn == CAP)
    return {"1": t == 1, "2": t == 2, "3-4": ((t == 3) | (t == 4)) & ~capped, "5+": (t >= 5) | capped}, cc, nn


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true", help=f"{nr.QUICK_CELLS} cells per sample, output in 0_counts_quick/")
    args = p.parse_args()
    t0 = time.time()
    out = nr.step_dir("0_counts", args.quick)
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    rng = np.random.default_rng(nr.SEED)

    tables, Ks, cell_pairs, n_cells = [], {}, {}, {}
    for sample in nr.SAMPLES:
        table, K, cp, nc = sample_counts(sample, args.quick, rng)
        s = nr.short(sample)
        tables.append(table); Ks[s] = K; cell_pairs[s] = cp; n_cells[s] = nc
        print(f"[{s}] {nc} cells of the analysis, {table.cells_detected.sum() / nc / len(table):.3f} of pairs non-empty ({time.time() - t0:.0f} s)", flush=True)
    genes = tables[0]["gene"].values
    assert all((t["gene"].values == genes).all() for t in tables), "gene order differs between samples"
    gene_table = pd.concat(tables, ignore_index=True)
    gene_table.to_csv(out + "gene_table.csv", index=False)
    shorts = [nr.short(s) for s in nr.SAMPLES]
    G = len(genes)
    N_GRID[-1] = G

    # consensus order: a gene is in the set of N iff its third-best rank is <= N
    ranks = np.column_stack([t["rank"].values for t in tables])
    entry = np.sort(ranks, axis=1)[:, MIN_SAMPLES - 1]
    ge2 = np.column_stack([t["cells_ge2"].values for t in tables])
    n_meas = (ge2 >= MIN_CELLS_GE2).sum(axis=1)
    measurable = n_meas >= MIN_SAMPLES
    det = np.column_stack([t["detection_rate"].values for t in tables])

    # known gene groups
    groups = {f"curated: {k}": v for k, v in su.CURATED_PROGRAMS.items()}
    groups.update({f"hallmark: {k.replace('HALLMARK_', '').replace('_', ' ').capitalize()}": v
                   for k, v in su.read_gmt(nr.UTILS_PATH + "hallmark_pathways_filtered.gmt").items()})
    groups["634 target genes"] = list(np.load(TARGET_634, allow_pickle=True))
    member = {k: np.isin(genes, v) for k, v in groups.items()}
    in_any = np.any([m for k, m in member.items() if k.startswith(("curated", "hallmark"))], axis=0)

    consensus = pd.DataFrame({"gene": genes, "entry_N": entry, "n_samples_measurable": n_meas, "measurable": measurable,
                              "median_detection_rate": np.median(det, axis=1), "median_cells_ge2": np.median(ge2, axis=1),
                              "median_mean_r_detected": np.median(np.column_stack([t["mean_r_detected"].values for t in tables]), axis=1),
                              "in_program_or_hallmark": in_any, "target_634": member["634 target genes"]})
    consensus.sort_values("entry_N").to_csv(out + "gene_consensus.csv", index=False)

    # per sample and N
    masks, cc, nn = total_bin_masks()
    rows_n, rows_comp = [], []
    for N in N_GRID:
        in_set = entry <= N
        marginal = in_set & (entry > N - WINDOW)
        for j, s in enumerate(shorts):
            K = Ks[s][in_set].sum(axis=0)
            pairs_all = n_cells[s] * in_set.sum()
            nonempty = K.sum()
            r_idx, c_idx, nz = cell_pairs[s]
            use = nz & in_set[c_idx]
            per_cell = np.bincount(r_idx[use], minlength=n_cells[s])
            use_det = in_set[c_idx]
            per_cell_det = np.bincount(r_idx[use_det], minlength=n_cells[s])
            rows_n.append({
                "N": N, "sample": s, "genes_in_set": in_set.sum(), "share_nonempty": nonempty / pairs_all,
                "share_cyto_only": K[1:, 0].sum() / pairs_all, "share_nuc_only": K[0, 1:].sum() / pairs_all,
                "share_two_sided": K[1:, 1:].sum() / pairs_all, "share_ge2_transcripts": (nonempty - K[1, 0] - K[0, 1]) / pairs_all,
                "share_nonzero_score": (nonempty - np.trace(K)) / pairs_all,
                "nonempty_cyto_only": K[1:, 0].sum() / nonempty, "nonempty_nuc_only": K[0, 1:].sum() / nonempty,
                "nonempty_two_sided": K[1:, 1:].sum() / nonempty, "nonempty_single_transcript": (K[1, 0] + K[0, 1]) / nonempty,
                "median_detected_genes_per_cell": np.median(per_cell_det), "median_nonzero_scores_per_cell": np.median(per_cell),
                "marginal_median_detection_rate": np.median(det[marginal, j]), "marginal_median_cells_ge2": np.median(ge2[marginal, j]),
                "measurable_genes_in_set": (in_set & measurable).sum(), "measurable_share_of_set": (in_set & measurable).mean() / in_set.mean(),
                "measurable_share_of_marginal": measurable[marginal].mean()})
            for b in TOTAL_BINS:
                m = masks[b]
                tot = K[m].sum()
                rows_comp.append({"N": N, "sample": s, "total_count": b, "share_of_nonempty": tot / nonempty,
                                  "cyto_only": K[m & (nn == 0)].sum() / tot, "nuc_only": K[m & (cc == 0)].sum() / tot,
                                  "two_sided": K[m & (cc > 0) & (nn > 0)].sum() / tot})
    by_n = pd.DataFrame(rows_n)
    by_n.to_csv(out + "by_n.csv", index=False)
    comp = pd.DataFrame(rows_comp)
    comp.to_csv(out + "pair_composition.csv", index=False)

    # retention of the known groups
    rows_ret = []
    for name, m in member.items():
        row = {"group": name, "genes_in_panel": int(m.sum()), "measurable": int((m & measurable).sum())}
        for N in N_GRID:
            row[f"kept_{N}"] = (m & (entry <= N)).sum() / max(m.sum(), 1)
        for N in N_GRID:
            row[f"measurable_kept_{N}"] = (m & measurable & (entry <= N)).sum() / max((m & measurable).sum(), 1)
        rows_ret.append(row)
    retention = pd.DataFrame(rows_ret)
    retention.to_csv(out + "retention.csv", index=False)

    # rule for N
    order = np.argsort(entry, kind="stable")
    meas_sorted = measurable[order].astype(float)
    rolling = pd.Series(meas_sorted).rolling(WINDOW).mean().values   # share measurable among the last WINDOW genes entered
    ok = np.where(rolling >= 0.5)[0]
    upper = int(ok.max() + 1) if len(ok) else 0
    upper_N = int(entry[order][upper - 1]) if upper else 0
    lower_groups = int(np.ceil(np.quantile(entry[in_any & measurable], KEEP_SHARE))) if (in_any & measurable).any() else 0
    lower_all = int(np.ceil(np.quantile(entry[measurable], KEEP_SHARE)))
    chosen = (upper_N // ROUND) * ROUND
    curve = pd.DataFrame({"genes_entered": np.arange(1, G + 1), "entry_N": entry[order], "measurable": measurable[order],
                          "share_measurable_last_window": rolling, "cumulative_measurable": np.cumsum(meas_sorted)})
    curve.to_csv(out + "measurable_curve.csv", index=False)

    # figure
    colors = dict(zip(shorts, plt.cm.tab10.colors))
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    a = ax[0, 0]
    for s, t in zip(shorts, tables):
        d = t.sort_values("rank")
        a.plot(d["rank"], d["detection_rate"].rolling(50, center=True, min_periods=1).median(), color=colors[s], label=s)
    a.set(yscale="log", xlabel="expression rank of the gene within the sample", ylabel="detection rate (share of tumor cells)",
          title="Detection rate by expression rank")
    a.legend(frameon=False)
    a = ax[0, 1]
    for s, t in zip(shorts, tables):
        d = t.sort_values("rank")
        a.plot(d["rank"], d["cells_ge2"].rolling(50, center=True, min_periods=1).median().clip(lower=1), color=colors[s], label=s)
    a.axhline(MIN_CELLS_GE2, color="k", ls="--", lw=0.8)
    a.set(yscale="log", xlabel="expression rank of the gene within the sample", ylabel="cells with two or more transcripts",
          title=f"Cells that inform on variation (dashed: {MIN_CELLS_GE2})")
    a = ax[0, 2]
    a.plot(curve["genes_entered"], curve["share_measurable_last_window"], color="k")
    a.axhline(0.5, color="k", ls="--", lw=0.8)
    a.axvline(chosen, color="crimson", lw=1, ls="--")
    a.axvline((entry <= nr.N_FILTER).sum(), color="crimson", lw=1)
    a.set(xlabel="genes in the set (order of entry)", ylabel=f"share measurable among the last {WINDOW} genes",
          title=f"Measurable genes at the margin\n(rule {chosen} dashed, used {nr.N_FILTER})")
    a = ax[1, 0]
    for s in shorts:
        d = by_n[by_n["sample"] == s]
        a.plot(d["N"], d["share_nonempty"], color=colors[s], label=s)
        a.plot(d["N"], d["share_ge2_transcripts"], color=colors[s], ls=":")
    a.set(yscale="log", xlabel="N", ylabel="share of (gene, cell) pairs", title="Non-empty pairs (solid) and pairs with >= 2 transcripts (dotted)")
    a = ax[1, 1]
    ref = by_n[by_n["N"] == nr.N_FILTER].set_index("sample").loc[shorts]
    bottom = np.zeros(len(shorts))
    for col, lab, colr in [("nonempty_cyto_only", "cytoplasm only", "#d95f02"), ("nonempty_nuc_only", "nucleus only", "#1b9e77"),
                           ("nonempty_two_sided", "both compartments", "#7570b3")]:
        a.bar(shorts, ref[col].values, bottom=bottom, label=lab, color=colr)
        bottom += ref[col].values
    a.set(ylabel="share of the non-empty pairs", title=f"Types of non-empty pairs (N = {int(ref['N'].iloc[0])})")
    a.set_ylim(0, 1.15)
    a.legend(frameon=False, loc="upper center", ncol=3, fontsize=8)
    a = ax[1, 2]
    show = [k for k in member if k.startswith("curated")] + ["634 target genes"]
    for k in show:
        row = retention[retention["group"] == k].iloc[0]
        a.plot(N_GRID, [row[f"kept_{N}"] for N in N_GRID], marker="o", ms=3, lw=1, label=k.replace("curated: ", ""))
    a.axvline(chosen, color="crimson", lw=1, ls="--")
    a.axvline(nr.N_FILTER, color="crimson", lw=1)
    a.set(xlabel="N", ylabel="share of the group's coding panel genes kept", title="Known gene groups kept by the filter")
    a.legend(frameon=False, fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(out + "counts_overview.png", dpi=150)
    plt.close(fig)

    # summary
    pd.set_option("display.width", 250)
    fmt = lambda d: d.to_string(index=False, float_format=lambda x: f"{x:.3f}")
    cols_n = ["N", "sample", "genes_in_set", "share_nonempty", "share_ge2_transcripts", "share_nonzero_score", "nonempty_cyto_only",
              "nonempty_nuc_only", "nonempty_two_sided", "nonempty_single_transcript", "median_nonzero_scores_per_cell",
              "marginal_median_detection_rate", "marginal_median_cells_ge2", "measurable_genes_in_set", "measurable_share_of_marginal"]
    lines = [f"NCD ratio analysis, step 0: counts check{' (QUICK RUN)' if args.quick else ''}", "",
             "Cells of the analysis (tumor cells with at least one nuclear transcript): " + ", ".join(f"{s} {n_cells[s]}" for s in shorts),
             f"Measurable gene: >= {MIN_CELLS_GE2} cells with two or more transcripts in >= {MIN_SAMPLES} of 5 samples: "
             f"{measurable.sum()} of the {G} protein-coding panel genes", "",
             "Rule for N",
             f"  upper bound (at least half of the last {WINDOW} genes entered are measurable): {upper_N}",
             f"  lower bound (keeps {KEEP_SHARE:.0%} of the measurable genes of the curated and hallmark groups): {lower_groups}",
             f"  for reference, N that keeps {KEEP_SHARE:.0%} of all measurable genes: {lower_all}",
             f"  N by the rule (upper bound rounded down to {ROUND}): {chosen}"
             + ("" if chosen >= lower_groups else "  [CONFLICT with the lower bound]"),
             f"  N used by the analysis (decision after review; 2000, 2500 and 3000 are all reasonable): {nr.N_FILTER} "
             f"-> {(entry <= nr.N_FILTER).sum()} genes, {((entry <= nr.N_FILTER) & measurable).sum()} measurable", "",
             "Per sample and N", fmt(by_n[cols_n]), "",
             "Non-empty pairs by total count (shares within the bin)", fmt(comp[comp["N"].isin([2000, 2500, 3000])]), "",
             "Known gene groups: share of the group's protein-coding panel genes kept",
             fmt(retention[["group", "genes_in_panel", "measurable"] + [f"kept_{N}" for N in N_GRID]]), "",
             "Known gene groups: share of the group's measurable genes kept",
             fmt(retention[["group", "measurable"] + [f"measurable_kept_{N}" for N in N_GRID]]), ""]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    print(f"done in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
