"""NCD ratio analysis, step 2: selection statistics and the candidate target gene lists
(plans/ncd_ratio_exploration.md, section 4).

Statistics per sample and gene, on the tumor cells in which the gene is detected (c + n >= 1), for the genes kept
by the filter of step 1:
  mean_r, mean_r_centred   mean of r = log2((c + 1) / (n + 1)); centred = minus the sample's median over the genes
  var_r                    variance of r (absolute variation)
  null_var_mean, _sd       the same variance under the count-preserving null (ncd_ratio_utils.null_counts: every
                           pair keeps its total, transcripts reassigned with the gene's cytoplasmic fraction in the
                           sample), mean and SD over N_NULL draws
  excess                   var_r / null_var_mean (size)
  z                        (var_r - null_var_mean) / null_var_sd (evidence)
  excess_t2                the same ratio on the pairs with exactly two transcripts (exact binomial null with the
                           gene's fraction on these pairs): compares all genes at the same count. A check, not used
                           for the ranking.
Single-transcript cells carry no information on variation (their variance is fixed by the mean), so the excess
rests on the cells with two or more transcripts (cells_ge2).

The statistics are computed by ncd_ratio_utils.pair_stats, the rankings by ncd_ratio_utils.entry_ranks.

Three ranked lists. A gene qualifies in a sample if it is in the sample's top K by the statistic; it is selected if
it qualifies in >= MIN_SAMPLES samples. Hence the rank of a gene is its third-best within-sample rank (entry_K), and
the top n genes by entry_K are the list of size n. Genes with the same entry_K are ordered by the median of the
statistic over the samples (larger first).
  V-null  statistic excess; a gene is ranked in a sample only if cells_ge2 >= MIN_CELLS_GE2 and z >= Z_MIN there
  V-abs   statistic var_r; all genes
  M       statistic |mean_r_centred|; all genes; the qualifying samples must agree on the direction
The choice between the lists is step 3 (3_audit.py): V-abs is the target list (ncd_ratio_utils.TARGET_LIST), M is a
reference only. List size: ncd_ratio_utils.LIST_SIZE (500; larger nested check LIST_SIZE_CHECK = 1000).
direction: sign of the median over samples of mean_r_centred (cytoplasmic if > 0); direction_m: the direction
under which the gene enters list M.

Run from this folder (AIVC-env, a few minutes):  python 2_selection.py [--quick]
Output: output/1_exploration/ncd_ratio/2_selection/ (emptied first; --quick writes to 2_selection_quick/)
  gene_stats.csv      per sample and gene
  gene_summary.csv    per gene: medians over samples, entry_K and rank under the three rules, direction
  list_v_null.csv, list_v_abs.csv, list_m.csv     the first TOP genes of each list, ranked
  list_overlap.csv, groups_in_lists.csv, summary.txt, selection_overview.png
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
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # shared subtype_utils.py
import ncd_ratio_utils as nr
import subtype_utils as su

TOP = nr.LIST_SIZE_POOL   # genes written to the list files
SHOW = nr.LIST_SIZE      # list size shown in the figure: the primary list
GROUP_SIZES = (nr.LIST_SIZE, nr.LIST_SIZE_CHECK)
SIZES = [250, 500, 750, 1000]
TARGET_634 = "../../../output/1_exploration/enrichment/overlap_genes_3plus_samples.npy"
LISTS = {"v_null": "V-null", "v_abs": "V-abs", "m": "M"}


def sample_stats(sample, kept, quick, rng):
    Xc, Xn, cells, genes = nr.load_counts(sample, nr.QUICK_CELLS if quick else None, rng)
    cols_keep = np.isin(genes, kept)
    assert (genes[cols_keep] == kept).all(), "filtered genes are not in panel order"
    _, cols, c, n = nr.pairs(Xc[:, cols_keep].tocsr(), Xn[:, cols_keep].tocsr())
    del Xc, Xn
    return pd.DataFrame({"sample": nr.short(sample), "gene": kept, "n_cells": len(cells),
                         **nr.pair_stats(cols, c, n, len(kept), len(cells), rng)})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true", help=f"{nr.QUICK_CELLS} cells per sample, output in 2_selection_quick/")
    args = p.parse_args()
    t0 = time.time()
    out = nr.step_dir("2_selection", args.quick)
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    rng = np.random.default_rng(nr.SEED)
    kept = nr.load_filtered_genes()
    shorts = [nr.short(s) for s in nr.SAMPLES]
    tables = []
    for sample in nr.SAMPLES:
        tables.append(sample_stats(sample, kept, args.quick, rng))
        t = tables[-1]
        print(f"[{nr.short(sample)}] {len(kept)} genes; z >= {nr.Z_MIN}: {(t.z >= nr.Z_MIN).sum()}; "
              f"max |null mean - observed mean| of r: {(t.null_mean_r - t.mean_r).abs().max():.4f} ({time.time() - t0:.0f} s)", flush=True)
    stats = pd.concat(tables, ignore_index=True)
    stats.to_csv(out + "gene_stats.csv", index=False)
    mat = lambda col: np.column_stack([t[col].values for t in tables])
    G = len(kept)

    # check against step 0 (full run only)
    if not args.quick:
        ref = pd.read_csv(nr.step_dir("0_counts") + "gene_table.csv").set_index(["sample", "gene"])
        ref = ref.loc[list(zip(stats["sample"], stats["gene"]))]
        assert np.allclose(ref["mean_r_detected"].values, stats["mean_r"].values) and np.allclose(ref["var_r_detected"].values, stats["var_r"].values)
        print("mean and variance of r reproduce step 0")

    excess, z, ge2, var, cen = mat("excess"), mat("z"), mat("cells_ge2"), mat("var_r"), mat("mean_r_centred")
    enough = ge2 >= nr.MIN_CELLS_GE2
    passing = enough & (z >= nr.Z_MIN)
    entry, direction_m = nr.entry_ranks(excess, z, ge2, var, cen)
    med = lambda a: np.nanmedian(a, axis=1)
    summary = pd.DataFrame({
        "gene": kept, "median_detection_rate": med(mat("detection_rate")), "median_cells_ge2": med(ge2),
        "median_cytoplasmic_fraction": med(mat("cytoplasmic_fraction")), "median_mean_r": med(mat("mean_r")),
        "median_mean_r_centred": med(cen), "n_samples_positive": (cen > 0).sum(axis=1),
        "median_var_r": med(var), "median_excess": med(excess), "median_z": med(z), "median_excess_t2": med(mat("excess_t2")),
        "n_samples_enough_cells": enough.sum(axis=1), "n_samples_z_pass": passing.sum(axis=1),
        "target_634": np.isin(kept, np.load(TARGET_634, allow_pickle=True))})
    tie = {"v_null": -summary["median_excess"].values, "v_abs": -summary["median_var_r"].values, "m": -np.abs(summary["median_mean_r_centred"].values)}
    for k in LISTS:
        summary[f"entry_K_{k}"] = entry[k]
        order = np.lexsort((tie[k], entry[k]))
        rk = np.full(G, np.inf)
        finite = np.isfinite(entry[k][order])
        rk[order[finite]] = np.arange(1, finite.sum() + 1)
        summary[f"rank_{k}"] = rk
    summary["direction_m"] = np.where(direction_m > 0, "cytoplasmic", "nuclear")
    summary["direction"] = np.where(summary["median_mean_r_centred"] > 0, "cytoplasmic", "nuclear")
    summary.to_csv(out + "gene_summary.csv", index=False)
    per_sample_cols = ["detection_rate", "cells_ge2", "mean_r_centred", "var_r", "excess", "z"]
    wide = pd.concat([pd.DataFrame({f"{c}_{s}": t[c].values for s, t in zip(shorts, tables)}) for c in per_sample_cols], axis=1)
    for k in LISTS:
        top = summary[f"rank_{k}"] <= TOP
        pd.concat([summary[top], wide[top]], axis=1).sort_values(f"rank_{k}").to_csv(out + f"list_{k}.csv", index=False)

    # overlaps and groups
    sets = lambda n: {**{LISTS[k]: set(kept[summary[f"rank_{k}"] <= n]) for k in LISTS}, "634 target genes": set(kept[summary["target_634"]])}
    rows = []
    for n in SIZES:
        s = sets(n)
        for a in s:
            rows.append({"size": n, "list": a, **{b: len(s[a] & s[b]) for b in s}})
    overlap = pd.DataFrame(rows)
    overlap.to_csv(out + "list_overlap.csv", index=False)
    groups = {f"curated: {k}": v for k, v in su.CURATED_PROGRAMS.items()}
    groups.update({f"hallmark: {k.replace('HALLMARK_', '').replace('_', ' ').capitalize()}": v
                   for k, v in su.read_gmt(nr.UTILS_PATH + "hallmark_pathways_filtered.gmt").items()})
    rows = []
    for name, g in groups.items():
        m = np.isin(kept, g)
        row = {"group": name, "filtered_genes": int(m.sum())}
        for n in GROUP_SIZES:
            for k in LISTS:
                row[f"{LISTS[k]}_{n}"] = int((m & (summary[f"rank_{k}"].values <= n)).sum())
        rows.append(row)
    group_table = pd.DataFrame(rows)
    group_table.to_csv(out + "groups_in_lists.csv", index=False)

    # figure
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    col = {"v_null": "#d95f02", "v_abs": "#1b9e77", "m": "#7570b3"}
    a = ax[0, 0]
    counts = [(enough[:, j] & (z[:, j] >= nr.Z_MIN)).sum() for j in range(5)]
    a.bar(shorts + ["3 of 5"], counts + [(passing.sum(axis=1) >= nr.MIN_SAMPLES).sum()], color=["#999999"] * 5 + [col["v_null"]])
    a.bar(shorts, [enough[:, j].sum() for j in range(5)], color="none", edgecolor="k", lw=0.8)
    a.set(ylabel="genes", title=f"Genes above the null (z >= {nr.Z_MIN:g}); outline: genes with enough cells")
    a = ax[0, 1]
    x, y = summary["median_detection_rate"], summary["median_excess"]
    a.scatter(x, y, s=4, color="#bbbbbb")
    top = summary["rank_v_null"] <= SHOW
    a.scatter(x[top], y[top], s=5, color=col["v_null"], label=f"V-null first {SHOW}")
    a.axhline(1, color="k", lw=0.6, ls="--")
    a.set(xscale="log", xlabel="detection rate (median over samples)", ylabel="excess: observed / null variance of r (median)",
          title="Excess over the null against expression")
    a.legend(frameon=False)
    a = ax[0, 2]
    for k in LISTS:
        e = np.sort(entry[k][np.isfinite(entry[k])])
        a.plot(np.arange(1, len(e) + 1), e, color=col[k], label=LISTS[k])
    for n in GROUP_SIZES:
        a.axvline(n, color="k", lw=0.6, ls=":")
    a.plot([0, G], [0, G], color="k", lw=0.5, ls="--")
    a.set(xlabel="list size n", ylabel="K: within-sample top K needed (in 3 of 5 samples)", title="Agreement of the samples (dashed: perfect agreement)")
    a.legend(frameon=False)
    a = ax[1, 0]
    x, y = summary["median_cytoplasmic_fraction"], summary["median_var_r"]
    a.scatter(x, y, s=4, color="#bbbbbb")
    top = summary["rank_v_abs"] <= SHOW
    a.scatter(x[top], y[top], s=5, color=col["v_abs"], label=f"V-abs first {SHOW}")
    a.set(xlabel="cytoplasmic fraction of the gene (median over samples)", ylabel="variance of r on detected cells (median)",
          title="Absolute variance against the balance of the gene")
    a.legend(frameon=False)
    a = ax[1, 1]
    bins = np.linspace(-1.2, 1.2, 61)
    a.hist(summary["median_mean_r_centred"], bins=bins, color="#bbbbbb", label="all filtered genes")
    for k in LISTS:
        a.hist(summary["median_mean_r_centred"][summary[f"rank_{k}"] <= SHOW], bins=bins, histtype="step", lw=1.5, color=col[k], label=f"{LISTS[k]} first {SHOW}")
    a.set(xlabel="centred mean score of the gene (median over samples); > 0 more cytoplasmic", ylabel="genes", title="Direction of the listed genes")
    a.legend(frameon=False)
    a = ax[1, 2]
    names = list(sets(SHOW).keys())
    o = overlap[overlap["size"] == SHOW].set_index("list").loc[names, names].values
    im = a.imshow(o, cmap="Blues")
    a.set_xticks(range(len(names))); a.set_xticklabels(names, rotation=20); a.set_yticks(range(len(names))); a.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            a.text(j, i, str(o[i, j]), ha="center", va="center", color="w" if o[i, j] > o.max() / 2 else "k")
    a.set(title=f"Genes shared by the lists (size {SHOW})")
    fig.tight_layout()
    fig.savefig(out + "selection_overview.png", dpi=150)
    plt.close(fig)

    # summary
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 200)
    fmt = lambda d, **k: d.to_string(float_format=lambda x: f"{x:.3f}", **k)
    per_sample = pd.DataFrame({
        "sample": shorts, "genes_enough_cells": enough.sum(axis=0), "genes_z_pass": passing.sum(axis=0),
        "median_excess_enough": [np.median(excess[enough[:, j], j]) for j in range(5)],
        "p90_excess_enough": [np.quantile(excess[enough[:, j], j], 0.9) for j in range(5)],
        "median_excess_t2": [np.nanmedian(mat("excess_t2")[:, j]) for j in range(5)],
        "rho_excess_detection": [spearmanr(excess[enough[:, j], j], mat("detection_rate")[enough[:, j], j])[0] for j in range(5)],
        "rho_excess_t2_detection": [spearmanr(mat("excess_t2")[:, j], mat("detection_rate")[:, j], nan_policy="omit")[0] for j in range(5)],
        "rho_var_detection": [spearmanr(var[:, j], mat("detection_rate")[:, j])[0] for j in range(5)],
        "median_mean_r": [np.median(mat("mean_r")[:, j]) for j in range(5)]})
    size_rows = []
    for n in SIZES:
        for k in LISTS:
            sel = summary[f"rank_{k}"] <= n
            size_rows.append({"size": n, "list": LISTS[k], "genes": int(sel.sum()), "K_needed": float(np.sort(entry[k])[min(n, G) - 1]),
                              "cytoplasmic": int((summary["direction"][sel] == "cytoplasmic").sum()), "nuclear": int((summary["direction"][sel] == "nuclear").sum()),
                              "median_detection_rate": summary["median_detection_rate"][sel].median(),
                              "median_excess": summary["median_excess"][sel].median(), "median_excess_t2": summary["median_excess_t2"][sel].median(),
                              "measurable_share": (summary["n_samples_enough_cells"][sel] >= nr.MIN_SAMPLES).mean()})
    size_table = pd.DataFrame(size_rows)
    rho_rank = {LISTS[k]: spearmanr(np.where(np.isfinite(entry[k]), entry[k], G + 1), summary["median_detection_rate"])[0] for k in LISTS}
    v = summary.sort_values("rank_v_null")
    curve = pd.DataFrame({"rank": [50, 100, 250, 500, 750, 1000, 1500, 2000]})
    curve["genes_available"] = int(np.isfinite(entry["v_null"]).sum())
    curve["median_excess_at_rank"] = [v["median_excess"].iloc[r - 1] if r <= np.isfinite(entry["v_null"]).sum() else np.nan for r in curve["rank"]]
    curve["median_z_at_rank"] = [v["median_z"].iloc[r - 1] if r <= np.isfinite(entry["v_null"]).sum() else np.nan for r in curve["rank"]]
    lines = [f"NCD ratio analysis, step 2: selection statistics and candidate lists{' (QUICK RUN)' if args.quick else ''}", "",
             f"Genes: {G} (filter of step 1). Null: {nr.N_NULL} draws per sample. Evidence threshold z >= {nr.Z_MIN:g}; "
             f"cells with two or more transcripts required: {nr.MIN_CELLS_GE2}.",
             f"Fewest detected cells of a gene in its third-best sample: {int(np.sort(mat('cells_detected'), axis=1)[:, -nr.MIN_SAMPLES].min())}", "",
             "Per sample (rho = Spearman over genes with the detection rate)", fmt(per_sample, index=False), "",
             f"V-null: genes above the null in >= {nr.MIN_SAMPLES} of 5 samples: {(passing.sum(axis=1) >= nr.MIN_SAMPLES).sum()} "
             f"(in 5: {(passing.sum(axis=1) == 5).sum()}, in 4: {(passing.sum(axis=1) == 4).sum()}, in 3: {(passing.sum(axis=1) == 3).sum()})",
             "V-null along its ranking", fmt(curve, index=False), "",
             "Lists by size (K_needed = within-sample top K that gives the size; direction by the centred mean score)", fmt(size_table, index=False), "",
             "Spearman of the gene's entry K with its detection rate (negative = the list favours highly expressed genes): "
             + ", ".join(f"{k} {r:.2f}" for k, r in rho_rank.items()), "",
             "Genes shared by the lists", fmt(overlap, index=False), ""]
    show = ["gene", "direction", "median_detection_rate", "median_cytoplasmic_fraction", "median_mean_r_centred", "median_var_r", "median_excess", "median_z", "median_excess_t2"]
    for k in LISTS:
        lines += [f"{LISTS[k]}: first 30 genes", fmt(summary.sort_values(f"rank_{k}")[show].head(30), index=False), ""]
    lines += [f"Known gene groups: genes in the lists (sizes {GROUP_SIZES[0]} and {GROUP_SIZES[1]})", fmt(group_table, index=False), ""]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    print(f"done in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
