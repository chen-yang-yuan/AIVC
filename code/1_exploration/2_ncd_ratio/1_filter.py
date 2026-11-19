"""NCD ratio analysis, step 1: gene filter (plans/ncd_ratio_exploration.md, sections 3 and 7).

Genes: the protein-coding panel genes; cells: tumor cells with at least one nuclear transcript (ncd_ratio_utils).
Rule: a gene is kept if it is in the top N_FILTER by mean count per cell (cytoplasmic + nuclear) in at least
MIN_SAMPLES of the five samples (ncd_ratio_utils.N_FILTER = 3000, MIN_SAMPLES = 3). Leakage is not considered.
The kept genes are the candidates of the target gene lists of step 2.

Input: the per-gene table of step 0 (0_counts/gene_table.csv, full run), so that the two steps agree by
construction. Run 0_counts.py first.

Run from this folder (AIVC-env, a few seconds):  python 1_filter.py
Output: output/1_exploration/ncd_ratio/1_filter/ (emptied first)
  filtered_genes.npy             the kept genes in panel order (read by ncd_ratio_utils.load_filtered_genes)
  gene_filter.csv                all protein-coding panel genes: kept or not, ranks and counts per sample
  dropped_sample_specific.csv    genes not kept although in the top N_FILTER of one or two samples
  kept_per_sample.csv            the kept genes per sample (detection, cells with two or more transcripts, direction)
  top_n_overlap.csv              genes shared by the per-sample top-N_FILTER lists
  groups_kept.csv                known gene groups (curated programs, hallmark sets, the 634 target genes) kept
  summary.txt, filter_overview.png
"""
import os
import shutil
import sys
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

TARGET_634 = "../../../output/1_exploration/enrichment/overlap_genes_3plus_samples.npy"
PER_SAMPLE = ["rank", "mean_count", "detection_rate", "cells_ge2", "cytoplasmic_fraction", "mean_r_detected"]
STRONG_RANKS = [500, 1000]   # a dropped gene is reported as strongly sample-specific below these ranks


def main():
    out = nr.step_dir("1_filter")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    N = nr.N_FILTER
    shorts = [nr.short(s) for s in nr.SAMPLES]
    table = pd.read_csv(nr.step_dir("0_counts") + "gene_table.csv")
    assert (table.groupby("sample")["n_cells"].first() > nr.QUICK_CELLS).all(), "0_counts/gene_table.csv is not from a full run"
    genes = table[table["sample"] == shorts[0]]["gene"].values
    wide = {}
    for s in shorts:
        t = table[table["sample"] == s].set_index("gene").loc[genes]
        for c in PER_SAMPLE:
            wide[f"{c}_{s}"] = t[c].values
    wide = pd.DataFrame(wide, index=genes)
    ranks = wide[[f"rank_{s}" for s in shorts]].values
    in_top = ranks <= N
    n_top = in_top.sum(axis=1)
    kept = n_top >= nr.MIN_SAMPLES
    ge2 = wide[[f"cells_ge2_{s}" for s in shorts]].values
    measurable = (ge2 >= nr.MIN_CELLS_GE2).sum(axis=1) >= nr.MIN_SAMPLES

    info = pd.DataFrame({"gene": genes, "kept": kept, "n_samples_top": n_top,
                         "samples_top": [",".join(np.array(shorts)[r]) for r in in_top],
                         "entry_N": np.sort(ranks, axis=1)[:, nr.MIN_SAMPLES - 1], "best_rank": ranks.min(axis=1),
                         "best_sample": np.array(shorts)[ranks.argmin(axis=1)], "measurable": measurable})
    full = pd.concat([info.reset_index(drop=True), wide.reset_index(drop=True)], axis=1)
    full.to_csv(out + "gene_filter.csv", index=False)
    np.save(out + "filtered_genes.npy", genes[kept])
    dropped = full[~kept & (n_top > 0)].sort_values("best_rank")
    dropped.to_csv(out + "dropped_sample_specific.csv", index=False)

    # known groups
    groups = {f"curated: {k}": v for k, v in su.CURATED_PROGRAMS.items()}
    groups.update({f"hallmark: {k.replace('HALLMARK_', '').replace('_', ' ').capitalize()}": v
                   for k, v in su.read_gmt(nr.UTILS_PATH + "hallmark_pathways_filtered.gmt").items()})
    groups["634 target genes"] = list(np.load(TARGET_634, allow_pickle=True))
    rows = []
    for name, g in groups.items():
        m = np.isin(genes, g)
        rows.append({"group": name, "panel_genes": m.sum(), "kept": (m & kept).sum(), "share_kept": (m & kept).sum() / max(m.sum(), 1),
                     "dropped_genes": ", ".join(genes[m & ~kept][np.argsort(ranks.min(axis=1)[m & ~kept])][:12])})
    group_table = pd.DataFrame(rows)
    group_table.to_csv(out + "groups_kept.csv", index=False)

    # per-sample description of the kept set
    by_n = pd.read_csv(nr.step_dir("0_counts") + "by_n.csv")
    by_n = by_n[by_n["N"] == N].set_index("sample")
    rows = []
    for j, s in enumerate(shorts):
        d, g2 = wide[f"detection_rate_{s}"].values[kept], ge2[kept, j]
        rows.append({"sample": s, "own_top_N_kept": (in_top[:, j] & kept).sum(), "kept_not_in_own_top_N": (kept & ~in_top[:, j]).sum(),
                     "detection_p10": np.quantile(d, 0.1), "detection_median": np.median(d), "detection_p90": np.quantile(d, 0.9),
                     "cells_ge2_p10": np.quantile(g2, 0.1), "cells_ge2_median": np.median(g2),
                     "genes_ge2_below_min": (g2 < nr.MIN_CELLS_GE2).sum(), "share_nonempty_pairs": by_n.loc[s, "share_nonempty"],
                     "cytoplasmic_fraction_median": np.median(wide[f"cytoplasmic_fraction_{s}"].values[kept]),
                     "mean_r_detected_median": np.median(wide[f"mean_r_detected_{s}"].values[kept]),
                     "genes_mean_r_positive": (wide[f"mean_r_detected_{s}"].values[kept] > 0).mean()})
    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(out + "kept_per_sample.csv", index=False)
    overlap = pd.DataFrame([[(in_top[:, a] & in_top[:, b]).sum() for b in range(5)] for a in range(5)], index=shorts, columns=shorts)
    overlap.to_csv(out + "top_n_overlap.csv")
    med_r = np.median(wide[[f"mean_r_detected_{s}" for s in shorts]].values, axis=1)
    sign_agree = ((wide[[f"mean_r_detected_{s}" for s in shorts]].values > 0).sum(axis=1))

    # figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    a = ax[0, 0]
    counts = np.bincount(n_top, minlength=6)
    a.bar(range(6), counts, color=["#bbbbbb"] * nr.MIN_SAMPLES + ["#1b9e77"] * (6 - nr.MIN_SAMPLES))
    for k, v in enumerate(counts):
        a.text(k, v + 20, str(v), ha="center", fontsize=9)
    a.set(xlabel=f"number of samples in which the gene is in the top {N}", ylabel="genes",
          title=f"Kept (green): {kept.sum()} of {len(genes)} genes")
    a = ax[0, 1]
    pos = np.arange(5)
    for off, mask, col, lab in [(-0.2, kept, "#1b9e77", "kept"), (0.2, ~kept, "#bbbbbb", "dropped")]:
        data = [np.log10(np.maximum(wide[f"detection_rate_{s}"].values[mask], 1e-5)) for s in shorts]
        bp = a.boxplot(data, positions=pos + off, widths=0.35, patch_artist=True, showfliers=False)
        for b in bp["boxes"]:
            b.set(facecolor=col)
        a.plot([], [], color=col, lw=6, label=lab)
    a.set_xticks(pos); a.set_xticklabels(shorts)
    a.set(ylabel="log10 detection rate (share of tumor cells)", title="Detection rate of kept and dropped genes")
    a.legend(frameon=False)
    a = ax[1, 0]
    best = ranks.min(axis=1)[~kept]
    a.hist(best, bins=np.arange(0, 5001, 100), color="#bbbbbb")
    a.axvline(N, color="crimson", lw=1)
    a.set(xlabel="best expression rank of the gene in any sample", ylabel="dropped genes",
          title=f"Dropped genes: {(best <= N).sum()} are in the top {N} of one or two samples")
    a = ax[1, 1]
    for j, s in enumerate(shorts):
        a.hist(wide[f"cytoplasmic_fraction_{s}"].values[kept], bins=np.linspace(0, 1, 41), histtype="step", lw=1.3, label=s)
    a.axvline(0.5, color="k", lw=0.6, ls="--")
    a.set(xlabel="cytoplasmic fraction of the gene's transcripts in the sample", ylabel="kept genes",
          title="Compartment balance of the kept genes")
    a.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out + "filter_overview.png", dpi=150)
    plt.close(fig)

    # summary
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
    fmt = lambda d, **k: d.to_string(float_format=lambda x: f"{x:.3f}", **k)
    strong = {r: int((~kept & (ranks.min(axis=1) <= r)).sum()) for r in STRONG_RANKS}
    show = ["gene", "n_samples_top", "samples_top", "best_rank", "best_sample"] + [f"rank_{s}" for s in shorts]
    lines = ["NCD ratio analysis, step 1: gene filter", "",
             f"Rule: in the top {N} by mean count per tumor cell in >= {nr.MIN_SAMPLES} of 5 samples",
             f"Kept: {kept.sum()} of {len(genes)} genes  (in the top {N} of 3 samples: {(n_top == 3).sum()}, of 4: {(n_top == 4).sum()}, of 5: {(n_top == 5).sum()})",
             f"Measurable for the variation list (>= {nr.MIN_CELLS_GE2} cells with two or more transcripts in >= {nr.MIN_SAMPLES} samples): "
             f"{(kept & measurable).sum()} kept, {(kept & ~measurable).sum()} kept but not measurable, {(~kept & measurable).sum()} measurable but dropped",
             f"Dropped: {(~kept).sum()} genes  (in the top {N} of no sample: {(n_top == 0).sum()}, of 1: {(n_top == 1).sum()}, of 2: {(n_top == 2).sum()})",
             "Dropped genes that are strongly expressed in one or two samples: "
             + ", ".join(f"{v} with a best rank <= {r}" for r, v in strong.items()), "",
             "Kept genes per sample", fmt(per_sample, index=False), "",
             f"Genes shared by the per-sample top-{N} lists", fmt(overlap), "",
             f"Direction of the kept genes: median over samples of the mean r on detected cells is positive (more cytoplasmic) for "
             f"{(med_r[kept] > 0).sum()} genes and negative for {(med_r[kept] < 0).sum()}; the sign agrees in all five samples for "
             f"{np.isin(sign_agree[kept], [0, 5]).sum()} genes", "",
             "Dropped genes with the best single-sample ranks (first 40)", fmt(dropped[show].head(40), index=False), "",
             "Known gene groups", fmt(group_table, index=False), ""]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    assert (nr.load_filtered_genes() == genes[kept]).all()
    print(f"-> {out}")


if __name__ == "__main__":
    main()
