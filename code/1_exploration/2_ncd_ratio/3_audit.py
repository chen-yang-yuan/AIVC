"""NCD ratio analysis, step 3: audit of the candidate target gene lists (plans/ncd_ratio_exploration.md, section 4.5).

Compares the two variation lists of step 2 at the size of the primary list (ncd_ratio_utils.LIST_SIZE = 500; the
larger nested list LIST_SIZE_CHECK = 1000 is reported next to it): V-null (variance of r above the counting null)
and V-abs (absolute variance of r). List M (absolute centred mean score) and the 634 target genes of the cytoplasmic
analysis are reference columns only.
V-null and V-abs share most genes, so the gene sets looked at are: shared, only V-null, only V-abs.

Aspects:
  overlap            shared genes; where the genes of one list rank in the other
  gene content       curated programs and hallmark sets, surface and secreted genes (name-based annotation, see
                     SURFACE / SECRETED below), hypergeometric test against the filtered genes; direction; expression
  reproducibility    Spearman of the within-sample statistic between the samples; samples in which a listed gene
                     qualifies
  stability          cells of every sample split in two random halves, statistics and lists recomputed per half
                     (cells with two or more transcripts required: MIN_CELLS_GE2 / 2)
  reliability        transcripts of every pair with two or more transcripts split at random in two halves, r computed
                     in each half; per gene the correlation of the two half-scores across cells, and the same on the
                     count-preserving null (reference; not exactly 0 because the halves share the pair's total)
  expression         the lists rebuilt with the statistic ranked within deciles of the detection rate (expression
                     level held fixed): how much of each list survives

No embedding or clustering here. The choice of the list is made in the review and recorded in
ncd_ratio_utils.TARGET_LIST (2026-11-18: V-abs). The script writes that ranking as target_genes.csv / target_genes.npy (LIST_SIZE_POOL genes, ranked; the first LIST_SIZE
genes are the primary list, flag in_primary_list), read by ncd_ratio_utils.load_target_genes.

Run from this folder (AIVC-env, about 5 min):  python 3_audit.py [--quick]
Needs: output of steps 1 and 2 (1_filter/filtered_genes.npy, 2_selection/gene_stats.csv, gene_summary.csv).
Output: output/1_exploration/ncd_ratio/3_audit/ (emptied first; --quick writes to 3_audit_quick/)
  set_profile.csv       one row per gene set: size, direction, expression, statistics, surface / secreted share,
                        reliability, stability
  overlap.csv, groups.csv, reproducibility.csv, stability.csv, expression_dependence.csv
  gene_audit.csv        per gene: set membership, annotation, reliability, half-split ranks
  half_stats.csv        per sample, half and gene: the statistics of the cell split
  genes_shared.csv, genes_only_v_null.csv, genes_only_v_abs.csv     the gene sets with their statistics, ranked
  target_genes.csv, target_genes.npy     the chosen ranking (ncd_ratio_utils.TARGET_LIST), first LIST_SIZE_POOL genes
  target_genes_union.csv, .npy           union of the V-abs and V-null lists of LIST_SIZE (column "set"); check of step 4
  summary.txt, audit_overview.png
"""
import argparse
import os
import re
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
from scipy.stats import hypergeom, spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # shared subtype_utils.py
import ncd_ratio_utils as nr
import subtype_utils as su

N, N_CHECK = nr.LIST_SIZE, nr.LIST_SIZE_CHECK
N_NULL_REL = 5   # null draws for the reference of the transcript-split reliability
N_BINS = 10      # detection-rate bins of the expression-matched rebuild
LISTS = {"v_null": "V-null", "v_abs": "V-abs", "m": "M"}
COL = {"v_null": "#d95f02", "v_abs": "#1b9e77", "m": "#7570b3", "shared": "#666666"}

# Name-based annotation of genes whose protein is on the cell surface or secreted. A proxy: no atlas annotation is in
# the repo. Patterns are matched on the whole gene symbol.
SURFACE = re.compile(
    r"^(CD\d+[A-Z]?\d*|ITG[AB]\d+[A-Z]?|CDH\d+|CLDN\d+|EPH[AB]\d+|EFN[AB]\d+|ERBB\d|FGFR\d|PDGFR[AB]|TNFRSF\d+[A-Z]?|HLA-.+|"
    r"IL\d+R[A-Z]?\d?|CXCR\d|CCR\d+|NOTCH\d|DLL\d|JAG\d|FZD\d+|PLXN[A-D]\d|SEMA4[A-G]|TSPAN\d+|PTPR[A-Z]\d?|CEACAM\d+|"
    r"NECTIN\d|TGFBR\d|IFNGR\d|IFNAR\d|NTRK\d|LRP[1568]|NRP[12]|SDC[1-4]|GPC[1-6]|DSG\d|DSC\d|MUC1|MUC4|MUC16|"
    r"EGFR|MET|KIT|RET|ALK|AXL|MERTK|TYRO3|ROR[12]|DDR[12]|INSR|IGF1R|FLT[134]|KDR|TEK|EPCAM|MSLN|FOLR1|TACSTD2|PVR|"
    r"L1CAM|NCAM1|ICAM1|VCAM1|PECAM1|PTPRC|PDCD1|CTLA4|LAG3|HAVCR2|TIGIT|B2M|PROM1|ALCAM|MCAM|BCAM|F11R|OCLN|TFRC|"
    r"LDLR|SLC2A1|CA9|ENG|THY1|NT5E|LGR[456]|ADAM10|ADAM17|MME|ANPEP|DPP4|FAP)$")
SECRETED = re.compile(
    r"^(CXCL\d+|CCL\d+|IL\d+[AB]?|IFN[ABG]\d*|TNFSF\d+[A-Z]?|TGFB[123]|WNT\d+[AB]?|FGF\d+|VEGF[A-D]|PDGF[A-D]|IGFBP\d|IGF[12]|"
    r"COL\d+A\d|MMP\d+|TIMP\d|SERPIN[A-Z]\d+|LAM[ABC]\d|BMP\d+[AB]?|ANGPT\d|ANGPTL\d|CSF\d|THBS\d|SFRP\d|TFF\d|SCGB\d[A-D]\d|"
    r"C1Q[ABC]|C1[RS]|C[23]|CF[BDHI]|FG[ABG]|SAA[12]|LOXL?\d?|CCN[1-6]|"
    r"TNF|EGF|AREG|EREG|HBEGF|TGFA|NRG1|HGF|KITLG|LIF|OSM|SPP1|POSTN|FN1|SPARC|TNC|VCAN|LUM|DCN|BGN|FBN1|ELN|DKK1|"
    r"GDF15|INHBA|MUC2|MUC5AC|MUC5B|MUC6|LCN2|SLPI|PI3|WFDC2|MSMB|KLK[23]|CST3|LYZ|CLU|A2M|ALB|TF|CRP|HP|GRN)$")


def matrix(stats, col, shorts, kept):
    """genes x samples matrix of one column of a per-sample table."""
    return np.column_stack([stats[stats["sample"] == s].set_index("gene").loc[kept, col].values for s in shorts])


def pair_corr(cols, a, b, G, min_pairs):
    """Per gene, the Pearson correlation of a and b over its pairs (nan below min_pairs or without variance)."""
    cnt = np.bincount(cols, minlength=G).astype(float)
    s = lambda w: np.bincount(cols, weights=w, minlength=G)
    m = np.maximum(cnt, 1)
    cov = s(a * b) / m - s(a) / m * s(b) / m
    va, vb = s(a * a) / m - (s(a) / m) ** 2, s(b * b) / m - (s(b) / m) ** 2
    ok = (cnt >= min_pairs) & (va > 1e-12) & (vb > 1e-12)
    return np.where(ok, cov / np.sqrt(np.maximum(va * vb, 1e-24)), np.nan), cnt


def sample_audit(sample, kept, quick, rng):
    """Cell split (statistics per half) and transcript split (reliability per gene) of one sample."""
    Xc, Xn, cells, genes = nr.load_counts(sample, nr.QUICK_CELLS if quick else None, rng)
    cols_keep = np.isin(genes, kept)
    assert (genes[cols_keep] == kept).all(), "filtered genes are not in panel order"
    rows, cols, c, n = nr.pairs(Xc[:, cols_keep].tocsr(), Xn[:, cols_keep].tocsr())
    del Xc, Xn
    G = len(kept)

    # cells split in two random halves
    half = rng.random(len(cells)) < 0.5
    halves = []
    for h in (0, 1):
        m = half[rows] == bool(h)
        st = nr.pair_stats(cols[m], c[m], n[m], G, int((half == bool(h)).sum()), rng, min_ge2=nr.MIN_CELLS_GE2 // 2)
        halves.append(pd.DataFrame({"sample": nr.short(sample), "half": "AB"[h], "gene": kept,
                                    **{k: st[k] for k in ["cells_detected", "cells_ge2", "mean_r_centred", "var_r", "excess", "z"]}}))

    # transcripts of every pair with two or more transcripts split in two random halves
    t = c + n
    p = np.bincount(cols, weights=c, minlength=G) / np.maximum(np.bincount(cols, weights=t, minlength=G), 1)
    m2 = t >= 2
    cols2, c2, n2 = cols[m2], c[m2], n[m2]

    def split_corr(cc, nn):
        ca, na = rng.binomial(cc, 0.5), rng.binomial(nn, 0.5)
        return pair_corr(cols2, nr.ratio_score(ca, na), nr.ratio_score(cc - ca, nn - na), G, nr.MIN_CELLS_GE2)

    rel, cnt = split_corr(c2, n2)
    null = np.nanmean([split_corr(*nr.null_counts(c2, n2, p[cols2], rng))[0] for _ in range(N_NULL_REL)], axis=0)
    reliability = pd.DataFrame({"sample": nr.short(sample), "gene": kept, "pairs_ge2": cnt, "reliability": rel, "reliability_null": null})
    return pd.concat(halves, ignore_index=True), reliability


def binned_ranks(values, eligible, detection):
    """Within-sample ranks of the genes when the statistic is ranked within N_BINS bins of the detection rate (rank
    by the within-bin percentile), so that the expression level does not decide the rank. Genes with the same
    percentile in different bins keep the panel order."""
    pct = np.full(values.shape, np.nan)
    for j in range(values.shape[1]):
        e = eligible[:, j] & np.isfinite(values[:, j])
        d = pd.DataFrame({"v": values[e, j], "bin": pd.qcut(detection[e, j], N_BINS, labels=False, duplicates="drop")})
        pct[e, j] = d.groupby("bin")["v"].rank(pct=True).values
    return nr.ranks_desc(pct, eligible)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help=f"{nr.QUICK_CELLS} cells per sample, output in 3_audit_quick/")
    args = ap.parse_args()
    t0 = time.time()
    out = nr.step_dir("3_audit", args.quick)
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    rng = np.random.default_rng(nr.SEED)
    kept = nr.load_filtered_genes()
    G = len(kept)
    shorts = [nr.short(s) for s in nr.SAMPLES]
    step2 = nr.step_dir("2_selection")
    summary = pd.read_csv(step2 + "gene_summary.csv")
    stats = pd.read_csv(step2 + "gene_stats.csv")
    assert (summary["gene"].values == kept).all()
    mat = lambda col: matrix(stats, col, shorts, kept)
    excess, z, ge2, var, cen, det = mat("excess"), mat("z"), mat("cells_ge2"), mat("var_r"), mat("mean_r_centred"), mat("detection_rate")
    passing = (ge2 >= nr.MIN_CELLS_GE2) & (z >= nr.Z_MIN)
    everything = np.ones_like(passing)
    stat = {"v_null": excess, "v_abs": var, "m": np.abs(cen)}
    elig = {"v_null": passing, "v_abs": everything, "m": everything}
    entry, direction_m = nr.entry_ranks(excess, z, ge2, var, cen)
    for k in LISTS:   # the ranking of step 2 is reproduced from its tables
        assert np.allclose(np.where(np.isfinite(entry[k]), entry[k], -1), summary[f"entry_K_{k}"].replace(np.inf, -1).fillna(-1).values)
    rank = {k: summary[f"rank_{k}"].values for k in LISTS}

    # ---- gene sets
    in_list = {k: rank[k] <= N for k in LISTS}
    sets = {"shared": in_list["v_null"] & in_list["v_abs"], "only V-null": in_list["v_null"] & ~in_list["v_abs"],
            "only V-abs": in_list["v_abs"] & ~in_list["v_null"], "V-null": in_list["v_null"], "V-abs": in_list["v_abs"],
            f"V-null first {N_CHECK}": rank["v_null"] <= N_CHECK, f"V-abs first {N_CHECK}": rank["v_abs"] <= N_CHECK,
            "M (reference)": in_list["m"], "634 target genes (reference)": summary["target_634"].values.astype(bool),
            "all filtered genes": np.ones(G, bool)}
    assert sets["shared"].sum() + sets["only V-null"].sum() == N and sets["shared"].sum() + sets["only V-abs"].sum() == N
    surface = np.array([bool(SURFACE.match(g)) for g in kept])
    secreted = np.array([bool(SECRETED.match(g)) for g in kept]) & ~surface

    # ---- overlap
    rows = []
    for size in sorted((N, N_CHECK)):
        a, b, m = rank["v_null"] <= size, rank["v_abs"] <= size, rank["m"] <= size
        rows.append({"size": size, "shared V-null / V-abs": int((a & b).sum()), "only V-null": int((a & ~b).sum()), "only V-abs": int((b & ~a).sum()),
                     "V-null in M": int((a & m).sum()), "V-abs in M": int((b & m).sum()),
                     "only V-null: median rank in V-abs": float(np.median(rank["v_abs"][a & ~b])),
                     "only V-abs: median rank in V-null": float(np.median(np.where(np.isfinite(rank["v_null"][b & ~a]), rank["v_null"][b & ~a], G))),
                     "only V-abs: not above the null in 3 samples": int((~np.isfinite(rank["v_null"][b & ~a])).sum())})
    overlap = pd.DataFrame(rows)
    overlap.to_csv(out + "overlap.csv", index=False)

    # ---- gene content: programs, surface, secreted
    groups = {"surface (name-based)": kept[surface], "secreted (name-based)": kept[secreted]}
    groups.update({f"curated: {k}": v for k, v in su.CURATED_PROGRAMS.items()})
    groups.update({f"hallmark: {k.replace('HALLMARK_', '').replace('_', ' ').capitalize()}": v
                   for k, v in su.read_gmt(nr.UTILS_PATH + "hallmark_pathways_filtered.gmt").items()})
    test_sets = ["V-null", "V-abs", "M (reference)", "shared", "only V-null", "only V-abs"]
    rows = []
    for name, g in groups.items():
        m = np.isin(kept, g)
        row = {"group": name, "filtered_genes": int(m.sum())}
        for s in test_sets:
            k = int((m & sets[s]).sum())
            row[s] = k
            row[f"p {s}"] = hypergeom.sf(k - 1, G, int(m.sum()), int(sets[s].sum())) if m.sum() else np.nan
        row["genes only V-null"] = " ".join(kept[m & sets["only V-null"]])
        row["genes only V-abs"] = " ".join(kept[m & sets["only V-abs"]])
        rows.append(row)
    group_table = pd.DataFrame(rows)
    group_table.to_csv(out + "groups.csv", index=False)

    # ---- reproducibility across samples
    rows = []
    for k in LISTS:
        rho = []
        for i in range(5):
            for j in range(i + 1, 5):
                e = elig[k][:, i] & elig[k][:, j]
                rho.append(spearmanr(stat[k][e, i], stat[k][e, j])[0])
        within = nr.ranks_desc(stat[k], elig[k])
        if k == "m":
            within = np.where(np.sign(cen) == direction_m[:, None], within, np.inf)
        K = np.sort(entry[k])[N - 1]
        q = (within[in_list[k]] <= K).sum(axis=1)
        rows.append({"list": LISTS[k], "statistic": {"v_null": "excess", "v_abs": "var_r", "m": "|centred mean r|"}[k],
                     "spearman between samples, mean": np.mean(rho), "min": np.min(rho), "max": np.max(rho),
                     "K needed": K, "listed genes in the top K of 5 samples": (q == 5).mean(), "of 4": (q == 4).mean(), "of 3": (q == 3).mean()})
        assert (q >= nr.MIN_SAMPLES).all(), "a listed gene qualifies in fewer than MIN_SAMPLES samples"
    repro = pd.DataFrame(rows)
    repro.to_csv(out + "reproducibility.csv", index=False)

    # ---- dependence on the expression level: lists rebuilt within detection-rate bins
    rows = []
    matched = {}
    for k in LISTS:
        rb = binned_ranks(stat[k], elig[k], det)
        eb = nr.third_best(rb)
        if k == "m":   # as in entry_ranks: the qualifying samples must agree on the direction
            eb = np.minimum(*[nr.third_best(np.where(np.sign(cen) == d, rb, np.inf)) for d in (1, -1)])
        matched[k] = nr.top_genes(eb, -np.nanmedian(stat[k], axis=1), N)
        rows.append({"list": LISTS[k],
                     "spearman statistic vs detection rate (mean over samples)": np.mean([spearmanr(stat[k][elig[k][:, j], j], det[elig[k][:, j], j])[0] for j in range(5)]),
                     "spearman entry K vs detection rate": spearmanr(np.where(np.isfinite(entry[k]), entry[k], G + 1), summary["median_detection_rate"])[0],
                     "median detection rate of the list": float(summary["median_detection_rate"][in_list[k]].median()),
                     "median detection rate, expression-matched list": float(summary["median_detection_rate"][matched[k]].median()),
                     f"genes kept by the expression-matched list (of {N})": int((matched[k] & in_list[k]).sum())})
    expr = pd.DataFrame(rows)
    expr.to_csv(out + "expression_dependence.csv", index=False)
    matched_overlap = int((matched["v_null"] & matched["v_abs"]).sum())

    # ---- new computations on the counts: cell split and transcript split
    half_tables, rel_tables = [], []
    for sample in nr.SAMPLES:
        h, r = sample_audit(sample, kept, args.quick, rng)
        half_tables.append(h); rel_tables.append(r)
        print(f"[{nr.short(sample)}] cell split and transcript split done ({time.time() - t0:.0f} s)", flush=True)
    half_stats = pd.concat(half_tables, ignore_index=True)
    half_stats.to_csv(out + "half_stats.csv", index=False)
    rel_all = pd.concat(rel_tables, ignore_index=True)
    rel = np.nanmedian(matrix(rel_all, "reliability", shorts, kept), axis=1)
    rel_null = np.nanmedian(matrix(rel_all, "reliability_null", shorts, kept), axis=1)

    min_half = nr.MIN_CELLS_GE2 // 2
    hs, half_in, half_rank = {}, {}, {}
    for h in "AB":
        d = half_stats[half_stats["half"] == h]
        hm = lambda col: matrix(d, col, shorts, kept)
        hs[h] = {"v_null": hm("excess"), "v_abs": hm("var_r"), "m": np.abs(hm("mean_r_centred")),
                 "elig": (hm("cells_ge2") >= min_half) & (hm("z") >= nr.Z_MIN)}
        e, _ = nr.entry_ranks(hm("excess"), hm("z"), hm("cells_ge2"), hm("var_r"), hm("mean_r_centred"), min_ge2=min_half)
        for k in LISTS:
            tie = -np.nanmedian(hs[h][k], axis=1)
            half_in[h, k, N], half_in[h, k, N_CHECK] = nr.top_genes(e[k], tie, N), nr.top_genes(e[k], tie, N_CHECK)
            half_rank[h, k] = e[k]
    rows = []
    for k in LISTS:
        rho = []
        for j in range(5):
            e = (hs["A"]["elig"][:, j] & hs["B"]["elig"][:, j]) if k == "v_null" else np.ones(G, bool)
            rho.append(spearmanr(hs["A"][k][e, j], hs["B"][k][e, j])[0])
        row = {"list": LISTS[k], "spearman between the halves, mean over samples": np.mean(rho), "min": np.min(rho)}
        for size in (N, N_CHECK):
            a, b = half_in["A", k, size], half_in["B", k, size]
            full = rank[k] <= size
            row[f"genes in both half lists (size {size})"] = int((a & b).sum())
            row[f"full list found in both halves (size {size})"] = int((a & b & full).sum())
        rows.append(row)
    stability = pd.DataFrame(rows)
    stability.to_csv(out + "stability.csv", index=False)
    both_halves = {k: half_in["A", k, N] & half_in["B", k, N] for k in LISTS}

    # ---- per gene and per set
    gene_audit = summary[["gene", "direction", "median_detection_rate", "median_cells_ge2", "median_cytoplasmic_fraction", "median_mean_r_centred",
                          "median_var_r", "median_excess", "median_z", "median_excess_t2", "n_samples_z_pass", "target_634",
                          "rank_v_null", "rank_v_abs", "rank_m"]].copy()
    gene_audit["set"] = np.select([sets["shared"], sets["only V-null"], sets["only V-abs"]], ["shared", "only V-null", "only V-abs"], "not listed")
    gene_audit["surface"], gene_audit["secreted"] = surface, secreted
    gene_audit["reliability"], gene_audit["reliability_null"] = rel, rel_null
    for k in ("v_null", "v_abs"):
        gene_audit[f"in_both_half_lists_{k}"] = both_halves[k]
        gene_audit[f"in_expression_matched_{k}"] = matched[k]
    gene_audit.to_csv(out + "gene_audit.csv", index=False)
    gene_audit[sets["shared"]].sort_values("rank_v_null").to_csv(out + "genes_shared.csv", index=False)
    gene_audit[sets["only V-null"]].sort_values("rank_v_null").to_csv(out + "genes_only_v_null.csv", index=False)
    gene_audit[sets["only V-abs"]].sort_values("rank_v_abs").to_csv(out + "genes_only_v_abs.csv", index=False)

    # ---- the target list: the chosen list (nr.TARGET_LIST), in the order of its ranking
    target = gene_audit[rank[nr.TARGET_LIST] <= nr.LIST_SIZE_POOL].sort_values(f"rank_{nr.TARGET_LIST}").copy()
    target.insert(1, "in_primary_list", target[f"rank_{nr.TARGET_LIST}"] <= nr.LIST_SIZE)
    assert len(target) == nr.LIST_SIZE_POOL and target["in_primary_list"].sum() == nr.LIST_SIZE
    target.to_csv(out + "target_genes.csv", index=False)
    np.save(out + "target_genes.npy", target["gene"].values.astype(str))
    # union of the two variation lists at the primary size (a check of step 4): the chosen list in its order, then
    # the genes only in the other list in the order of that list
    assert nr.TARGET_LIST in ("v_abs", "v_null")
    other = "v_null" if nr.TARGET_LIST == "v_abs" else "v_abs"
    first = gene_audit[in_list[nr.TARGET_LIST]].sort_values(f"rank_{nr.TARGET_LIST}")
    second = gene_audit[in_list[other] & ~in_list[nr.TARGET_LIST]].sort_values(f"rank_{other}")
    union = pd.concat([first, second], ignore_index=True)
    assert len(union) == 2 * N - sets["shared"].sum() and union["gene"].is_unique
    union.to_csv(out + "target_genes_union.csv", index=False)
    np.save(out + "target_genes_union.npy", union["gene"].values.astype(str))

    rows = []
    for name, m in sets.items():
        own = "v_null" if "V-null" in name else "v_abs" if "V-abs" in name else None
        rows.append({
            "set": name, "genes": int(m.sum()), "cytoplasmic": int((summary["direction"][m] == "cytoplasmic").sum()),
            "nuclear": int((summary["direction"][m] == "nuclear").sum()),
            "median detection rate": summary["median_detection_rate"][m].median(),
            "median cytoplasmic fraction": summary["median_cytoplasmic_fraction"][m].median(),
            "median var r": summary["median_var_r"][m].median(), "median excess": summary["median_excess"][m].median(),
            "median excess at two transcripts": summary["median_excess_t2"][m].median(),
            "above the null in 5 samples": (summary["n_samples_z_pass"][m] == 5).mean(),
            "in the 634 target genes": int(summary["target_634"][m].sum()),
            "surface": int((surface & m).sum()), "secreted": int((secreted & m).sum()),
            "median reliability": np.nanmedian(rel[m]), "median reliability of the null": np.nanmedian(rel_null[m]),
            "genes with reliability": int(np.isfinite(rel[m]).sum()),
            "found in both half lists": both_halves[own][m].mean() if own and "first" not in name else np.nan,
            "kept by the expression-matched list": matched[own][m].mean() if own and "first" not in name else np.nan})
    profile = pd.DataFrame(rows)
    profile.to_csv(out + "set_profile.csv", index=False)

    # ---- figure
    fig, ax = plt.subplots(2, 4, figsize=(24, 11))
    three = ["shared", "only V-null", "only V-abs"]
    three_col = [COL["shared"], COL["v_null"], COL["v_abs"]]
    a = ax[0, 0]
    x, y = rank["v_null"].copy(), rank["v_abs"].copy()
    x[~np.isfinite(x)] = G
    a.scatter(x, y, s=3, color="#dddddd")
    for s, c_ in zip(three, three_col):
        a.scatter(x[sets[s]], y[sets[s]], s=5, color=c_, label=f"{s} ({sets[s].sum()})")
    a.axvline(N, color="k", lw=0.6, ls=":"); a.axhline(N, color="k", lw=0.6, ls=":")
    a.set(xlabel="rank in V-null (genes not above the null at the right edge)", ylabel="rank in V-abs", title=f"The two rankings (lists of {N})")
    a.legend(frameon=False, markerscale=3)
    a = ax[0, 1]
    data = [np.log10(summary["median_detection_rate"][sets[s]]) for s in three + ["M (reference)"]]
    bp = a.boxplot(data, labels=three + ["M (ref.)"], patch_artist=True, showfliers=False)
    for patch, c_ in zip(bp["boxes"], three_col + [COL["m"]]):
        patch.set_facecolor(c_); patch.set_alpha(0.6)
    a.set(ylabel="log10 detection rate (median over samples)", title="Expression level of the gene sets")
    a = ax[0, 2]
    bins = np.linspace(0, 1, 41)
    for s, c_ in zip(three + ["M (reference)"], three_col + [COL["m"]]):
        a.hist(summary["median_cytoplasmic_fraction"][sets[s]], bins=bins, histtype="step", lw=1.6, color=c_, density=True, label=s)
    a.set(xlabel="cytoplasmic fraction of the gene (median over samples)", ylabel="density", title="Balance of the gene sets")
    a.legend(frameon=False)
    a = ax[0, 3]
    show = group_table[(group_table["filtered_genes"] >= 8) & (group_table["group"].str.startswith(("surface", "secreted", "curated")))].copy()
    yy = np.arange(len(show))[::-1]
    for i, (s, c_) in enumerate(zip(["V-null", "V-abs", "M (reference)"], [COL["v_null"], COL["v_abs"], COL["m"]])):
        a.barh(yy + (1 - i) * 0.27, show[s] / (show["filtered_genes"] * sets[s].sum() / G), height=0.27, color=c_, label=s)
    a.axvline(1, color="k", lw=0.6, ls="--")
    a.set_yticks(yy); a.set_yticklabels([f"{g} ({n_})" for g, n_ in zip(show["group"], show["filtered_genes"])], fontsize=8)
    a.set(xlabel="genes in the list / expected by chance", title="Gene content (filtered genes in brackets)")
    a.legend(frameon=False)
    a = ax[1, 0]
    w = 0.38
    xx = np.arange(3)
    a.bar(xx - w / 2, repro["spearman between samples, mean"], w, color="#999999", label="between samples (mean of 10 pairs)")
    a.bar(xx + w / 2, stability["spearman between the halves, mean over samples"], w, color="#444444", label="between cell halves (mean of 5 samples)")
    a.set_xticks(xx); a.set_xticklabels([f"{LISTS[k]}\n{s}" for k, s in zip(LISTS, repro["statistic"])])
    a.set(ylabel="Spearman over genes", title="Reproducibility of the selection statistic", ylim=(0, 1.25))
    a.legend(frameon=False, loc="upper left")
    a = ax[1, 1]
    a.bar(xx - w / 2, [stability[f"full list found in both halves (size {N})"][i] / N for i in range(3)], w, color="#999999", label="found in both cell-half lists")
    a.bar(xx + w / 2, expr[f"genes kept by the expression-matched list (of {N})"] / N, w, color="#444444", label="kept when expression is held fixed")
    a.set_xticks(xx); a.set_xticklabels(list(LISTS.values()))
    a.set(ylabel=f"share of the list of {N}", title="Stability of the lists", ylim=(0, 1.25))
    a.legend(frameon=False, loc="upper left")
    a = ax[1, 2]
    data = [rel[sets[s]][np.isfinite(rel[sets[s]])] for s in three + ["M (reference)", "all filtered genes"]]
    bp = a.boxplot(data, labels=three + ["M (ref.)", "all"], patch_artist=True, showfliers=False)
    for patch, c_ in zip(bp["boxes"], three_col + [COL["m"], "#dddddd"]):
        patch.set_facecolor(c_); patch.set_alpha(0.6)
    for i, s in enumerate(three + ["M (reference)", "all filtered genes"]):
        a.plot([i + 0.6, i + 1.4], [np.nanmedian(rel_null[sets[s]])] * 2, color="r", lw=1.5, label="median under the null" if i == 0 else None)
    a.axhline(0, color="k", lw=0.5)
    a.set(ylabel="correlation of the two half-scores across cells", title="Reliability of the per-cell score (transcript split)")
    a.legend(frameon=False)
    a = ax[1, 3]
    a.scatter(np.log10(summary["median_detection_rate"]), rel - rel_null, s=3, color="#dddddd")
    for s, c_ in zip(three, three_col):
        a.scatter(np.log10(summary["median_detection_rate"][sets[s]]), (rel - rel_null)[sets[s]], s=5, color=c_, label=s)
    a.axhline(0, color="k", lw=0.5)
    a.set(xlabel="log10 detection rate (median over samples)", ylabel="reliability minus its null", title="Reliability against expression")
    a.legend(frameon=False, markerscale=3, loc="lower left")
    fig.tight_layout()
    fig.savefig(out + "audit_overview.png", dpi=150)
    plt.close(fig)

    # ---- summary
    pd.set_option("display.width", 320); pd.set_option("display.max_colwidth", 200); pd.set_option("display.max_columns", 50)
    num = lambda v: f"{v:.0f}" if float(v).is_integer() else f"{v:.3g}" if abs(v) < 0.001 else f"{v:.3f}"
    fmt = lambda d, **k: d.to_string(float_format=num, **k)
    gcols = ["group", "filtered_genes"] + test_sets + ["p V-null", "p V-abs", "p only V-null", "p only V-abs"]
    gshow = group_table[group_table["filtered_genes"] >= 8][gcols]
    names = lambda m, by: " ".join(summary.loc[m].sort_values(by)["gene"].head(80))
    lines = [f"NCD ratio analysis, step 3: audit of the candidate lists{' (QUICK RUN)' if args.quick else ''}", "",
             f"Primary list size {N} (larger nested check: {N_CHECK}). Candidates: V-null and V-abs. M and the 634 target genes are references.",
             f"Target list (decision after the review, ncd_ratio_utils.TARGET_LIST): {LISTS[nr.TARGET_LIST]} -> target_genes.csv, target_genes.npy ({nr.LIST_SIZE_POOL} ranked genes; the first {nr.LIST_SIZE} are the primary list)",
             f"Surface / secreted: name-based annotation, {surface.sum()} / {secreted.sum()} of the {G} filtered genes.", "",
             "Overlap", fmt(overlap.T, header=False), "",
             "Gene sets", fmt(profile.set_index("set").T), "",
             "Gene groups: genes in the sets, and hypergeometric p of enrichment against the filtered genes (groups with >= 8 filtered genes)",
             fmt(gshow, index=False), "",
             "Reproducibility across samples", fmt(repro, index=False), "",
             f"Stability under a split of the cells (two halves per sample; cells with two or more transcripts required: {min_half})",
             fmt(stability.T, header=False), "",
             "Dependence on the expression level", fmt(expr.T, header=False),
             f"V-null and V-abs after expression matching share {matched_overlap} of {N} genes", "",
             "Surface genes, only V-null: " + " ".join(kept[surface & sets["only V-null"]]),
             "Surface genes, only V-abs:  " + " ".join(kept[surface & sets["only V-abs"]]),
             "Surface genes, shared:      " + " ".join(kept[surface & sets["shared"]]),
             "Secreted genes, only V-null: " + " ".join(kept[secreted & sets["only V-null"]]),
             "Secreted genes, only V-abs:  " + " ".join(kept[secreted & sets["only V-abs"]]),
             "Secreted genes, shared:      " + " ".join(kept[secreted & sets["shared"]]), "",
             "Only V-null, first 80 by its rank:", names(sets["only V-null"], "rank_v_null"), "",
             "Only V-abs, first 80 by its rank:", names(sets["only V-abs"], "rank_v_abs"), ""]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    print(f"done in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
