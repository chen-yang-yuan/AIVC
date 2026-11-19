"""Shared helpers of the NCD ratio analysis (plans/ncd_ratio_exploration.md).

Score of gene g in cell i, c = cytoplasmic count, n = nuclear count:

    r = log2((c + 1) / (n + 1))

r = 0 for an empty pair and for a balanced pair, r > 0 more cytoplasmic, r < 0 more nuclear.

Contents: the constants of the analysis (gene filter, null, list size, chosen list), loading of the counts and of
the outputs of earlier steps, the score, the count-preserving null, the selection statistics (pair_stats) and the
ranking rules of the candidate lists (ranks_desc, third_best, entry_ranks, top_genes), and for step 4 the embedding
constants and the per-sample correction of the score (correct_within_sample).

The compartment matrices of processed_data/ hold the tumor cells (rows, order of
cell_ids.npy) x the 5001 panel genes (columns, order of gene_ids.npy) with all decoded transcripts (no Q20 filter).

Cells of the analysis (decided 2026-11-18, plan document section 10.8): the tumor cells with at least one nuclear
transcript over all panel genes. A cell with zero nuclear transcripts (no segmented nucleus, or a small nucleus with
no transcript in it) has n = 0 for every gene, so every detected gene would score positive. load_counts drops them.

Genes of the analysis (decided 2026-11-18): the protein-coding genes of the panel (Ensembl biotype protein_coding in
data/_utils/gene_biotype.csv, written by gene_biotype.py; 4972 of the 5001 panel genes). load_counts drops the
other columns, so every rule of the analysis (expression rank, top N, the lists) is stated among coding genes.
"""
import os

import numpy as np
from scipy import sparse

SAMPLES = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate"]
OUT_ROOT = "../../../output/1_exploration/ncd_ratio/"
UTILS_PATH = "../../../data/_utils/"
SEED = 0
QUICK_CELLS = 20000
N_FILTER = 3000   # gene filter: genes in the top N_FILTER by expression in tumor cells in >= MIN_SAMPLES samples (step 0)
MIN_SAMPLES = 3    # a gene must qualify in at least this many of the five samples
MIN_CELLS_GE2 = 200   # measurable gene: cells with two or more transcripts required per sample (in >= MIN_SAMPLES samples)
N_NULL = 20       # draws of the count-preserving null (step 2)
Z_MIN = 3.0       # evidence required for the variation list: (observed - null mean) / null SD of the variance of r
LIST_SIZE = 500          # size of the primary target gene list: the first genes of the chosen ranking (decided 2026-11-18)
LIST_SIZE_CHECK = 1000   # nested larger list, a check of step 4 (250 is a possible further check, not planned); the second
                         # check of step 4 is the union of the V-abs and V-null lists of LIST_SIZE (load_target_genes)
LIST_SIZE_POOL = 1000    # genes of each ranking written to the list files
TARGET_LIST = "v_abs"   # the list chosen in the audit of step 3 (V-abs: absolute variance of r); decided 2026-11-18
# step 4 (4_embedding.py): pooled embedding and clustering. Parameters as in the cytoplasmic analysis
# (1_cytoplasmic/2_embedding.py, subtype_utils.RES_GRID) so that the two can be read side by side.
LISTS = ["primary", "larger", "union"]     # gene lists of load_target_genes
CORRECTIONS = ["zscore", "depth"]           # per-sample correction of r before pooling (correct_within_sample)
RUNS = ["data", "null", "half_a", "half_b"]   # observed counts, count-preserving null, the two transcript halves
N_COMPS_PCA = 50
N_PCS = 30          # PCs used for the neighbours
N_NEIGHBORS = 30
SCALE_MAX = 10.0    # corrected scores are clipped at +- SCALE_MAX
RES_GRID = [0.5, 0.8, 1.0, 1.2, 1.5]   # clustering resolutions
L2 = float(np.log2(3.0))   # |r| of a pair with two transcripts in one compartment
BIG = 1 << 21   # packs (nuclear, cytoplasmic) counts of one pair into one integer


def short(sample):
    return sample.replace("Xenium_5K_", "")


def step_dir(name, quick=False):
    """Output folder of a step; a quick run writes to <name>_quick/."""
    return os.path.join(OUT_ROOT, name + ("_quick" if quick else "")) + "/"


def coding_genes():
    """Names of the protein-coding panel genes (Ensembl biotype protein_coding; data/_utils/gene_biotype.csv)."""
    import pandas as pd
    t = pd.read_csv(UTILS_PATH + "gene_biotype.csv")
    return t.loc[t["biotype"] == "protein_coding", "gene"].astype(str).values


def load_counts(sample, quick_cells=None, rng=None, all_cells=False, all_genes=False):
    """Cytoplasmic and nuclear count matrices (csr, cells x genes), cell ids and gene names of one sample.
    Cells: the tumor cells with at least one nuclear transcript (see the module docstring); all_cells=True returns
    every tumor cell. Genes: the protein-coding panel genes, in panel order; all_genes=True returns the whole panel.
    The two flags are used by check_upstream.py only. The cell rule counts the transcripts of all panel genes.
    quick_cells: random subsample of the cells for development runs, drawn after the removal."""
    processed = f"../../../data/{sample}/processed_data/"
    genes = np.load(processed + "gene_ids.npy", allow_pickle=True).astype(str)
    cells = np.load(processed + "cell_ids.npy", allow_pickle=True).astype(str)
    Xc = sparse.load_npz(processed + "cytoplasmic_expression_matrix.npz").tocsr()
    Xn = sparse.load_npz(processed + "nuclear_expression_matrix.npz").tocsr()
    assert Xc.shape == Xn.shape == (len(cells), len(genes)), "matrices do not match tumor cells x panel"
    if not all_cells:
        keep = np.asarray(Xn.sum(axis=1)).ravel() > 0
        cells, Xc, Xn = cells[keep], Xc[keep], Xn[keep]
    if not all_genes:
        coding = np.isin(genes, coding_genes())
        genes, Xc, Xn = genes[coding], Xc[:, coding].tocsr(), Xn[:, coding].tocsr()
    if quick_cells and len(cells) > quick_cells:
        pick = np.sort(rng.choice(len(cells), quick_cells, replace=False))
        cells, Xc, Xn = cells[pick], Xc[pick], Xn[pick]
    return Xc, Xn, cells, genes


def pairs(Xc, Xn):
    """Non-empty (cell, gene) pairs of one sample: row index, column index, cytoplasmic count c, nuclear count n."""
    assert max(Xc.max(), Xn.max()) < BIG
    P = (Xc.astype(np.int64) + Xn.astype(np.int64) * BIG).tocoo()
    keep = P.data > 0
    return P.row[keep], P.col[keep], P.data[keep] % BIG, P.data[keep] // BIG


def ratio_score(c, n):
    """The NCD ratio score log2((c + 1) / (n + 1))."""
    return np.log2((c + 1.0) / (n + 1.0))


def load_filtered_genes():
    """The genes kept by the gene filter (step 1, 1_filter.py), in panel order."""
    return np.load(step_dir("1_filter") + "filtered_genes.npy", allow_pickle=True).astype(str)


def load_target_genes(which="primary"):
    """A target gene list written by step 3 (3_audit.py), as gene names. Decided 2026-11-18:
      "primary"  the first LIST_SIZE (500) genes of the chosen ranking TARGET_LIST (V-abs); the target of the analysis
      "larger"   the first LIST_SIZE_CHECK (1000) genes of the same ranking; nested check of step 4
      "union"    the primary list plus the first LIST_SIZE genes of the V-null ranking (630 genes); check of step 4
                 of whether the choice between V-abs and V-null matters
    An integer returns that many first genes of the chosen ranking (at most LIST_SIZE_POOL)."""
    if which == "union":
        return np.load(step_dir("3_audit") + "target_genes_union.npy", allow_pickle=True).astype(str)
    size = {"primary": LIST_SIZE, "larger": LIST_SIZE_CHECK}.get(which, which)
    assert isinstance(size, int) and size <= LIST_SIZE_POOL, which
    return np.load(step_dir("3_audit") + "target_genes.npy", allow_pickle=True).astype(str)[:size]


def null_counts(c, n, p, rng):
    """Count-preserving null of the pairs of one sample: the c + n transcripts of every pair are reassigned to
    the compartments with probability p of the cytoplasm (p: per pair, the cytoplasmic fraction of its gene in the
    sample). Returns the null cytoplasmic and nuclear counts."""
    t = c + n
    c0 = rng.binomial(t, p)
    return c0, t - c0


def pair_stats(cols, c, n, G, n_cells, rng, min_ge2=MIN_CELLS_GE2):
    """Selection statistics of the G genes from the non-empty pairs of one set of cells (cols: gene index of the pair,
    c, n: its counts), on the cells in which the gene is detected. Columns as in gene_stats.csv of step 2
    (2_selection.py documents them). n_cells: number of cells of the set; min_ge2: pairs with exactly two transcripts
    required for excess_t2."""
    t = c + n
    count = lambda w=None, m=None: np.bincount(cols if m is None else cols[m], weights=None if w is None else (w if m is None else w[m]), minlength=G)
    det = count()
    p = count(c) / np.maximum(count(t), 1)

    def moments(cc, nn):
        r = ratio_score(cc, nn)
        mean = np.bincount(cols, weights=r, minlength=G) / np.maximum(det, 1)
        return mean, np.bincount(cols, weights=r * r, minlength=G) / np.maximum(det, 1) - mean ** 2

    mean_r, var_r = moments(c, n)
    null_mean, null_var = [], []
    pp = p[cols]
    for _ in range(N_NULL):
        m0, v0 = moments(*null_counts(c, n, pp, rng))
        null_mean.append(m0); null_var.append(v0)
    null_mean, null_var = np.array(null_mean), np.array(null_var)

    # pairs with exactly two transcripts: observed variance of r against the exact binomial variance
    t2 = t == 2
    n2 = count(m=t2)
    both_c, both_n = count(m=t2 & (c == 2)), count(m=t2 & (n == 2))
    obs2 = L2 ** 2 * (both_c + both_n) / np.maximum(n2, 1) - (L2 * (both_c - both_n) / np.maximum(n2, 1)) ** 2
    p2 = (2 * both_c + (n2 - both_c - both_n)) / np.maximum(2 * n2, 1)
    exp2 = L2 ** 2 * (p2 ** 2 + (1 - p2) ** 2) - (L2 * (p2 ** 2 - (1 - p2) ** 2)) ** 2
    n1 = count(m=t == 1)
    out = {
        "cells_detected": det, "detection_rate": det / n_cells,
        "cells_ge2": det - n1, "cytoplasmic_fraction": p, "mean_r": mean_r, "mean_r_centred": mean_r - np.median(mean_r),
        "var_r": var_r, "null_var_mean": null_var.mean(axis=0), "null_var_sd": null_var.std(axis=0, ddof=1),
        "null_mean_r": null_mean.mean(axis=0), "pairs_t2": n2,
        "excess_t2": np.where((n2 >= min_ge2) & (exp2 > 0), obs2 / np.maximum(exp2, 1e-12), np.nan)}
    out["excess"] = out["var_r"] / out["null_var_mean"]
    out["z"] = (out["var_r"] - out["null_var_mean"]) / out["null_var_sd"]
    return out


def correct_within_sample(R, log_depth, correction):
    """Per-sample correction of the score matrix R (cells of ONE sample x genes, dense; empty pairs are 0) before
    pooling. Returns the corrected matrix (float32) and per gene: mean and SD of r, slope on depth, SD after the
    depth regression.
      "zscore"  z-score per gene over the cells of the sample
      "depth"   first the residual of a per-gene OLS of r on an intercept and log_depth (log10 of the cell's total
                transcripts + 1; the formula of regress_out_covariate in 1_cytoplasmic/2_embedding.py), then z-score
    The z-score uses the SD with ddof = 1 and is clipped at +- SCALE_MAX. A gene without variance in the sample
    becomes 0."""
    assert correction in CORRECTIONS
    X = np.asarray(R, dtype=np.float32)
    mean_r, sd_r = X.mean(axis=0, dtype=np.float64), X.std(axis=0, ddof=1, dtype=np.float64)
    slope = np.zeros(X.shape[1])
    X = X - mean_r.astype(np.float32)
    if correction == "depth":
        d = np.asarray(log_depth, dtype=np.float64)
        d = d - d.mean()
        slope = (d @ X.astype(np.float64)) / (d @ d)
        X = X - np.outer(d, slope).astype(np.float32)
    sd = X.std(axis=0, ddof=1, dtype=np.float64)
    ok = sd > 1e-12
    X = np.where(ok, X / np.where(ok, sd, 1.0).astype(np.float32), 0.0).astype(np.float32)
    np.clip(X, -SCALE_MAX, SCALE_MAX, out=X)
    return X, {"mean_r": mean_r, "sd_r": sd_r, "depth_slope": slope, "sd_corrected": sd}


def third_best(rank_matrix):
    """Rank of a gene over the samples: its MIN_SAMPLES-th best within-sample rank (inf if ranked in fewer)."""
    return np.sort(rank_matrix, axis=1)[:, MIN_SAMPLES - 1]


def ranks_desc(values, eligible):
    """Within-sample rank (1 = largest value) among the eligible genes; inf for the others.
    values, eligible: genes x samples."""
    out = np.full(values.shape, np.inf)
    for j in range(values.shape[1]):
        e = eligible[:, j] & np.isfinite(values[:, j])
        order = np.argsort(-values[e, j], kind="stable")
        rk = np.empty(e.sum())
        rk[order] = np.arange(1, e.sum() + 1)
        out[e, j] = rk
    return out


def entry_ranks(excess, z, ge2, var, cen, min_ge2=MIN_CELLS_GE2):
    """Entry K (third-best within-sample rank) of every gene under the three rules of step 2, from the genes x samples
    matrices of the statistics. Returns {"v_null", "v_abs", "m"} and the direction of M (+1 cytoplasmic, -1 nuclear).
      V-null  excess; ranked in a sample only if cells_ge2 >= min_ge2 and z >= Z_MIN there
      V-abs   var_r; all genes
      M       |mean_r_centred|; all genes; the qualifying samples must agree on the direction"""
    passing = (ge2 >= min_ge2) & (z >= Z_MIN)
    everything = np.ones_like(passing)
    entry = {"v_null": third_best(ranks_desc(excess, passing)), "v_abs": third_best(ranks_desc(var, everything))}
    rank_m = ranks_desc(np.abs(cen), everything)
    by_dir = {d: third_best(np.where(np.sign(cen) == d, rank_m, np.inf)) for d in (1, -1)}
    entry["m"] = np.minimum(by_dir[1], by_dir[-1])
    return entry, np.where(by_dir[1] <= by_dir[-1], 1, -1)


def top_genes(entry, tie, size):
    """Boolean mask of the first `size` genes by entry K (ties broken by `tie`, smaller first; unranked genes excluded)."""
    order = np.lexsort((tie, entry))
    order = order[np.isfinite(entry[order])][:size]
    mask = np.zeros(len(entry), bool)
    mask[order] = True
    return mask
