"""Artifact controls: abundance-bin percentiles, expression-matched control genes and random gene sets,
and stratified checks (segmentation method, area-ratio bins)."""

import numpy as np
import pandas as pd

import config as C


def expression_bins(mean_expr, n_bins=None):
    """Quantile bins over genes with mean > 0 (-1 for zero-expression genes)."""
    n_bins = n_bins or C.N_EXPR_BINS
    mean_expr = np.asarray(mean_expr, dtype=float)
    out = np.full(len(mean_expr), -1, dtype=np.int64)
    ok = mean_expr > 0
    r = pd.Series(mean_expr[ok]).rank(method="first").values
    out[ok] = np.minimum((r - 1) * n_bins // ok.sum(), n_bins - 1).astype(np.int64)
    return out


def bin_percentile(values, bins):
    """Percentile rank (0-100) of each gene's statistic among genes in the same abundance bin."""
    values = np.asarray(values, dtype=float)
    bins = np.asarray(bins)
    out = np.full(len(values), np.nan)
    for b in np.unique(bins[bins >= 0]):
        m = (bins == b) & np.isfinite(values)
        if m.sum() < 2:
            continue
        r = pd.Series(values[m]).rank(pct=True).values * 100
        out[m] = r
    return out


def add_bin_percentiles(df, cols, mean_col="mean_total", n_bins=None):
    bins = expression_bins(df[mean_col].values, n_bins)
    df = df.copy()
    df["expr_bin"] = bins
    for c in cols:
        if c in df:
            df[c + "_pct"] = bin_percentile(df[c].values, bins)
            df[c + "_abs_pct"] = bin_percentile(np.abs(df[c].values), bins)
    return df


def matched_control_genes(target, panel, mean_expr, seed=None):
    """Non-target genes matched one-to-one on log mean expression (greedy nearest neighbour, no replacement)."""
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    panel = list(panel)
    pos = {g: i for i, g in enumerate(panel)}
    logmu = np.log1p(np.asarray(mean_expr, dtype=float))
    tset = set(target)
    cand = np.array([i for i, g in enumerate(panel) if g not in tset])
    cand_val = logmu[cand]
    taken = np.zeros(len(cand), dtype=bool)
    out = []
    for j in rng.permutation(len(target)):
        d = np.abs(cand_val - logmu[pos[target[j]]])
        d[taken] = np.inf
        pick = int(np.argmin(d))
        taken[pick] = True
        out.append(panel[cand[pick]])
    return sorted(out)


def matched_random_sets(member_idx, bins, n_sets=None, seed=None, exclude_idx=None):
    """Random gene sets matched on the abundance-bin profile of a pathway (indices into the panel)."""
    n_sets = n_sets or C.N_MATCHED_SETS
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    bins = np.asarray(bins)
    member_idx = np.asarray(member_idx)
    excl = set(member_idx.tolist()) | set([] if exclude_idx is None else list(exclude_idx))
    pools = {b: np.array([i for i in np.where(bins == b)[0] if i not in excl]) for b in np.unique(bins[member_idx])}
    need = pd.Series(bins[member_idx]).value_counts()
    sets = []
    for _ in range(n_sets):
        s = []
        for b, cnt in need.items():
            pool = pools[b]
            if len(pool) == 0:
                continue
            s.extend(rng.choice(pool, size=min(cnt, len(pool)), replace=False).tolist())
        sets.append(np.array(sorted(s)))
    return sets


def set_null_percentile(stat_obs, stat_null):
    """Percentile of the observed set statistic among matched random sets (higher = more extreme)."""
    stat_null = np.asarray(stat_null, dtype=float)
    stat_null = stat_null[np.isfinite(stat_null)]
    if len(stat_null) == 0 or not np.isfinite(stat_obs):
        return np.nan
    return float((stat_null < stat_obs).mean() * 100)
