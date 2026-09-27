"""Compartment statistics: depths, leave-one-gene-out offsets, per-gene localisation (log odds ratio),
reliability against binomial counting noise, gene classes, rarefaction, half-splits, spatial tiles and
split-half concordance. Everything is vectorised over genes on the sparse nonzero pattern.

Definitions (per cell i, gene g): k_ig = cytoplasmic count, n_ig = k_ig + nuclear count (binomial trials),
M_i = cytoplasmic depth, Q_i = nuclear depth over the panel. The leave-one-gene-out offset is the logit of the
cell's cytoplasmic fraction over all OTHER genes, so a gene's own log odds ratio (beta_g) measures how its
localisation departs from the cell's global nucleus/cytoplasm balance.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.special import expit

import config as C

BIG = 1 << 21
EPS = 1e-9


# ------------------------------------------------------------------ depths
def cell_depths(Xn, Xc):
    q = np.asarray(Xn.sum(1)).ravel().astype(np.float64)
    m = np.asarray(Xc.sum(1)).ravel().astype(np.float64)
    tot = q + m
    return pd.DataFrame({"nuc_depth": q, "cyto_depth": m, "in_cell_depth": tot,
                         "nuc_frac": np.where(tot > 0, q / np.maximum(tot, 1), np.nan),
                         "d_match": np.minimum(q, m)})


def logit(p):
    return np.log(p) - np.log1p(-p)


# ------------------------------------------------------------------ paired sparse structure
class Pair:
    """Nuclear and cytoplasmic counts on one shared nonzero pattern (entries where n_ig > 0)."""

    def __init__(self, Xn, Xc, M=None, Q=None):
        """M, Q: optional per-cell cytoplasmic / nuclear depths to use for the offset (e.g. full-panel depths when
        Xn, Xc are pathway aggregates)."""
        P = (Xn.astype(np.int64) * BIG + Xc.astype(np.int64)).tocsr()
        P.sum_duplicates()
        P.eliminate_zeros()
        self.shape = P.shape
        self.indptr = P.indptr
        self.cols = P.indices.astype(np.int64)
        self.nuc = (P.data // BIG).astype(np.float64)
        self.k = (P.data % BIG).astype(np.float64)
        self.n = self.nuc + self.k
        self.rows = np.repeat(np.arange(P.shape[0]), np.diff(P.indptr)).astype(np.int64)
        self.M = np.asarray(Xc.sum(1)).ravel().astype(np.float64) if M is None else np.asarray(M, dtype=np.float64)
        self.Q = np.asarray(Xn.sum(1)).ravel().astype(np.float64) if Q is None else np.asarray(Q, dtype=np.float64)
        self.n_cells, self.n_genes = P.shape

    def loo_offset(self, clip=None):
        """logit of the cytoplasmic fraction over all other genes, per (cell, gene) entry."""
        clip = clip or C.OFFSET_CLIP
        m = self.M[self.rows] - self.k
        q = self.Q[self.rows] - (self.n - self.k)
        off = logit((m + 0.5) / (m + q + 1.0))
        return np.clip(off, -clip, clip)

    def subset_rows(self, mask):
        """Restrict to a subset of cells (boolean mask) -- entries and depths."""
        keep = mask[self.rows]
        new = Pair.__new__(Pair)
        new.shape = (int(mask.sum()), self.n_genes)
        remap = np.cumsum(mask) - 1
        new.rows = remap[self.rows[keep]]
        new.cols = self.cols[keep]
        new.nuc = self.nuc[keep]; new.k = self.k[keep]; new.n = self.n[keep]
        new.M = self.M[mask]; new.Q = self.Q[mask]
        new.n_cells, new.n_genes = new.shape
        new.indptr = None
        return new

    def gene_mean_total(self):
        return np.bincount(self.cols, weights=self.n, minlength=self.n_genes) / self.n_cells


def aggregate(X, M):
    """Sum counts over gene sets: (cells x genes) @ (genes x sets) -> csr (cells x sets)."""
    return sparse.csr_matrix(sparse.csr_matrix(X, dtype=np.int64) @ sparse.csr_matrix(M, dtype=np.int64)).astype(np.int32)


def pair_from_sets(Xn, Xc, M, depth_M, depth_Q):
    """A Pair whose 'genes' are gene sets (pathways or matched random sets), offsets from full-panel depths."""
    return Pair(aggregate(Xn, M), aggregate(Xc, M), M=depth_M, Q=depth_Q)


# ------------------------------------------------------------------ per-gene log odds ratio
def binom_dev_entries(k, n, p):
    p = np.clip(p, EPS, 1 - EPS)
    t1 = np.where(k > 0, k * (np.log(np.maximum(k, EPS)) - np.log(n * p)), 0.0)
    r = n - k
    t2 = np.where(r > 0, r * (np.log(np.maximum(r, EPS)) - np.log(n * (1 - p))), 0.0)
    return 2.0 * (t1 + t2)


def fit_intercepts(pair, offset, groups=None, n_groups=None, n_iter=30):
    """Offset-binomial intercepts per (group, gene) by damped Newton, vectorised over all cells.

    groups: int array per ENTRY (default: a single group). Returns beta (n_groups x G), plus per-(group, gene)
    trial sums N and success sums K.
    """
    G = pair.n_genes
    if groups is None:
        groups = np.zeros(len(pair.cols), dtype=np.int64)
        n_groups = 1
    idx = groups * G + pair.cols
    size = n_groups * G
    K = np.bincount(idx, weights=pair.k, minlength=size)
    N = np.bincount(idx, weights=pair.n, minlength=size)
    beta = np.zeros(size)
    # initialise at the pooled empirical logit corrected for the mean offset
    off_sum = np.bincount(idx, weights=offset * pair.n, minlength=size)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(N > 0, logit((K + 0.5) / (N + 1.0)) - off_sum / np.maximum(N, 1), 0.0)
    for _ in range(n_iter):
        eta = offset + beta[idx]
        p = expit(eta)
        g = np.bincount(idx, weights=pair.k - pair.n * p, minlength=size)
        h = np.bincount(idx, weights=pair.n * p * (1 - p), minlength=size)
        step = np.where(h > 1e-12, g / np.maximum(h, 1e-12), 0.0)
        step = np.clip(step, -1.5, 1.5)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-7:
            break
    beta = np.where(N > 0, beta, np.nan)
    return beta.reshape(n_groups, G), K.reshape(n_groups, G), N.reshape(n_groups, G)


def per_gene_logor(pair, offset, genes):
    """Intercept-only offset-binomial per gene with quasi-binomial SE and dispersion."""
    G = pair.n_genes
    beta, K, N = fit_intercepts(pair, offset)
    beta, K, N = beta[0], K[0], N[0]
    p = expit(offset + np.nan_to_num(beta)[pair.cols])
    dev = np.bincount(pair.cols, weights=binom_dev_entries(pair.k, pair.n, p), minlength=G)
    pearson = np.bincount(pair.cols, weights=(pair.k - pair.n * p) ** 2 / np.maximum(pair.n * p * (1 - p), EPS),
                          minlength=G)
    n_units = np.bincount(pair.cols, minlength=G).astype(float)
    h = np.bincount(pair.cols, weights=pair.n * p * (1 - p), minlength=G)
    phi = np.where(n_units > 1, pearson / np.maximum(n_units - 1, 1), np.nan)
    se = np.sqrt(np.maximum(phi, 1.0) / np.maximum(h, EPS))
    # null (offset only) deviance for a per-gene "localisation departs from the global balance" statistic
    dev0 = np.bincount(pair.cols, weights=binom_dev_entries(pair.k, pair.n, expit(offset)), minlength=G)
    df = pd.DataFrame({"gene": genes, "K_cyto": K, "N_total": N, "n_units": n_units,
                       "mean_total": N / pair.n_cells, "cyto_frac_pooled": np.where(N > 0, K / np.maximum(N, 1), np.nan),
                       "beta": beta, "se": se, "ci_lo": beta - 1.96 * se, "ci_hi": beta + 1.96 * se,
                       "phi": phi, "dev_null": dev0, "dev_fit": dev,
                       "dev_expl_offset": np.where(dev0 > 0, (dev0 - dev) / np.maximum(dev0, EPS), np.nan)})
    return df


def per_gene_reliability(pair, offset, beta, n_min=None, B=None, seed=None, min_units=None):
    """Fraction of the per-unit localisation variance that exceeds binomial counting noise.

    y = logit((k+0.5)/(n-k+0.5)) - offset on units with n >= n_min. Parametric null: k_sim ~ Binom(n, expit(offset +
    beta_g)); rel_param = 1 - mean_b Var(y_sim) / Var(y_obs). Analytic cross-check: 1/(n p (1-p)) with per-unit p.
    """
    n_min = n_min or C.N_MIN_LOC
    B = B or C.B_NULL
    seed = C.SEED if seed is None else seed
    min_units = min_units or C.MIN_UNITS
    G = pair.n_genes
    ok = pair.n >= n_min
    cols, k, n, off = pair.cols[ok], pair.k[ok], pair.n[ok], offset[ok]
    cnt = np.bincount(cols, minlength=G).astype(float)

    def var_by_gene(y):
        s1 = np.bincount(cols, weights=y, minlength=G)
        s2 = np.bincount(cols, weights=y * y, minlength=G)
        with np.errstate(invalid="ignore", divide="ignore"):
            v = (s2 - s1 ** 2 / np.maximum(cnt, 1)) / np.maximum(cnt - 1, 1)
        return np.where(cnt > 1, v, np.nan)

    y = logit((k + 0.5) / (n + 1.0)) - off
    var_obs = var_by_gene(y)
    p_hat = expit(off + np.nan_to_num(beta)[cols])
    rng = np.random.default_rng(seed)
    var_sim = np.zeros(G)
    for _ in range(B):
        ks = rng.binomial(n.astype(np.int64), p_hat).astype(np.float64)
        var_sim += var_by_gene(logit((ks + 0.5) / (n + 1.0)) - off)
    var_sim /= B
    pe = np.clip((k + 0.5) / (n + 1.0), 1e-3, 1 - 1e-3)
    samp = np.bincount(cols, weights=1.0 / (n * pe * (1 - pe)), minlength=G) / np.maximum(cnt, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rel_param = 1.0 - var_sim / var_obs
        rel_analytic = 1.0 - samp / var_obs
    bad = cnt < min_units
    rel_param[bad] = np.nan
    rel_analytic[bad] = np.nan
    return pd.DataFrame({"n_units_rel": cnt, "var_obs": var_obs, "var_noise_param": var_sim,
                         "var_noise_analytic": samp, "rel_param": rel_param, "rel_analytic": rel_analytic,
                         "rel_param_clip": np.clip(rel_param, 0, 1)})


def classify_genes(df, thr=None, min_counts=None):
    thr = thr or C.LOG2
    min_counts = min_counts or C.MIN_TOTAL_COUNTS_GENE
    cls = np.full(len(df), "ambiguous", dtype=object)
    b, lo, hi = df["beta"].values, df["ci_lo"].values, df["ci_hi"].values
    cls[(b < -thr) & (hi < 0)] = "nuclear-retained"
    cls[(b > thr) & (lo > 0)] = "cytoplasm-enriched"
    cls[(np.abs(b) <= thr)] = "balanced"
    cls[(df["N_total"].values < min_counts) | ~np.isfinite(b)] = "low-coverage"
    return pd.Series(cls, index=df.index)


# ------------------------------------------------------------------ rarefaction and half-splits
def _downsample(X, target, seed):
    """Hypergeometric thinning of each row of X to `target[i]` transcripts (rows with fewer are unchanged)."""
    a = ad.AnnData(X=sparse.csr_matrix(X, dtype=np.int64))
    sc.pp.downsample_counts(a, counts_per_cell=np.asarray(target).astype(np.int64), random_state=seed)
    out = sparse.csr_matrix(a.X).astype(np.int32)
    out.eliminate_zeros()
    return out


def rarefy_matched(Xn, Xc, d, seed=None):
    """Downsample nuclear, cytoplasmic and total counts to the same per-cell depth d = min(Q_i, M_i)."""
    seed = C.SEED if seed is None else seed
    d = np.asarray(d).astype(np.int64)
    nuc_m = _downsample(Xn, d, seed)
    cyto_m = _downsample(Xc, d, seed + 1)
    tot = (Xn + Xc).tocsr()
    total_m = _downsample(tot, d, seed + 2)
    return nuc_m, cyto_m, total_m, tot.astype(np.int32)


def thin_half(X, seed=None):
    """Two complementary halves by binomial thinning (p = 0.5): disjoint reads, equal expected depth."""
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    X = sparse.csr_matrix(X, dtype=np.int64)
    h1 = X.copy()
    h1.data = rng.binomial(X.data, 0.5)
    h2 = X - h1
    h1.eliminate_zeros(); h2.eliminate_zeros()
    return h1.astype(np.int32), h2.astype(np.int32)


# ------------------------------------------------------------------ spatial tiles
def spatial_tiles(xy, tile_um, min_cells=None):
    """Non-overlapping grid tiles; tiles with < min_cells cells get id -1."""
    min_cells = min_cells or C.TILE_MIN_CELLS
    x0, y0 = np.nanmin(xy[:, 0]), np.nanmin(xy[:, 1])
    ix = np.floor((xy[:, 0] - x0) / tile_um).astype(np.int64)
    iy = np.floor((xy[:, 1] - y0) / tile_um).astype(np.int64)
    raw = ix * (iy.max() + 1) + iy
    uniq, inv, cnt = np.unique(raw, return_inverse=True, return_counts=True)
    keep = cnt >= min_cells
    new_id = np.full(len(uniq), -1, dtype=np.int64)
    new_id[keep] = np.arange(keep.sum())
    return new_id[inv]


def choose_tile_um(xy, target_cells=None, candidates=None, min_cells=None):
    target_cells = target_cells or C.TILE_TARGET_CELLS
    candidates = candidates or C.TILE_CANDIDATES_UM
    best, best_gap = candidates[0], np.inf
    for t in candidates:
        tid = spatial_tiles(xy, t, min_cells=1)
        med = np.median(np.bincount(tid[tid >= 0]))
        gap = abs(np.log(med / target_cells))
        if gap < best_gap:
            best, best_gap = t, gap
    return best


def tile_indicator(tile_id):
    ok = tile_id >= 0
    n_t = int(tile_id.max()) + 1 if ok.any() else 0
    rows = np.where(ok)[0]
    return sparse.csr_matrix((np.ones(len(rows)), (rows, tile_id[ok])), shape=(len(tile_id), n_t))


def pool(X, T):
    """Sum rows of X within tiles: (tiles x genes)."""
    return sparse.csr_matrix(T.T @ sparse.csr_matrix(X)).astype(np.int64)


# ------------------------------------------------------------------ concordance
def _lognorm(X, depth=None, scale=None):
    """log1p(count / depth * scale) with sparsity preserved; depth defaults to the row sum."""
    X = sparse.csr_matrix(X, dtype=np.float64)
    if depth is None:
        depth = np.asarray(X.sum(1)).ravel()
    depth = np.maximum(np.asarray(depth, dtype=np.float64), 1.0)
    scale = scale or np.median(depth[depth > 1])
    Y = sparse.diags(scale / depth) @ X
    Y = sparse.csr_matrix(Y)
    Y.data = np.log1p(Y.data)
    return Y


def col_corr(A, B, w=None):
    """Column-wise (weighted) Pearson correlation of two sparse matrices with identical shape."""
    A = sparse.csr_matrix(A); B = sparse.csr_matrix(B)
    n = A.shape[0]
    if w is None:
        w = np.ones(n)
    W = w.sum()
    sa = np.asarray((sparse.diags(w) @ A).sum(0)).ravel() / W
    sb = np.asarray((sparse.diags(w) @ B).sum(0)).ravel() / W
    saa = np.asarray((sparse.diags(w) @ A.multiply(A)).sum(0)).ravel() / W
    sbb = np.asarray((sparse.diags(w) @ B.multiply(B)).sum(0)).ravel() / W
    sab = np.asarray((sparse.diags(w) @ A.multiply(B)).sum(0)).ravel() / W
    va, vb, cab = saa - sa ** 2, sbb - sb ** 2, sab - sa * sb
    with np.errstate(invalid="ignore", divide="ignore"):
        r = cab / np.sqrt(va * vb)
    r[(va <= 0) | (vb <= 0)] = np.nan
    return r


def spearman_brown(r):
    with np.errstate(invalid="ignore"):
        return 2 * r / (1 + r)


def concordance(nuc_m, cyto_m, genes, seed=None, rel_min=None, w=None, gene_mask=None):
    """Split-half concordance of nuclear vs cytoplasmic expression across units (cells or tiles).

    r_nc: cross-compartment correlation at matched depth; r_nn, r_cc: correlation of complementary halves within a
    compartment (Spearman-Brown corrected to full depth: rel_n, rel_c); r_true = r_nc / sqrt(rel_n * rel_c).
    """
    seed = C.SEED if seed is None else seed
    rel_min = rel_min or C.REL_MIN
    if gene_mask is not None:
        nuc_m = sparse.csr_matrix(nuc_m)[:, gene_mask]
        cyto_m = sparse.csr_matrix(cyto_m)[:, gene_mask]
        genes = np.asarray(genes)[gene_mask]
    n1, n2 = thin_half(nuc_m, seed)
    c1, c2 = thin_half(cyto_m, seed + 7)
    scale = np.median(np.asarray(sparse.csr_matrix(nuc_m).sum(1)).ravel())
    Vn, Vc = _lognorm(nuc_m, scale=scale), _lognorm(cyto_m, scale=scale)
    Vn1, Vn2 = _lognorm(n1, scale=scale / 2), _lognorm(n2, scale=scale / 2)
    Vc1, Vc2 = _lognorm(c1, scale=scale / 2), _lognorm(c2, scale=scale / 2)
    r_nc = col_corr(Vn, Vc, w)
    r_nn = col_corr(Vn1, Vn2, w)
    r_cc = col_corr(Vc1, Vc2, w)
    r_nc_half = 0.5 * (col_corr(Vn1, Vc1, w) + col_corr(Vn2, Vc2, w))
    rel_n, rel_c = spearman_brown(r_nn), spearman_brown(r_cc)
    with np.errstate(invalid="ignore", divide="ignore"):
        r_true = r_nc / np.sqrt(rel_n * rel_c)
    r_true[~((rel_n >= rel_min) & (rel_c >= rel_min))] = np.nan
    tot = sparse.csr_matrix(nuc_m) + sparse.csr_matrix(cyto_m)
    mean_unit = np.asarray(tot.mean(0)).ravel()
    return pd.DataFrame({"gene": genes, "mean_per_unit": mean_unit, "n_units": nuc_m.shape[0],
                         "r_nc": r_nc, "r_nc_half": r_nc_half, "r_nn_half": r_nn, "r_cc_half": r_cc,
                         "rel_n": rel_n, "rel_c": rel_c, "r_true": r_true})


def bootstrap_tiles(stat_fn, tile_ids, n_boot=None, seed=None):
    """Percentile CI of a vector statistic under tile (block) resampling. stat_fn(weights) -> array."""
    n_boot = n_boot or C.N_BOOT
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    tiles, inv = np.unique(tile_ids, return_inverse=True)
    out = []
    for _ in range(n_boot):
        mult = np.bincount(rng.integers(0, len(tiles), len(tiles)), minlength=len(tiles)).astype(float)
        out.append(stat_fn(mult[inv]))
    out = np.asarray(out, dtype=float)
    return np.nanpercentile(out, 2.5, axis=0), np.nanpercentile(out, 97.5, axis=0)
