"""Association of per-gene compartment quantities with categorical input axes, vectorised over genes.

Three models per gene and axis (levels l):
  frac : k_ig ~ Binomial(n_ig, expit(offset_ig + beta_lg))   cytoplasmic fraction, leave-one-gene-out offset
  cyto : k_ig ~ Poisson(M_i * rate_lg)                       cytoplasmic count, cytoplasmic depth as exposure
  nuc  : (n_ig - k_ig) ~ Poisson(Q_i * rate_lg)              nuclear count, nuclear depth as exposure
Deviance explained = 1 - D_axis / D_null (null = one level). Poisson fits are closed form; the binomial fit is a
per-(level, gene) Newton iteration on group sums. Label permutations give the overfitting baseline.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit

import config as C
import compartment as CP

EPS = 1e-12


# ------------------------------------------------------------------ Poisson with log exposure (closed form)
def poisson_axis(X, exposure, labels):
    """X: cells x genes counts (csr); exposure: cells; labels: ints 0..L-1 (-1 = excluded).

    Returns rate (L x G), dev_null (G), dev_axis (G), E (L) exposure sums, K (L x G).
    """
    X = sparse.csr_matrix(X, dtype=np.float64)
    labels = np.asarray(labels)
    ok = (labels >= 0) & (exposure > 0)
    X, exposure, labels = X[ok], exposure[ok], labels[ok]
    L = labels.max() + 1
    Cm = sparse.csr_matrix((np.ones(len(labels)), (np.arange(len(labels)), labels)), shape=(len(labels), L))
    K = np.asarray((Cm.T @ X).todense())                       # L x G
    E = np.asarray(Cm.T @ exposure).ravel()                     # L
    Kg = K.sum(0)
    Eg = E.sum()
    # sum_nz k log k and sum_nz k log exposure per gene
    coo = X.tocoo()
    klogk = np.bincount(coo.col, weights=coo.data * np.log(coo.data), minlength=X.shape[1])
    kloge = np.bincount(coo.col, weights=coo.data * np.log(exposure[coo.row]), minlength=X.shape[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        term_axis = np.where(K > 0, K * np.log(K / E[:, None]), 0.0).sum(0)
        term_null = np.where(Kg > 0, Kg * np.log(Kg / Eg), 0.0)
    dev_axis = 2 * (klogk - kloge - term_axis)
    dev_null = 2 * (klogk - kloge - term_null)
    with np.errstate(divide="ignore"):
        rate = np.log((K + 0.5) / E[:, None])
    return rate, dev_null, dev_axis, E, K


# ------------------------------------------------------------------ offset-binomial with a categorical axis
def binom_axis(pair, offset, labels_cell, n_iter=25):
    """labels_cell: ints per cell (-1 = excluded). Returns beta (L x G), dev_null, dev_axis, N (L x G)."""
    labels_cell = np.asarray(labels_cell)
    lab_e = labels_cell[pair.rows]
    ok = lab_e >= 0
    sub = _subset_entries(pair, ok)
    off = offset[ok]
    lab_e = lab_e[ok]
    L = labels_cell.max() + 1
    beta, K, N = CP.fit_intercepts(sub, off, groups=lab_e, n_groups=L)
    beta0, _, _ = CP.fit_intercepts(sub, off)
    G = pair.n_genes
    idx = lab_e * G + sub.cols
    p_axis = expit(off + np.nan_to_num(beta.ravel())[idx])
    p_null = expit(off + np.nan_to_num(beta0.ravel())[sub.cols])
    dev_axis = np.bincount(sub.cols, weights=CP.binom_dev_entries(sub.k, sub.n, p_axis), minlength=G)
    dev_null = np.bincount(sub.cols, weights=CP.binom_dev_entries(sub.k, sub.n, p_null), minlength=G)
    return beta, dev_null, dev_axis, N


class _Entries:
    pass


def _subset_entries(pair, mask):
    s = _Entries()
    s.cols = pair.cols[mask]; s.k = pair.k[mask]; s.n = pair.n[mask]
    s.n_genes = pair.n_genes
    return s


# ------------------------------------------------------------------ effect sizes
def weighted_level_sd(effect, weights):
    """Count-weighted SD of level effects (L x G with L weights, or L x G weights)."""
    effect = np.asarray(effect, dtype=float)
    w = np.asarray(weights, dtype=float)
    if w.ndim == 1:
        w = np.repeat(w[:, None], effect.shape[1], axis=1)
    w = np.where(np.isfinite(effect), w, 0.0)
    e = np.where(np.isfinite(effect), effect, 0.0)
    W = w.sum(0)
    mu = (w * e).sum(0) / np.maximum(W, EPS)
    var = (w * (e - mu[None, :]) ** 2).sum(0) / np.maximum(W, EPS)
    return np.where(W > 0, np.sqrt(var), np.nan)


def dev_explained(dev_null, dev_axis):
    with np.errstate(divide="ignore", invalid="ignore"):
        out = 1.0 - dev_axis / dev_null
    out[~(dev_null > 0)] = np.nan
    return out


# ------------------------------------------------------------------ driver
def run_axis(pair, offset, Xn, Xc, labels, gene_idx, genes, axis_name, n_perm=None, seed=None):
    """Three models for one axis on the genes in gene_idx; long table with permutation baselines."""
    n_perm = C.N_PERM if n_perm is None else n_perm
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    Xn_s = sparse.csr_matrix(Xn)[:, gene_idx]
    Xc_s = sparse.csr_matrix(Xc)[:, gene_idx]
    sub_pair = _pair_gene_subset(pair, gene_idx)
    sub_offset = offset_for_subset(offset, sub_pair)

    def fit(lab):
        rate_c, dn_c, da_c, E_c, _ = poisson_axis(Xc_s, pair.M, lab)
        rate_n, dn_n, da_n, E_n, _ = poisson_axis(Xn_s, pair.Q, lab)
        beta_f, dn_f, da_f, N_f = binom_axis(sub_pair, sub_offset, lab)
        return {"cyto": (dev_explained(dn_c, da_c), weighted_level_sd(rate_c, E_c)),
                "nuc": (dev_explained(dn_n, da_n), weighted_level_sd(rate_n, E_n)),
                "frac": (dev_explained(dn_f, da_f), weighted_level_sd(beta_f, N_f))}

    obs = fit(labels)
    perm = [fit(_permute_valid(labels, rng)) for _ in range(n_perm)]
    rows = []
    for model in ("cyto", "nuc", "frac"):
        de, sd = obs[model]
        de_p = np.nanmean([p[model][0] for p in perm], axis=0) if perm else np.full(len(de), np.nan)
        sd_p = np.nanmean([p[model][1] for p in perm], axis=0) if perm else np.full(len(de), np.nan)
        rows.append(pd.DataFrame({"gene": np.asarray(genes)[gene_idx], "axis": axis_name, "model": model,
                                  "n_levels": int(labels.max() + 1), "dev_expl": de, "dev_expl_perm": de_p,
                                  "dev_expl_excess": de - de_p, "sd_levels": sd, "sd_levels_perm": sd_p}))
    return pd.concat(rows, ignore_index=True)


def _permute_valid(labels, rng):
    out = labels.copy()
    ok = labels >= 0
    out[ok] = rng.permutation(labels[ok])
    return out


def _pair_gene_subset(pair, gene_idx):
    """Restrict a Pair's entries to a gene subset, renumbering columns 0..len(gene_idx)-1."""
    gene_idx = np.asarray(gene_idx)
    remap = np.full(pair.n_genes, -1, dtype=np.int64)
    remap[gene_idx] = np.arange(len(gene_idx))
    keep = remap[pair.cols] >= 0
    s = CP.Pair.__new__(CP.Pair)
    s.rows = pair.rows[keep]; s.cols = remap[pair.cols[keep]]
    s.k = pair.k[keep]; s.n = pair.n[keep]; s.nuc = pair.nuc[keep]
    s.M, s.Q = pair.M, pair.Q
    s.n_cells, s.n_genes = pair.n_cells, len(gene_idx)
    s.shape = (s.n_cells, s.n_genes)
    s.indptr = None
    s._offset_mask = keep
    return s


def offset_for_subset(offset, sub_pair):
    return offset[sub_pair._offset_mask]


def cross_labels(a, b):
    """Cross-classification of two categorical axes (-1 propagates)."""
    a, b = np.asarray(a), np.asarray(b)
    Lb = b.max() + 1
    out = a * Lb + b
    out[(a < 0) | (b < 0)] = -1
    # compress to 0..L-1
    u, inv = np.unique(out[out >= 0], return_inverse=True)
    res = np.full(len(out), -1, dtype=np.int64)
    res[out >= 0] = inv
    return res


def wide_summary(long):
    """axis x model dev_expl_excess per gene, plus the cytoplasm-specific contrasts."""
    w = long.pivot_table(index="gene", columns=["axis", "model"], values="dev_expl_excess")
    w.columns = [f"{a}__{m}" for a, m in w.columns]
    for a in long["axis"].unique():
        if f"{a}__cyto" in w and f"{a}__nuc" in w:
            w[f"{a}__cyto_minus_nuc"] = w[f"{a}__cyto"] - w[f"{a}__nuc"]
    return w.reset_index()
