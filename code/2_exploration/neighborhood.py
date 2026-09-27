"""Neighbourhood features for tumor cells from all cells of the slide: Gaussian-kernel cell-type composition,
crowding, boundary vs core, distance to the nearest non-tumor cell, niches and quantile bins."""

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors

import config as C

CHUNK = 20_000


def kernel_geometry(xy_query, xy_all, self_idx, sigma=None, cutoff=None, radii=None):
    """Row-normalised Gaussian kernel weights (query x all, self excluded) and per-query counts within radii."""
    sigma = sigma or C.SIGMA
    cutoff = cutoff or C.KERNEL_CUT
    radii = radii or C.RADII_COUNT
    nn = NearestNeighbors(radius=cutoff, algorithm="kd_tree").fit(xy_all)
    n_q = len(xy_query)
    ri, rj, rw = [], [], []
    mass = np.zeros(n_q)
    counts = {r: np.zeros(n_q, dtype=np.int32) for r in radii}
    for s in range(0, n_q, CHUNK):
        e = min(s + CHUNK, n_q)
        dist, ind = nn.radius_neighbors(xy_query[s:e], return_distance=True)
        for a, (dd, ii) in enumerate(zip(dist, ind)):
            g = s + a
            keep = ii != self_idx[g]
            dd, ii = dd[keep], ii[keep]
            for r in radii:
                counts[r][g] = int(np.count_nonzero(dd <= r))
            if len(ii) == 0:
                continue
            w = np.exp(-(dd ** 2) / (2.0 * sigma ** 2))
            mass[g] = w.sum()
            ri.append(np.full(len(ii), g, dtype=np.int64)); rj.append(ii.astype(np.int64)); rw.append(w)
    if ri:
        ri, rj, rw = np.concatenate(ri), np.concatenate(rj), np.concatenate(rw)
    else:
        ri = rj = np.array([], dtype=np.int64); rw = np.array([])
    W = sparse.coo_matrix((rw, (ri, rj)), shape=(n_q, len(xy_all))).tocsr()
    rs = np.asarray(W.sum(axis=1)).ravel()
    inv = np.where(rs > 0, 1.0 / np.maximum(rs, 1e-12), 0.0)
    W = sparse.diags(inv) @ W
    meta = pd.DataFrame({"nbhd_mass": mass})
    for r in radii:
        meta[f"n_{int(r)}"] = counts[r]
    return sparse.csr_matrix(W), meta


def composition(W, types_all, vocab=None):
    vocab = vocab or C.CELL_TYPES_18
    types_all = np.asarray(types_all)
    cols = []
    for t in vocab:
        ind = (types_all == t).astype(float)
        cols.append(np.asarray(W @ ind).ravel())
    return pd.DataFrame(np.vstack(cols).T, columns=[f"comp_{t}" for t in vocab])


def nbhd_features(ds_obs_all, tumor_cell_ids, sigma=None, cutoff=None, radii=None):
    """All neighbourhood features for the tumor cells (rows ordered like tumor_cell_ids)."""
    radii = radii or C.RADII_COUNT
    xy_all = ds_obs_all[["global_x", "global_y"]].values.astype(np.float64)
    types = ds_obs_all["cell_type_merged"].astype(str).values
    pos = pd.Series(np.arange(len(ds_obs_all)), index=ds_obs_all["cell_id"].astype(str).values)
    self_idx = pos.loc[np.asarray(tumor_cell_ids).astype(str)].values
    xy_q = xy_all[self_idx]
    W, meta = kernel_geometry(xy_q, xy_all, self_idx, sigma, cutoff, radii)
    comp = composition(W, types)
    is_tum = types == C.TUMOR_TYPE
    is_imm = np.isin(types, C.IMMUNE_TYPES)
    r0 = int(radii[0])
    # tumor / immune counts within the inner radius
    nn = NearestNeighbors(radius=radii[0]).fit(xy_all)
    n_tum = np.zeros(len(xy_q), dtype=np.int32); n_imm = np.zeros(len(xy_q), dtype=np.int32)
    for s in range(0, len(xy_q), CHUNK):
        e = min(s + CHUNK, len(xy_q))
        ind = nn.radius_neighbors(xy_q[s:e], return_distance=False)
        for a, ii in enumerate(ind):
            ii = ii[ii != self_idx[s + a]]
            n_tum[s + a] = is_tum[ii].sum(); n_imm[s + a] = is_imm[ii].sum()
    # distance to the nearest non-tumor cell
    non = np.where(~is_tum)[0]
    if len(non):
        d_non, _ = NearestNeighbors(n_neighbors=1).fit(xy_all[non]).kneighbors(xy_q)
        d_non = d_non.ravel()
    else:
        d_non = np.full(len(xy_q), np.nan)
    feat = pd.concat([meta, comp], axis=1)
    feat[f"n_tumor_{r0}"] = n_tum
    feat[f"n_immune_{r0}"] = n_imm
    feat["log1p_crowding"] = np.log1p(n_tum)
    feat["tumor_frac"] = np.where(meta[f"n_{r0}"] > 0, n_tum / np.maximum(meta[f"n_{r0}"], 1), np.nan)
    feat["immune_frac"] = np.where(meta[f"n_{r0}"] > 0, n_imm / np.maximum(meta[f"n_{r0}"], 1), np.nan)
    feat["immune_kernel"] = comp[[f"comp_{t}" for t in C.IMMUNE_TYPES]].sum(axis=1).values
    feat["stromal_kernel"] = comp[[f"comp_{t}" for t in C.STROMAL_TYPES]].sum(axis=1).values
    feat["dist_nontumor"] = d_non
    feat["log1p_dist_nontumor"] = np.log1p(d_non)
    return feat


def niche_labels(comp, k=None, seed=None):
    k = k or C.N_NICHES
    seed = C.SEED if seed is None else seed
    Z = np.sqrt(np.clip(comp.values.astype(float), 0, None))
    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(Z)
    lab = km.labels_
    order = np.argsort(-np.bincount(lab, minlength=k))
    remap = np.empty(k, dtype=int); remap[order] = np.arange(k)
    return remap[lab]


def quantile_bins(x, n=None):
    """Rank-based quantile bins 0..n-1 (NaN -> -1)."""
    n = n or C.N_BINS
    x = np.asarray(x, dtype=float)
    out = np.full(len(x), -1, dtype=np.int64)
    ok = np.isfinite(x)
    r = pd.Series(x[ok]).rank(method="first").values
    out[ok] = np.minimum((r - 1) * n // ok.sum(), n - 1).astype(np.int64)
    return out
