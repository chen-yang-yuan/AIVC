"""Identical clustering pipeline for nuclear / cytoplasmic / total matrices, partition comparison, cytoplasm-only
clusters, spatial coherence and cluster-level feature enrichment."""

import random

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import igraph as ig
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

import config as C


# ------------------------------------------------------------------ pipeline
def shared_gene_mask(nuc_m, cyto_m, frac_min=None):
    """Genes detected in >= frac_min of cells in either compartment (one gene set for every matrix)."""
    frac_min = frac_min or C.DETECT_FRAC_MIN
    n = nuc_m.shape[0]
    dn = np.asarray((sparse.csr_matrix(nuc_m) > 0).sum(0)).ravel() / n
    dc = np.asarray((sparse.csr_matrix(cyto_m) > 0).sum(0)).ravel() / n
    return (dn >= frac_min) | (dc >= frac_min)


def leiden_igraph(conn, resolution, seed=None, n_iterations=-1):
    """Leiden (modularity) on a symmetric sparse connectivity matrix via python-igraph."""
    seed = C.SEED if seed is None else seed
    random.seed(seed)
    tri = sparse.triu(sparse.csr_matrix(conn), k=1).tocoo()
    g = ig.Graph(n=conn.shape[0], edges=list(zip(tri.row.tolist(), tri.col.tolist())), directed=False)
    g.es["weight"] = tri.data.tolist()
    part = g.community_leiden(objective_function="modularity", weights="weight",
                              resolution=resolution, n_iterations=n_iterations)
    labels = np.asarray(part.membership, dtype=np.int64)
    # relabel by decreasing size for stable display
    order = np.argsort(-np.bincount(labels))
    remap = np.empty_like(order); remap[order] = np.arange(len(order))
    return remap[labels]


def cluster_matrix(X, gene_mask, resolutions=None, seed=None, n_pcs=None, knn=None, verbose=False):
    """normalize_total (median depth) -> log1p -> PCA -> kNN -> Leiden at each resolution."""
    resolutions = resolutions or C.RESOLUTIONS
    seed = C.SEED if seed is None else seed
    n_pcs = n_pcs or C.PCA_COMPS
    knn = knn or C.KNN
    a = ad.AnnData(X=sparse.csr_matrix(X)[:, gene_mask].astype(np.float32))
    depth = np.asarray(a.X.sum(axis=1)).ravel()
    sc.pp.normalize_total(a, target_sum=float(np.median(depth[depth > 0])))
    sc.pp.log1p(a)
    n_pcs = min(n_pcs, a.n_vars - 1, a.n_obs - 1)
    sc.pp.pca(a, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    pcs = a.obsm["X_pca"].astype(np.float32)
    labels = cluster_from_pcs(pcs, resolutions, knn, seed, verbose)
    return labels, pcs


def cluster_from_pcs(pcs, resolutions=None, knn=None, seed=None, verbose=False):
    """kNN graph on a given embedding -> Leiden at each resolution (the same path cluster_matrix uses)."""
    resolutions = resolutions or C.RESOLUTIONS
    seed = C.SEED if seed is None else seed
    knn = knn or C.KNN
    a = ad.AnnData(obs=pd.DataFrame(index=[str(i) for i in range(pcs.shape[0])]))
    a.obsm["X_pca"] = np.asarray(pcs, dtype=np.float32)
    sc.pp.neighbors(a, n_neighbors=knn, n_pcs=pcs.shape[1], use_rep="X_pca", random_state=seed)
    labels = {}
    for r in resolutions:
        labels[r] = leiden_igraph(a.obsp["connectivities"], r, seed)
        if verbose:
            print(f"      res {r}: {labels[r].max() + 1} clusters")
    return labels


# ------------------------------------------------------------------ partition comparison
def compare_partitions(a, b):
    return {"ari": float(adjusted_rand_score(a, b)), "nmi": float(normalized_mutual_info_score(a, b)),
            "n_a": int(len(np.unique(a))), "n_b": int(len(np.unique(b)))}


def jaccard_matrix(a, b):
    """Jaccard overlap of every cluster of a with every cluster of b (rows: a, cols: b)."""
    ct = pd.crosstab(pd.Series(a, name="a"), pd.Series(b, name="b"))
    inter = ct.values.astype(float)
    sa = inter.sum(1, keepdims=True)
    sb = inter.sum(0, keepdims=True)
    return pd.DataFrame(inter / (sa + sb - inter), index=ct.index, columns=ct.columns)


def cluster_recovery(ref, other):
    """For each cluster of `ref`, the best-matching cluster of `other` by Jaccard and F1."""
    ct = pd.crosstab(pd.Series(ref, name="ref"), pd.Series(other, name="other")).values.astype(float)
    sa, sb = ct.sum(1, keepdims=True), ct.sum(0, keepdims=True)
    jac = ct / (sa + sb - ct)
    f1 = 2 * ct / (sa + sb)
    rows = []
    for i in range(ct.shape[0]):
        j = int(np.argmax(f1[i]))
        rows.append({"ref_cluster": i, "size": int(sa[i, 0]), "best_other": j,
                     "best_f1": float(f1[i, j]), "best_jaccard": float(jac[i].max())})
    return pd.DataFrame(rows)


def all_pairwise(labels, names=None):
    names = names or list(labels.keys())
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = compare_partitions(labels[a], labels[b])
            d.update({"a": a, "b": b})
            rows.append(d)
    return pd.DataFrame(rows)[["a", "b", "ari", "nmi", "n_a", "n_b"]]


def cyto_only_clusters(lab, thr=None, thr_repro=None):
    """Cytoplasmic clusters absent from the nuclear and total partitions but reproduced in both cyto halves."""
    thr = thr or C.JACCARD_CYTO_ONLY
    thr_repro = thr_repro or C.JACCARD_REPRO
    cy = lab["cyto_m"]
    rows = []
    jn = jaccard_matrix(cy, lab["nuc_m"]).max(axis=1)
    jt = jaccard_matrix(cy, lab["total_m"]).max(axis=1)
    jf = jaccard_matrix(cy, lab["total_full"]).max(axis=1)
    j1 = jaccard_matrix(cy, lab["cyto_h1"]).max(axis=1)
    j2 = jaccard_matrix(cy, lab["cyto_h2"]).max(axis=1)
    sizes = pd.Series(cy).value_counts().sort_index()
    for c in sizes.index:
        rows.append({"cyto_cluster": int(c), "size": int(sizes[c]),
                     "max_jaccard_nuc": float(jn[c]), "max_jaccard_total_m": float(jt[c]),
                     "max_jaccard_total_full": float(jf[c]),
                     "repro_h1": float(j1[c]), "repro_h2": float(j2[c])})
    df = pd.DataFrame(rows)
    df["cyto_only"] = ((df["max_jaccard_nuc"] < thr) & (df["max_jaccard_total_m"] < thr)
                       & (df["max_jaccard_total_full"] < thr)
                       & (df["repro_h1"] > thr_repro) & (df["repro_h2"] > thr_repro))
    df["novel_vs_nuc"] = (df["max_jaccard_nuc"] < thr) & (df["repro_h1"] > thr_repro) & (df["repro_h2"] > thr_repro)
    return df


# ------------------------------------------------------------------ spatial coherence
def knn_indices(xy, k=None):
    k = k or C.PURITY_K
    nn = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    _, ind = nn.kneighbors(xy)
    return ind[:, 1:]


def spatial_purity(labels, nn_ind, n_perm=5, seed=None):
    """Per-cluster fraction of spatial neighbours sharing the label, and the same under label permutation."""
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)

    def purity(lab):
        same = (lab[nn_ind] == lab[:, None]).mean(axis=1)
        return pd.Series(same).groupby(lab).mean()

    obs = purity(labels)
    perm = pd.concat([purity(rng.permutation(labels)) for _ in range(n_perm)], axis=1)
    sizes = pd.Series(labels).value_counts().sort_index()
    return pd.DataFrame({"cluster": obs.index, "size": sizes.loc[obs.index].values, "purity": obs.values,
                         "purity_perm": perm.mean(axis=1).loc[obs.index].values,
                         "purity_excess": (obs - perm.mean(axis=1)).loc[obs.index].values})


def cluster_feature_enrichment(labels, features, n_perm=20, seed=None):
    """Mean of each feature per cluster and its z-score against label permutation."""
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    F = features.values.astype(float)
    L = labels.max() + 1
    cnt = np.bincount(labels, minlength=L).astype(float)

    def means(lab):
        return np.vstack([np.bincount(lab, weights=F[:, j], minlength=L) for j in range(F.shape[1])]).T / np.maximum(cnt, 1)[:, None]

    obs = means(labels)
    perms = np.stack([means(rng.permutation(labels)) for _ in range(n_perm)])
    mu, sd = perms.mean(0), perms.std(0) + 1e-12
    z = (obs - mu) / sd
    long = []
    for c in range(L):
        for j, f in enumerate(features.columns):
            long.append({"cluster": c, "size": int(cnt[c]), "feature": f, "mean": obs[c, j],
                         "mean_perm": mu[c, j], "z": z[c, j]})
    return pd.DataFrame(long)


def cluster_depth(labels, depth):
    s = pd.DataFrame({"cluster": labels, "depth": depth}).groupby("cluster")["depth"]
    return pd.DataFrame({"cluster": s.median().index, "median_depth": s.median().values, "size": s.size().values})
