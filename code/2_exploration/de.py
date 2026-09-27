"""Compartment-specific differential expression across a partition: pseudobulk log fold changes per cluster in
the nuclear and the cytoplasmic half-matrices (equal depth), tile block-bootstrap CIs for their difference,
Wilcoxon top-k lists for the conventional overlap panel, and the cytoplasm-only / nucleus-only calls."""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

import config as C
import compartment as CP


def _cluster_indicator(labels):
    labels = np.asarray(labels)
    L = labels.max() + 1
    return sparse.csr_matrix((np.ones(len(labels)), (np.arange(len(labels)), labels)), shape=(len(labels), L)), L


def pseudobulk_logfc(X, labels, w=None):
    """log2 rate ratio (cluster vs rest) per gene with 0.5 pseudo-counts; rate = counts / depth."""
    X = sparse.csr_matrix(X, dtype=np.float64)
    Cm, L = _cluster_indicator(labels)
    if w is None:
        w = np.ones(X.shape[0])
    Wd = sparse.diags(w)
    K_in = np.asarray((Cm.T @ (Wd @ X)).todense())          # L x G
    depth = np.asarray(X.sum(1)).ravel() * w
    D_in = np.asarray(Cm.T @ depth).ravel()                   # L
    K_all, D_all = K_in.sum(0, keepdims=True), D_in.sum()
    K_out, D_out = K_all - K_in, D_all - D_in
    lfc = np.log2((K_in + 0.5) / np.maximum(D_in, 1)[:, None]) - np.log2((K_out + 0.5) / np.maximum(D_out, 1)[:, None])
    return lfc, K_in, D_in


def paired_logfc(Xn_h, Xc_h, labels, genes, tile_ids, gene_mask=None, n_boot=None, seed=None):
    """lfc per cluster in nuclear and cytoplasmic halves, delta = lfc_cyto - lfc_nuc with tile bootstrap CI."""
    n_boot = n_boot or C.N_BOOT
    Xn_h = sparse.csr_matrix(Xn_h); Xc_h = sparse.csr_matrix(Xc_h)
    if gene_mask is not None:
        Xn_h, Xc_h, genes = Xn_h[:, gene_mask], Xc_h[:, gene_mask], np.asarray(genes)[gene_mask]
    lfc_n, K_n, D_n = pseudobulk_logfc(Xn_h, labels)
    lfc_c, K_c, D_c = pseudobulk_logfc(Xc_h, labels)
    delta = lfc_c - lfc_n

    def stat(w):
        a, _, _ = pseudobulk_logfc(Xc_h, labels, w)
        b, _, _ = pseudobulk_logfc(Xn_h, labels, w)
        return (a - b).ravel()

    ok = tile_ids >= 0
    lo, hi = CP.bootstrap_tiles(lambda w: stat(w * ok), tile_ids, n_boot=n_boot, seed=seed)
    L, G = delta.shape
    df = pd.DataFrame({"cluster": np.repeat(np.arange(L), G), "gene": np.tile(genes, L),
                       "lfc_nuc": lfc_n.ravel(), "lfc_cyto": lfc_c.ravel(), "delta": delta.ravel(),
                       "delta_lo": lo, "delta_hi": hi,
                       "K_nuc_in": K_n.ravel(), "K_cyto_in": K_c.ravel()})
    return df


def call_de(paired, lfc_de=None, lfc_null=None):
    lfc_de = lfc_de or C.LOGFC_DE
    lfc_null = lfc_null or C.LOGFC_NULL
    ci_excl = (paired["delta_lo"] > 0) | (paired["delta_hi"] < 0)
    out = paired.copy()
    out["cyto_only"] = (out["lfc_cyto"].abs() > lfc_de) & (out["lfc_nuc"].abs() < lfc_null) & ci_excl
    out["nuc_only"] = (out["lfc_nuc"].abs() > lfc_de) & (out["lfc_cyto"].abs() < lfc_null) & ci_excl
    out["shared"] = (out["lfc_nuc"].abs() > lfc_de) & (out["lfc_cyto"].abs() > lfc_de) & \
                    (np.sign(out["lfc_nuc"]) == np.sign(out["lfc_cyto"]))
    return out


def wilcoxon_topk(X, genes, labels, gene_mask=None, k=None):
    """Top-k Wilcoxon markers per cluster on log-normalised data (conventional overlap panel)."""
    k = k or C.TOPK_DE
    X = sparse.csr_matrix(X)
    if gene_mask is not None:
        X, genes = X[:, gene_mask], np.asarray(genes)[gene_mask]
    a = ad.AnnData(X=X.astype(np.float32))
    a.var_names = [str(g) for g in genes]
    a.obs["cl"] = pd.Categorical([str(l) for l in labels])
    depth = np.asarray(a.X.sum(1)).ravel()
    sc.pp.normalize_total(a, target_sum=float(np.median(depth[depth > 0])))
    sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, "cl", method="wilcoxon", n_genes=k)
    names = a.uns["rank_genes_groups"]["names"]
    return {c: [str(x) for x in names[c]] for c in a.obs["cl"].cat.categories}


def topk_overlap(lists_a, lists_b):
    rows = []
    for c in lists_a:
        A, B = set(lists_a[c]), set(lists_b.get(c, []))
        rows.append({"cluster": c, "n_shared": len(A & B), "n_a_only": len(A - B), "n_b_only": len(B - A),
                     "jaccard": len(A & B) / max(len(A | B), 1)})
    return pd.DataFrame(rows)


def cluster_centroids(X, labels, gene_mask=None):
    """Depth-normalised log centroid per cluster (for cross-sample cluster recurrence)."""
    X = sparse.csr_matrix(X, dtype=np.float64)
    if gene_mask is not None:
        X = X[:, gene_mask]
    Cm, L = _cluster_indicator(labels)
    K = np.asarray((Cm.T @ X).todense())
    D = K.sum(1, keepdims=True)
    return np.log1p(K / np.maximum(D, 1) * 1e4)


def stratified_paired_logfc(Xn_h, Xc_h, labels, strata, genes, tile_ids, gene_mask=None, n_boot=None, seed=None):
    """Within-stratum cluster-vs-rest log2 FC in each compartment, strata-weighted by cluster depth, with a tile
    block-bootstrap CI for delta = lfc_cyto - lfc_nuc. Cells with stratum -1 are excluded."""
    n_boot = n_boot or C.N_BOOT
    labels = np.asarray(labels); strata = np.asarray(strata)
    ok = (strata >= 0) & (labels >= 0)
    Xn_h = sparse.csr_matrix(Xn_h, dtype=np.float64)[ok]; Xc_h = sparse.csr_matrix(Xc_h, dtype=np.float64)[ok]
    labels, strata, tile_ids = labels[ok], strata[ok], np.asarray(tile_ids)[ok]
    if gene_mask is not None:
        Xn_h, Xc_h, genes = Xn_h[:, gene_mask], Xc_h[:, gene_mask], np.asarray(genes)[gene_mask]
    L, S = labels.max() + 1, strata.max() + 1
    joint = strata * L + labels
    J = sparse.csr_matrix((np.ones(len(joint)), (np.arange(len(joint)), joint)), shape=(len(joint), S * L))
    dn = np.asarray(Xn_h.sum(axis=1)).ravel(); dc = np.asarray(Xc_h.sum(axis=1)).ravel()

    def lfc(X, depth, w):
        K = np.asarray((J.T @ (sparse.diags(w) @ X)).todense()).reshape(S, L, -1)   # S x L x G
        D = np.asarray(J.T @ (w * depth)).ravel().reshape(S, L)                      # S x L
        Ks, Ds = K.sum(1, keepdims=True), D.sum(1, keepdims=True)
        r_in = np.log2((K + 0.5) / np.maximum(D, 1)[:, :, None])
        r_out = np.log2((Ks - K + 0.5) / np.maximum(Ds - D, 1)[:, :, None])
        wgt = D / np.maximum(D.sum(0, keepdims=True), 1e-12)                          # strata weights per cluster
        return (wgt[:, :, None] * (r_in - r_out)).sum(0)                              # L x G

    ones = np.ones(len(labels))
    lfc_n, lfc_c = lfc(Xn_h, dn, ones), lfc(Xc_h, dc, ones)
    delta = lfc_c - lfc_n
    lo, hi = CP.bootstrap_tiles(lambda w: (lfc(Xc_h, dc, w) - lfc(Xn_h, dn, w)).ravel(), tile_ids, n_boot=n_boot, seed=seed)
    G = delta.shape[1]
    return pd.DataFrame({"cluster": np.repeat(np.arange(L), G), "gene": np.tile(genes, L),
                         "lfc_nuc": lfc_n.ravel(), "lfc_cyto": lfc_c.ravel(), "delta": delta.ravel(),
                         "delta_lo": lo, "delta_hi": hi})
