"""Stress axes of tumor cells and the analysis helpers of 6_stress_substates.ipynb.

Design: plans/stress_substate_exploration.md. Every spatial axis is computed within one sample from all cells of
the slide (coordinates in um, obs["cell_type_merged"]); expression axes use whole-cell log-normalised counts of the
tumor cells; thresholds (tertiles, z-scores) are set on the pooled cells by the notebook.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.linear_model import Ridge
from sklearn.metrics import adjusted_mutual_info_score, roc_auc_score
from sklearn.neighbors import KDTree, NearestNeighbors

# ----------------------------------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------------------------------
SAMPLES = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate"]
SHORT = {s: s.replace("Xenium_5K_", "") for s in SAMPLES}
TUMOR = "Malignant cell"
EXCLUDED_TYPES = ["Unknown", "Mixed"]           # not a reference cell anywhere, not in the tumor-fraction denominator

# Coarse reference sets (plans/stress_substate_exploration.md section 3.1)
REF_SETS = {"endothelial": ["Endothelial cell"],
            "fibroblast": ["Fibroblast (CAF)"],
            "myeloid": ["Myeloid cell", "Dendritic cell"],
            "tcell": ["CD4+ T cell", "CD8+ T cell", "T cell"],
            "immune": ["CD4+ T cell", "CD8+ T cell", "T cell", "B cell", "Dendritic cell", "Myeloid cell", "Mast cell"]}

CONTACT_R = 50.0      # um; counts of reference cells, tumor fraction
VASC_R2 = 100.0       # um; second endothelial count
CROWD_R = 30.0        # um; n_cells_30um
K_CROWD = 10          # nearest other cells for the crowding score
KERNEL_R = 100.0      # um; cutoff of the Gaussian kernels
SIGMA = 50.0          # um; kernel width
BLOCK = 200.0         # um; spatial blocks for bootstrap / cross-validation folds
CHUNK = 20_000

# Label-free neighbour-expression kernels (markers present on the panel; stromal is weak: no collagens / DCN / LUM)
MARKERS = {"vascular": ["PECAM1", "CDH5", "CLDN5", "KDR", "ENG"],
           "cytotoxic": ["CD3E", "CD2", "CD8A", "GZMA", "GZMB", "GZMK", "PRF1"],
           "myeloid": ["CD68", "CD163", "CD14", "ITGAM"],
           "stromal": ["PDGFRB", "FAP"]}
GZM = ["GZMA", "GZMB", "GZMK"]
HSP11 = ["HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPH1", "DNAJB1", "HSPB1", "HSPD1", "HSPE1"]
HYPOXIA_TF8 = ["HIF1A", "EPAS1", "NFE2L2", "CREB1", "RELA", "RELB", "NFKB1", "NFKB2"]

# Axis catalogue: name -> (tier, group kind, direction note). "count" axes get Away/Close groups, the others
# pooled tertiles. Direction: what the High tertile (or Close) means for the tumor cell.
AXES = {
    # primary tier (expression-free)
    "vascular_dist": ("primary", "tertile", "far from vessels (hypoxia proxy)"),
    "vascular_n100": ("primary", "count", "vessel-rich within 100 um"),
    "crowding": ("primary", "tertile", "crowded (inverse mean distance to 10 nearest cells)"),
    "n_cells_30um": ("primary", "tertile", "many cells within 30 um"),
    "tumor_frac": ("primary", "tertile", "tumor core (labelled cells within 50 um)"),
    "nontumor_dist": ("primary", "tertile", "far from any labelled non-tumor cell"),
    "tcell_n50": ("primary", "count", "T cells within 50 um"),
    "tcell_dist": ("primary", "tertile", "far from T cells"),
    "tcell_gzm_kernel": ("primary", "tertile", "GZM-weighted T-cell attack (label-restricted)"),
    "immune_n50": ("primary", "count", "immune cells within 50 um"),
    "myeloid_n50": ("primary", "count", "myeloid cells within 50 um"),
    "myeloid_dist": ("primary", "tertile", "far from myeloid cells"),
    "fibroblast_n50": ("primary", "count", "fibroblasts within 50 um"),
    "fibroblast_dist": ("primary", "tertile", "far from fibroblasts"),
    "cell_area": ("primary", "tertile", "large cell"),
    "nucleus_ratio": ("primary", "tertile", "large nucleus-to-cell area ratio"),
    "k_vascular": ("primary", "tertile", "label-free vascular marker kernel"),
    "k_cytotoxic": ("primary", "tertile", "label-free cytotoxic marker kernel"),
    "k_myeloid": ("primary", "tertile", "label-free myeloid marker kernel"),
    "k_stromal": ("primary", "tertile", "label-free stromal marker kernel (weak)"),
    # secondary tier (whole-cell expression of the tumor cell)
    "hypoxia": ("secondary", "tertile", "HALLMARK_HYPOXIA score"),
    "hypoxia_tf": ("secondary", "tertile", "old 8-TF rule (2nd highest TF)"),
    "heat_shock": ("secondary", "tertile", "11 HSP genes"),
    "tnfa_nfkb": ("secondary", "tertile", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    "ifng": ("secondary", "tertile", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    "sg_markers": ("secondary", "tertile", "SG marker genes (fraction in SGs > 0.4)"),
    "granule_count": ("secondary", "tertile", "granules per cell (log1p)"),
    "granule_sg_transcripts": ("secondary", "tertile", "SG-marker transcripts in granules per cell (log1p)"),
}
GROUP_NAMES = {"vascular_dist": ("Near", "Mid", "Far"), "tumor_frac": ("Boundary", "Mid", "Core")}

# Labels that are spatial by construction: excluded from the axes that define them (analysis B)
LABEL_AXIS_EXCLUSIONS = {"Stroma border": ["fibroblast_n50", "fibroblast_dist", "k_stromal"],
                         "Myeloid border": ["myeloid_n50", "myeloid_dist", "k_myeloid"]}


# ----------------------------------------------------------------------------------------------------------------
# Spatial axes (one sample)
# ----------------------------------------------------------------------------------------------------------------
def _radius_counts(xy_ref, xy_query, r):
    """Number of reference cells within r um of each query cell."""
    tree = KDTree(xy_ref)
    out = np.zeros(len(xy_query), dtype=np.int64)
    for s in range(0, len(xy_query), CHUNK):
        out[s:s + CHUNK] = tree.query_radius(xy_query[s:s + CHUNK], r=r, count_only=True)
    return out


def label_axes(obs_all, tumor_idx):
    """Label-based spatial axes of the tumor cells (rows = tumor_idx positions of obs_all): per reference set the
    1-NN distance and the count within CONTACT_R (endothelial also within VASC_R2); crowding and n_cells_30um
    over all cells (self excluded); tumor fraction and distance to the nearest non-tumor cell over labelled cells
    (Unknown / Mixed dropped). A reference set absent from the sample gives count 0 and distance NaN."""
    xy_all = obs_all[["global_x", "global_y"]].values.astype(np.float64)
    types = obs_all["cell_type_merged"].astype(str).values
    xy_t = xy_all[tumor_idx]
    out = pd.DataFrame(index=obs_all.index[tumor_idx])
    for name, cts in REF_SETS.items():
        ref = np.isin(types, cts)
        if ref.sum() == 0:
            print(f"  {name}: no reference cell in this sample -> count 0, distance NaN")
            out[f"{name}_dist"] = np.nan
            out[f"{name}_n50"] = 0
            if name == "endothelial":
                out["endothelial_n100"] = 0
            continue
        out[f"{name}_dist"] = NearestNeighbors(n_neighbors=1).fit(xy_all[ref]).kneighbors(xy_t)[0].ravel()
        out[f"{name}_n50"] = _radius_counts(xy_all[ref], xy_t, CONTACT_R)
        if name == "endothelial":
            out["endothelial_n100"] = _radius_counts(xy_all[ref], xy_t, VASC_R2)
    d = NearestNeighbors(n_neighbors=K_CROWD + 1).fit(xy_all).kneighbors(xy_t)[0][:, 1:]   # drop self
    out["crowding"] = 1.0 / np.maximum(d.mean(axis=1), 1e-6)
    out["n_cells_30um"] = _radius_counts(xy_all, xy_t, CROWD_R) - 1
    labelled = ~np.isin(types, EXCLUDED_TYPES)
    tumor = types == TUMOR
    n_lab = _radius_counts(xy_all[labelled], xy_t, CONTACT_R) - 1
    n_tum = _radius_counts(xy_all[tumor], xy_t, CONTACT_R) - 1
    out["tumor_frac"] = np.where(n_lab > 0, n_tum / np.maximum(n_lab, 1), np.nan)
    non = labelled & ~tumor
    out["nontumor_dist"] = NearestNeighbors(n_neighbors=1).fit(xy_all[non]).kneighbors(xy_t)[0].ravel()
    return out


def kernel_score(xy_query, xy_src, m_src):
    """log1p of the Gaussian-kernel (SIGMA, cutoff KERNEL_R) weighted sum of m_src over the source cells around
    each query cell, and the number of source cells within the cutoff. Abundance times expression."""
    tree = KDTree(xy_src)
    score = np.zeros(len(xy_query))
    n = np.zeros(len(xy_query), dtype=np.int64)
    for s in range(0, len(xy_query), CHUNK):
        ind, dist = tree.query_radius(xy_query[s:s + CHUNK], r=KERNEL_R, return_distance=True)
        lens = np.array([len(i) for i in ind])
        if lens.sum() == 0:
            continue
        rows = np.repeat(np.arange(len(ind)), lens)
        ii = np.concatenate(ind)
        dd = np.concatenate(dist)
        w = np.exp(-(dd ** 2) / (2.0 * SIGMA ** 2))
        score[s:s + CHUNK] = np.bincount(rows, weights=w * m_src[ii], minlength=len(ind))
        n[s:s + CHUNK] = lens
    return np.log1p(score), n


def marker_mean(adata_lognorm, genes):
    """Mean log-normalised expression of the genes present on the panel (all cells)."""
    present = [g for g in genes if g in adata_lognorm.var_names]
    X = adata_lognorm[:, present].X
    return np.asarray(X.mean(axis=1)).ravel()


def kernel_axes(adata_lognorm, tumor_idx):
    """The label-restricted GZM T-cell kernel and the four label-free marker kernels over non-tumor cells."""
    xy_all = adata_lognorm.obs[["global_x", "global_y"]].values.astype(np.float64)
    types = adata_lognorm.obs["cell_type_merged"].astype(str).values
    xy_t = xy_all[tumor_idx]
    out = pd.DataFrame(index=adata_lognorm.obs.index[tumor_idx])
    tcell = np.isin(types, REF_SETS["tcell"])
    if tcell.sum() > 0:
        m = marker_mean(adata_lognorm[tcell], GZM)
        out["tcell_gzm_kernel"], _ = kernel_score(xy_t, xy_all[tcell], m)
    else:
        out["tcell_gzm_kernel"] = 0.0
    non = types != TUMOR
    for name, genes in MARKERS.items():
        m = marker_mean(adata_lognorm[non], genes)
        out[f"k_{name}"], n_non = kernel_score(xy_t, xy_all[non], m)
    out["nontumor_n100"] = n_non
    return out


# ----------------------------------------------------------------------------------------------------------------
# Expression axes (tumor cells, whole-cell log-normalised)
# ----------------------------------------------------------------------------------------------------------------
def read_gmt(path):
    sets = {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def expression_axes(ad_tumor, gene_sets):
    """sc.tl.score_genes for every signature (canonical gene list as it stands on the panel) plus the 8-TF
    hypoxia rule (second-highest log expression of the eight TFs)."""
    out = pd.DataFrame(index=ad_tumor.obs.index)
    for name, genes in gene_sets.items():
        glist = [g for g in genes if g in ad_tumor.var_names]
        sc.tl.score_genes(ad_tumor, glist, score_name="_tmp", ctrl_size=50, random_state=0)
        out[name] = ad_tumor.obs["_tmp"].values
    del ad_tumor.obs["_tmp"]
    tf = [g for g in HYPOXIA_TF8 if g in ad_tumor.var_names]
    X = ad_tumor[:, tf].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    out["hypoxia_tf"] = np.sort(X, axis=1)[:, -2]
    return out


def granule_axes(granule_adata, tumor_cell_ids, sg_genes):
    """Granules per tumor cell and SG-marker transcripts inside its granules (both log1p; 0 without granules)."""
    cid = granule_adata.obs["cell_id"].astype(str).values
    present = [g for g in sg_genes if g in granule_adata.var_names]
    sg_tx = np.asarray(granule_adata[:, present].X.sum(axis=1)).ravel()
    per_cell = pd.DataFrame({"cell_id": cid, "n": 1, "sg": sg_tx}).groupby("cell_id").sum()
    idx = pd.Index(tumor_cell_ids.astype(str))
    n = per_cell["n"].reindex(idx).fillna(0).values
    sg = per_cell["sg"].reindex(idx).fillna(0).values
    return pd.DataFrame({"granule_count": np.log1p(n), "granule_sg_transcripts": np.log1p(sg)})


# ----------------------------------------------------------------------------------------------------------------
# Pooled groups and blocks
# ----------------------------------------------------------------------------------------------------------------
def pooled_groups(values, kind, names=("Low", "Mid", "High")):
    """Tertile groups on the pooled cells (rank-based, NaN -> 'NA') or Away/Close for counts."""
    v = np.asarray(values, dtype=float)
    out = np.full(len(v), "NA", dtype=object)
    if kind == "count":
        ok = np.isfinite(v)
        out[ok] = np.where(v[ok] > 0, "Close", "Away")
        return out
    ok = np.isfinite(v)
    r = pd.Series(v[ok]).rank(method="first").values
    tert = np.minimum((r - 1) * 3 // ok.sum(), 2).astype(int)
    out[ok] = np.array(names, dtype=object)[tert]
    return out


def spatial_blocks(sample, x, y):
    bx = np.floor(np.asarray(x) / BLOCK).astype(int)
    by = np.floor(np.asarray(y) / BLOCK).astype(int)
    return pd.Series(sample).astype(str).values + "_" + bx.astype(str) + "_" + by.astype(str)


def block_resample(blocks, rng):
    """Index array of one spatial block bootstrap resample: blocks drawn with replacement within each sample."""
    codes, uniq = pd.factorize(blocks)
    sample_of_block = pd.Series(uniq).str.split("_").str[0].values
    members = pd.Series(np.arange(len(codes))).groupby(codes).apply(np.asarray)
    idx = []
    for s in np.unique(sample_of_block):
        bl = np.where(sample_of_block == s)[0]
        draw = rng.choice(bl, size=len(bl), replace=True)
        idx.extend(members[b] for b in draw)
    return np.concatenate(idx)


# ----------------------------------------------------------------------------------------------------------------
# Analysis A: stress profile of the labels
# ----------------------------------------------------------------------------------------------------------------
def dummies(codes):
    n_cat = codes.max() + 1
    D = np.zeros((len(codes), n_cat), dtype=np.float32)
    D[np.arange(len(codes)), codes] = 1.0
    return D


def r2_matrix(Y, X):
    """R2 of every column of Y on the design X (intercept included in X)."""
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    ss_res = (resid ** 2).sum(axis=0)
    ss_tot = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
    return 1.0 - ss_res / np.maximum(ss_tot, 1e-12)


def r2_beyond_sample(Y, sample_codes, label_codes):
    """R2(Y ~ sample + label) - R2(Y ~ sample), every column of Y at once."""
    Ds = dummies(sample_codes)
    Dl = dummies(label_codes)[:, 1:]          # drop one level, Ds spans the intercept
    return r2_matrix(Y, np.hstack([Ds, Dl])) - r2_matrix(Y, Ds)


def epsilon2(y, label_codes):
    """Kruskal-Wallis epsilon squared (rank-based effect size of the label on y)."""
    ok = np.isfinite(y)
    groups = [y[ok][label_codes[ok] == k] for k in np.unique(label_codes[ok])]
    H = stats.kruskal(*groups).statistic
    n = ok.sum()
    return (H - len(groups) + 1) / (n - len(groups))


# ----------------------------------------------------------------------------------------------------------------
# Analysis B: correspondence of stress groups and labels
# ----------------------------------------------------------------------------------------------------------------
def cramers_v(a, b):
    tab = pd.crosstab(a, b).values
    chi2 = stats.chi2_contingency(tab, correction=False)[0]
    n = tab.sum()
    return np.sqrt(chi2 / (n * (min(tab.shape) - 1)))


def ami(a, b):
    return adjusted_mutual_info_score(a, b)


def log2_or_table(label_codes, high, n_labels):
    """log2 odds ratio of being High (vs Low) inside each label vs outside, with 0.5 continuity correction."""
    a = np.bincount(label_codes[high], minlength=n_labels).astype(float)        # High in label
    b = np.bincount(label_codes[~high], minlength=n_labels).astype(float)       # Low in label
    c = high.sum() - a                                                            # High outside
    d = (~high).sum() - b                                                         # Low outside
    return np.log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))


def mantel_haenszel_log2_or(label_codes, high, sample_codes, n_labels):
    """Sample-stratified (Mantel-Haenszel) log2 OR per label."""
    out = np.full(n_labels, np.nan)
    for k in range(n_labels):
        num = den = 0.0
        for s in np.unique(sample_codes):
            m = sample_codes == s
            n = m.sum()
            a = np.sum(m & (label_codes == k) & high); b = np.sum(m & (label_codes == k) & ~high)
            c = np.sum(m & (label_codes != k) & high); d = np.sum(m & (label_codes != k) & ~high)
            num += a * d / n
            den += b * c / n
        out[k] = np.log2((num + 1e-9) / (den + 1e-9))
    return out


# ----------------------------------------------------------------------------------------------------------------
# Analysis C: information content of the PCs
# ----------------------------------------------------------------------------------------------------------------
def standardize(X):
    mu, sd = X.mean(axis=0), X.std(axis=0)
    return (X - mu) / np.maximum(sd, 1e-8)


def cv_predictions(X, y, folds):
    """Out-of-fold ridge predictions (alpha = 1) for every cell. folds: integer fold id per cell (spatial-block
    folds or the sample id for leave-one-sample-out). Returns the prediction array (NaN where y is NaN)."""
    pred = np.full(len(y), np.nan)
    ok = np.isfinite(y)
    for f in np.unique(folds):
        test = (folds == f) & ok
        train = (folds != f) & ok
        model = Ridge(alpha=1.0).fit(X[train], y[train])
        pred[test] = model.predict(X[test])
    return pred


def r2_of(y, pred, mask=None):
    m = np.isfinite(y) & np.isfinite(pred)
    if mask is not None:
        m &= mask
    if m.sum() < 50:
        return np.nan
    return 1.0 - ((y[m] - pred[m]) ** 2).sum() / max(((y[m] - y[m].mean()) ** 2).sum(), 1e-12)


def auroc_of(y, pred, mask=None):
    m = np.isfinite(y) & np.isfinite(pred)
    if mask is not None:
        m &= mask
    if m.sum() < 50 or len(np.unique(y[m])) < 2:
        return np.nan
    return roc_auc_score(y[m], pred[m])


def block_folds(blocks, n_folds, rng):
    """Fold id per cell from a random assignment of spatial blocks to folds (GroupKFold on shuffled blocks)."""
    codes, uniq = pd.factorize(blocks)
    perm = rng.permutation(len(uniq))
    fold_of_block = np.empty(len(uniq), dtype=int)
    fold_of_block[perm] = np.arange(len(uniq)) % n_folds
    return fold_of_block[codes]


def morans_i(z, W):
    """Moran's I of z on the (row-standardised) graph W; z is centred first. NaN cells are dropped."""
    ok = np.isfinite(z)
    if ok.sum() < len(z):
        W = W[ok][:, ok]
        z = z[ok]
    z = z - z.mean()
    rs = np.asarray(W.sum(axis=1)).ravel()
    Wn = sparse.diags(np.where(rs > 0, 1.0 / np.maximum(rs, 1e-12), 0.0)) @ W
    S0 = Wn.sum()
    return (len(z) / S0) * (z @ (Wn @ z)) / max(z @ z, 1e-12)


# ----------------------------------------------------------------------------------------------------------------
# Analysis D: program scores vs axes
# ----------------------------------------------------------------------------------------------------------------
def rank_cols(M):
    return np.apply_along_axis(lambda v: stats.rankdata(v, nan_policy="omit"), 0, M)


def residualize_on_sample(M, sample_codes):
    out = M.astype(float).copy()
    for s in np.unique(sample_codes):
        m = sample_codes == s
        out[m] -= np.nanmean(out[m], axis=0)
    return out


def spearman_matrix(A, B, sample_codes=None):
    """Spearman correlations between the columns of A and the columns of B (pooled; sample-adjusted partial
    version when sample_codes is given: ranks residualised on the sample means)."""
    RA, RB = rank_cols(A), rank_cols(B)
    if sample_codes is not None:
        RA, RB = residualize_on_sample(RA, sample_codes), residualize_on_sample(RB, sample_codes)
    RA = np.nan_to_num(RA - np.nanmean(RA, axis=0))
    RB = np.nan_to_num(RB - np.nanmean(RB, axis=0))
    num = RA.T @ RB
    den = np.sqrt((RA ** 2).sum(axis=0))[:, None] * np.sqrt((RB ** 2).sum(axis=0))[None, :]
    return num / np.maximum(den, 1e-12)


# ----------------------------------------------------------------------------------------------------------------
# Composite stress scores and explicit stress subtypes (approach 2)
# ----------------------------------------------------------------------------------------------------------------
EXPR_SIG = ["hypoxia", "heat_shock", "tnfa_nfkb", "ifng", "sg_markers"]     # whole-cell signature scores
EXPR_OTHER = ["hypoxia_tf", "granule_count", "granule_sg_transcripts"]      # rule score and granule burden


def axis_set(ax):
    """'env' for the 20 expression-free axes, 'signature' for the 5 whole-cell signature scores, 'other' for the
    8-TF rule and the two granule-burden axes."""
    if AXES[ax][0] == "primary":
        return "env"
    if ax in EXPR_SIG:
        return "signature"
    assert ax in EXPR_OTHER, ax
    return "other"


# Environmental composites: member axis -> sign (+1: the axis increases with the stress)
COMPOSITES_ENV = {
    "Vascular deprivation": {"vascular_dist": 1, "vascular_n100": -1, "k_vascular": -1},
    "Crowding": {"crowding": 1, "n_cells_30um": 1},
    "Core isolation": {"tumor_frac": 1, "nontumor_dist": 1},
    "Immune attack": {"tcell_n50": 1, "tcell_gzm_kernel": 1, "immune_n50": 1, "k_cytotoxic": 1, "tcell_dist": -1},
    "Myeloid contact": {"myeloid_n50": 1, "k_myeloid": 1, "myeloid_dist": -1},
    "Stromal contact": {"fibroblast_n50": 1, "k_stromal": 1, "fibroblast_dist": -1},
}

# Expression composites: the four signature programmes. The SG-marker score and the granule-burden axes are
# stress axes (approach 1) but not subtyping features (user decision 2026-11-08).
COMPOSITES_EXPR = {
    "Hypoxia (expr)": {"hypoxia": 1},
    "Heat shock": {"heat_shock": 1},
    "Inflammatory (TNFA/NFKB)": {"tnfa_nfkb": 1},
    "Interferon (IFNG)": {"ifng": 1},
}


ARCHETYPE_THR = 0.5     # standardised composite score above which a cell is assigned to its dominant archetype


def composite_scores(stress, composites):
    """Mean of the sign-aligned pooled z of the member axes (NaN ignored, remaining NaN -> median), re-standardised."""
    out = {}
    for name, members in composites.items():
        M = np.column_stack([sign * stress[f"z_{ax}"].values.astype(float) for ax, sign in members.items()])
        v = np.nanmean(M, axis=1)
        v = np.where(np.isfinite(v), v, np.nanmedian(v))
        out[name] = (v - v.mean()) / v.std()
    return pd.DataFrame(out)


def dominant_archetype(C, thr=ARCHETYPE_THR):
    """Name of the composite with the largest standardised score, or 'Unstressed' when none exceeds thr."""
    M = C.values
    names = np.array(C.columns, dtype=object)
    best = M.argmax(axis=1)
    return np.where(M.max(axis=1) > thr, names[best], "Unstressed")


def enrichment_table(a_codes, n_a, b_codes, n_b):
    """Contingency counts and log2(observed / expected under independence), 0.5 continuity correction."""
    tab = np.zeros((n_a, n_b))
    np.add.at(tab, (a_codes, b_codes), 1.0)
    n = max(tab.sum(), 1.0)
    exp = tab.sum(axis=1, keepdims=True) * tab.sum(axis=0, keepdims=True) / n
    return tab, np.log2((tab + 0.5) / (exp + 0.5))


def ami_from_table(tab):
    """Adjusted mutual information from a contingency table (same value as sklearn on the labels)."""
    from sklearn.metrics.cluster import mutual_info_score
    try:
        from sklearn.metrics.cluster._expected_mutual_info_fast import expected_mutual_information
    except ImportError:                                   # older / newer sklearn layouts
        from sklearn.metrics.cluster.supervised import expected_mutual_information
    tab = np.asarray(tab)
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    n = tab.sum()
    mi = mutual_info_score(None, None, contingency=tab)
    emi = expected_mutual_information(tab.astype(np.int64), int(n))
    pa, pb = tab.sum(axis=1) / n, tab.sum(axis=0) / n
    ha, hb = -(pa * np.log(pa)).sum(), -(pb * np.log(pb)).sum()
    denom = 0.5 * (ha + hb) - emi
    return (mi - emi) / denom if abs(denom) > 1e-12 else 0.0


def cramers_v_from_table(tab):
    tab = np.asarray(tab)
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    chi2 = stats.chi2_contingency(tab, correction=False)[0]
    return np.sqrt(chi2 / (tab.sum() * (min(tab.shape) - 1)))


def label_correspondence(a, b, boot_idx, sample_codes, sample_names, min_abs=0.5, min_same=3):
    """Correspondence of two labellings of the same cells (a = substates, b = stress subtypes): contingency
    table, per-cell log2 enrichment with block-bootstrap CI and per-sample sign agreement, Cramer's V and AMI
    with bootstrap CIs. Returns (contingency DataFrame, long DataFrame, dict of summary statistics)."""
    ac, al = pd.factorize(a)
    bc, bl = pd.factorize(b)
    al, bl = list(al), list(bl)
    tab, enr = enrichment_table(ac, len(al), bc, len(bl))
    boots, v_b, ami_b = [], [], []
    for idx in boot_idx:
        t, e = enrichment_table(ac[idx], len(al), bc[idx], len(bl))
        boots.append(e)
        v_b.append(cramers_v_from_table(t))
        ami_b.append(ami_from_table(t))
    boots = np.array(boots)
    lo, hi = np.percentile(boots, 2.5, axis=0), np.percentile(boots, 97.5, axis=0)
    per_s = {}
    for si, sname in enumerate(sample_names):
        m = sample_codes == si
        per_s[sname] = enrichment_table(ac[m], len(al), bc[m], len(bl))[1] if m.sum() > 0 else np.full_like(enr, np.nan)
    rows = []
    row_tot, col_tot = tab.sum(axis=1), tab.sum(axis=0)
    for i, la in enumerate(al):
        for j, lb in enumerate(bl):
            signs = [np.sign(per_s[s][i, j]) for s in sample_names if np.isfinite(per_s[s][i, j]) and per_s[s][i, j] != 0]
            same = int(np.sum(np.array(signs) == np.sign(enr[i, j])))
            rows.append({"substate": la, "stress_subtype": lb, "n": int(tab[i, j]),
                         "pct_of_substate": 100 * tab[i, j] / max(row_tot[i], 1), "pct_of_stress_subtype": 100 * tab[i, j] / max(col_tot[j], 1),
                         "log2_enrichment": enr[i, j], "ci_low": lo[i, j], "ci_high": hi[i, j],
                         **{f"log2_enr_{s}": per_s[s][i, j] for s in sample_names},
                         "n_samples_same_sign": same,
                         "reproducible": bool(abs(enr[i, j]) > min_abs and same >= min_same and (lo[i, j] > 0 or hi[i, j] < 0))})
    cont = pd.DataFrame(tab.astype(int), index=al, columns=bl)
    summary = {"cramers_v": cramers_v_from_table(tab), "v_ci_low": np.percentile(v_b, 2.5), "v_ci_high": np.percentile(v_b, 97.5),
               "ami": ami_from_table(tab), "ami_ci_low": np.percentile(ami_b, 2.5), "ami_ci_high": np.percentile(ami_b, 97.5),
               "n_reproducible": int(sum(r["reproducible"] for r in rows)), "v_boot": np.array(v_b), "ami_boot": np.array(ami_b)}
    return cont, pd.DataFrame(rows), summary


def within_sample_z(stress, axes, sample_codes):
    """Each axis standardised within sample (NaN kept)."""
    Z = stress[[f"z_{a}" for a in axes]].values.astype(np.float64).copy()
    for s in np.unique(sample_codes):
        m = sample_codes == s
        mu, sd = np.nanmean(Z[m], axis=0), np.nanstd(Z[m], axis=0)
        Z[m] = (Z[m] - mu) / np.maximum(sd, 1e-8)
    return Z
