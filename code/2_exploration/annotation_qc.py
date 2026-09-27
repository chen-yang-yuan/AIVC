"""Step 0b: annotation QC report. Marker-set scores for all cells, class profiles, tumor-set purity flags,
profiles of the manually annotated malignant clusters, and a cross-dataset transfer of the non-epithelial labels
(classifier trained on the three 10x-labelled datasets). Report only: nothing here changes the tumor definition."""

import pickle

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.linear_model import LogisticRegression

import config as C
import io_utils as IO


def load_all_cells(ds):
    """All cells: log1p-CPM sparse matrix over the panel plus the label columns."""
    a = sc.read_h5ad(IO.hydrate(f"{C.DATA_DIR}{ds}/intermediate_data/adata.h5ad"))
    X = sparse.csr_matrix(a.X, dtype=np.float32)
    depth = np.asarray(X.sum(axis=1)).ravel()
    Y = sparse.diags(np.where(depth > 0, 1e4 / np.maximum(depth, 1), 0.0)) @ X
    Y = sparse.csr_matrix(Y); Y.data = np.log1p(Y.data)
    obs = a.obs.copy()
    for col in ("cell_type", "cluster_labels", "cell_type_merged", "segmentation_method"):
        if col in obs:
            obs[col] = obs[col].astype(str)
    obs["cell_id"] = obs["cell_id"].astype(str)
    return Y, obs.reset_index(drop=True), np.asarray(a.var_names).astype(str)


def marker_scores(Y, genes, ds):
    """Mean log1p-CPM over each marker set (missing genes dropped); lineage = this cancer type, off_lineage = others."""
    pos = {g: i for i, g in enumerate(genes)}
    sets = dict(C.MARKER_SETS)
    sets["lineage"] = C.LINEAGE_MARKERS[ds]
    others = [g for d, gl in C.LINEAGE_MARKERS.items() if d != ds for g in gl if g not in C.LINEAGE_MARKERS[ds]]
    sets["off_lineage"] = sorted(set(others))
    out = {}
    for name, gl in sets.items():
        idx = [pos[g] for g in gl if g in pos]
        out[name] = np.asarray(Y[:, idx].mean(axis=1)).ravel() if idx else np.full(Y.shape[0], np.nan)
    df = pd.DataFrame(out)
    df.attrs["sets"] = {k: [g for g in v if g in pos] for k, v in sets.items()}
    return df


def provenance(ds, obs):
    manual = ds in C.MANUAL_ANNOTATION
    fine = obs["cluster_labels"] if manual else obs["cell_type"]
    ct = obs["cell_type_merged"]
    n_mal_clusters = None
    if manual:
        m = pickle.load(open(IO.hydrate(C.UTILS_DIR + "cell_type_dict_manual.pkl"), "rb"))
        n_mal_clusters = len(m[ds]["Malignant cell"])
    return {"dataset": C.short(ds), "source": "manual graph-cluster annotation" if manual else "10x supervised labels",
            "n_fine_labels": int(fine.nunique()), "n_categories_present": int(ct.nunique()),
            "n_cells": int(len(obs)), "tumor_frac": float((ct == C.TUMOR_TYPE).mean()),
            "unresolved_frac": float(ct.isin(["Unknown", "Mixed"]).mean()),
            "has_nonmalignant_epithelium": bool((ct == "Epithelial cell (non-malignant)").any()),
            "n_malignant_clusters": n_mal_clusters if manual else 1}


def class_profile(scores, obs):
    g = scores.groupby(obs["cell_type_merged"].values)
    mean = g.mean(); p10 = g.quantile(0.10); p90 = g.quantile(0.90); n = g.size()
    out = mean.add_suffix("_mean").join(p10.add_suffix("_p10")).join(p90.add_suffix("_p90"))
    out.insert(0, "n_cells", n)
    return out.reset_index().rename(columns={"index": "cell_type_merged"})


def purity_flags(scores, obs, ds):
    ct = obs["cell_type_merged"].values
    mal = ct == C.TUMOR_TYPE
    stroma = np.isin(ct, C.STROMAL_TYPES)
    immune = np.isin(ct, C.IMMUNE_TYPES)
    s = scores
    q = lambda col, mask, p: float(np.nanpercentile(s.loc[mask, col], p)) if mask.sum() > 0 else np.nan
    lin_neg = s["lineage"].values < q("lineage", mal, 5)
    if C.EPITHELIAL_TUMOR.get(ds, True):
        lin_neg &= s["epithelial"].values < q("epithelial", mal, 5)
    flags = pd.DataFrame({
        "lineage_negative": lin_neg,
        "immune_like": s["immune"].values > q("immune", stroma, 95),
        "endothelial_like": s["endothelial"].values > q("endothelial", immune | (ct == "Fibroblast (CAF)"), 95),
        "off_lineage_high": s["off_lineage"].values > q("off_lineage", stroma, 95),
        "proliferating": s["proliferation"].values > q("proliferation", stroma, 95),
    })
    summary = {"dataset": C.short(ds), "n_malignant": int(mal.sum())}
    for c in flags.columns:
        summary[f"frac_{c}"] = float(flags.loc[mal, c].mean())
    return flags, summary


def tumor_cluster_profiles(scores, obs, ds):
    """Manual datasets: score profile of every graph cluster called malignant, with z-scores across those clusters."""
    if ds not in C.MANUAL_ANNOTATION:
        return pd.DataFrame()
    mal = obs["cell_type_merged"].values == C.TUMOR_TYPE
    cl = obs.loc[mal, "cluster_labels"].values
    g = scores[mal].groupby(cl)
    prof = g.mean(); prof.insert(0, "n_cells", g.size())
    cols = [c for c in prof.columns if c != "n_cells"]
    w = prof["n_cells"].values / prof["n_cells"].sum()
    mu = (prof[cols].values * w[:, None]).sum(0)
    sd = np.sqrt(((prof[cols].values - mu) ** 2 * w[:, None]).sum(0)) + 1e-9
    z = pd.DataFrame((prof[cols].values - mu) / sd, index=prof.index, columns=[f"z_{c}" for c in cols])
    prof = prof.join(z)
    flag = np.zeros(len(prof), dtype=bool)
    if ds == "Xenium_5K_Prostate":
        flag = (prof["z_basal"] > 1) | (prof["z_lineage"] < -1)
    elif ds == "Xenium_5K_Skin":
        flag = prof["z_lineage"] < -1
    elif ds == "Xenium_5K_LC":
        flag = (prof["z_lineage"] < -1) | (prof["z_epithelial"] < -1)
    prof["benign_like"] = np.asarray(flag)
    prof["contaminated"] = (prof["z_immune"] > 1) | (prof["z_endothelial"] > 1)
    return prof.reset_index().rename(columns={"index": "cluster"})


# ------------------------------------------------------------------ cross-dataset transfer of non-epithelial labels
def _training_set(dss, max_per_class_ds=None, seed=None, exclude_ds=None):
    max_per = max_per_class_ds or C.TRANSFER_MAX_PER_CLASS_DS
    rng = np.random.default_rng(C.SEED if seed is None else seed)
    Xs, ys, dsl = [], [], []
    for ds in dss:
        if exclude_ds and ds == exclude_ds:
            continue
        Y, obs, genes = load_all_cells(ds)
        ct = obs["cell_type_merged"].values.copy()
        ct[np.isin(ct, [C.TUMOR_TYPE, "Epithelial cell (non-malignant)"])] = C.TRANSFER_EPITHELIAL_CLASS
        keep_idx = []
        for c in np.unique(ct):
            if c in C.TRANSFER_EXCLUDE:
                continue
            idx = np.where(ct == c)[0]
            if len(idx) > max_per:
                idx = rng.choice(idx, max_per, replace=False)
            keep_idx.append(idx)
        idx = np.sort(np.concatenate(keep_idx))
        Xs.append(Y[idx]); ys.append(ct[idx]); dsl.append(np.full(len(idx), C.short(ds)))
    return sparse.vstack(Xs).tocsr(), np.concatenate(ys), np.concatenate(dsl)


def fit_transfer(dss=None, seed=None):
    dss = dss or C.TRANSFER_TRAIN
    X, y, _ = _training_set(dss, seed=seed)
    clf = LogisticRegression(C=0.05, max_iter=300, solver="saga", tol=1e-3, class_weight="balanced")
    clf.fit(X, y)
    return clf


def lodo_reference(dss=None, seed=None):
    """Leave-one-dataset-out accuracy within the training datasets (per class), the reference for the transfer."""
    dss = dss or C.TRANSFER_TRAIN
    rows = []
    for held in dss:
        X, y, _ = _training_set(dss, seed=seed, exclude_ds=held, max_per_class_ds=3000)
        clf = LogisticRegression(C=0.05, max_iter=200, solver="saga", tol=1e-3, class_weight="balanced").fit(X, y)
        Xt, yt, _ = _training_set([held], seed=seed, max_per_class_ds=3000)
        pred = clf.predict(Xt)
        for c in np.unique(yt):
            m = yt == c
            rows.append({"held_out": C.short(held), "class": c, "n": int(m.sum()), "accuracy": float((pred[m] == c).mean())})
    return pd.DataFrame(rows)


def apply_transfer(clf, Y, obs, ds, conf=None):
    conf = conf or C.TRANSFER_CONF
    ct = obs["cell_type_merged"].values
    P = clf.predict_proba(Y)
    pred = clf.classes_[P.argmax(1)]; pmax = P.max(1)
    rows = []
    for c in np.unique(ct):
        m = ct == c
        target = C.TRANSFER_EPITHELIAL_CLASS if c in (C.TUMOR_TYPE, "Epithelial cell (non-malignant)") else c
        r = {"dataset": C.short(ds), "label": c, "n": int(m.sum()),
             "agreement": float((pred[m] == target).mean()) if target in clf.classes_ else np.nan,
             "frac_confident": float((pmax[m] >= conf).mean())}
        top = pd.Series(pred[m][pmax[m] >= conf]).value_counts(normalize=True)
        r["top_prediction"] = top.index[0] if len(top) else None
        r["top_prediction_frac"] = float(top.iloc[0]) if len(top) else np.nan
        for grp, types in (("immune", C.IMMUNE_TYPES), ("stromal", C.STROMAL_TYPES), ("endothelial", ["Endothelial cell", "Lymphatic endothelial cell"])):
            r[f"frac_confident_{grp}"] = float((np.isin(pred[m], types) & (pmax[m] >= conf)).mean())
        rows.append(r)
    return pd.DataFrame(rows), pred, pmax
