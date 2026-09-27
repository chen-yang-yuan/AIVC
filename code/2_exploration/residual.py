"""Step 3b: cytoplasmic structure beyond nuclear identity.

Residualise log-normalised cytoplasmic expression (half 1) on a nuclear design (nuclear half-1 PCs with an
errors-in-variables correction, nuclear-cluster one-hot, depth spline, morphology and segmentation covariates);
test the residual scores and residual clusters for spatial coherence and boundary/immune association under nulls
that keep the nuclear subtype (stratified permutation), an empirical leakage floor (nuclear half 2 residualised the
same way) and tile block bootstraps; apply the pre-registered eligibility rule for the Fig. 1 panel.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import SplineTransformer

import config as C
import compartment as CP
import association as AS


# ------------------------------------------------------------------ dense log-normalised matrices
def lognorm_dense(X, gene_mask, target_sum):
    X = sparse.csr_matrix(X, dtype=np.float32)[:, gene_mask]
    depth = np.asarray(X.sum(axis=1)).ravel()
    Y = sparse.diags(np.where(depth > 0, target_sum / np.maximum(depth, 1), 0.0)) @ X
    Y = sparse.csr_matrix(Y)
    Y.data = np.log1p(Y.data)
    return np.asarray(Y.todense(), dtype=np.float32)


def fit_pca(Y, n_pcs, seed=None):
    seed = C.SEED if seed is None else seed
    n_pcs = min(n_pcs, Y.shape[1] - 1, Y.shape[0] - 1)
    p = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed).fit(Y)
    return p.transform(Y).astype(np.float32), p.components_.astype(np.float32), p.mean_.astype(np.float32), p.singular_values_


def project(Y, loadings, mean):
    return ((Y - mean[None, :]) @ loadings.T).astype(np.float32)


def pc_reliability(scores_a, scores_b):
    """Correlation of PC scores between two complementary halves (reliability of a half-depth score)."""
    a = scores_a - scores_a.mean(0); b = scores_b - scores_b.mean(0)
    num = (a * b).sum(0)
    den = np.sqrt((a ** 2).sum(0) * (b ** 2).sum(0))
    return np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)


# ------------------------------------------------------------------ design and residualisation
def design_matrix(pcs, rel, cluster_labels, log_depth, log_area, ratio, seg_code, rel_min=None, n_knots=None):
    """[intercept | kept PCs | cluster one-hot (drop first) | depth spline | log area | ratio | seg one-hot]."""
    rel_min = C.RESID_PC_REL_MIN if rel_min is None else rel_min
    n_knots = n_knots or C.RESID_SPLINE_KNOTS
    keep = np.where(rel >= rel_min)[0]
    cols, names, is_pc, pc_rel = [np.ones((len(log_depth), 1))], ["intercept"], [False], [np.nan]
    P = pcs[:, keep] - pcs[:, keep].mean(0)
    cols.append(P); names += [f"pc{k + 1}" for k in keep]; is_pc += [True] * len(keep); pc_rel += list(rel[keep])
    L = int(cluster_labels.max()) + 1
    for l in range(1, L):
        cols.append((cluster_labels == l).astype(float)[:, None]); names.append(f"cl{l}"); is_pc.append(False); pc_rel.append(np.nan)
    sp = SplineTransformer(n_knots=n_knots, degree=3, include_bias=False).fit_transform(log_depth[:, None])
    cols.append(sp); names += [f"depth_s{j}" for j in range(sp.shape[1])]; is_pc += [False] * sp.shape[1]; pc_rel += [np.nan] * sp.shape[1]
    for v, nm in ((log_area, "log_area"), (ratio, "ratio")):
        v = np.where(np.isfinite(v), v, np.nanmedian(v))
        cols.append((v - v.mean())[:, None]); names.append(nm); is_pc.append(False); pc_rel.append(np.nan)
    S = int(seg_code.max()) + 1
    for s in range(1, S):
        cols.append((seg_code == s).astype(float)[:, None]); names.append(f"seg{s}"); is_pc.append(False); pc_rel.append(np.nan)
    X = np.hstack(cols).astype(np.float64)
    info = pd.DataFrame({"column": names, "is_pc": is_pc, "pc_rel": pc_rel})
    return X, info


def residualize(Y, X, info, eiv=True, cap=0.5, rel_res_min=0.2):
    """R = Y - X B with an errors-in-variables correction on the nuclear PC columns.

    Two stages: (1) OLS of Y and of the PC columns on the non-PC covariates (intercept, cluster one-hot, depth
    spline, morphology, segmentation); (2) regression of the Y residual on the PC residuals. The measurement-noise
    variance of PC k is (1 - rel_k) var(PC_k) (independent of the covariates), so the reliability of the
    residualised PC is rel_res_k = 1 - noise_k / var(resid PC_k). PCs with rel_res_k < rel_res_min are dropped (they
    are noise after partialling); the correction subtracts min(noise_k, cap x var(resid PC_k)) so that the corrected
    variance keeps at least (1 - cap) of the residual variance and coefficients cannot explode. info gains
    'rel_res', 'eiv_shrink' (fraction of the noise variance subtracted) and 'used'.
    """
    Y = Y.astype(np.float64)
    pc = info["is_pc"].values
    Z, P = X[:, ~pc], X[:, pc]
    ZtZ = Z.T @ Z + 1e-10 * np.trace(Z.T @ Z) / Z.shape[1] * np.eye(Z.shape[1])
    Bz_y = np.linalg.solve(ZtZ, Z.T @ Y)
    Ry = Y - Z @ Bz_y
    info["rel_res"] = np.nan; info["eiv_shrink"] = np.nan; info["used"] = ~pc
    if P.shape[1] == 0:
        info["eiv_applied"] = False
        B = np.zeros((X.shape[1], Y.shape[1])); B[~pc] = Bz_y
        return Ry.astype(np.float32), B.astype(np.float32)
    Bz_p = np.linalg.solve(ZtZ, Z.T @ P)
    Rp = P - Z @ Bz_p
    n = X.shape[0]
    rel = info.loc[pc, "pc_rel"].values
    noise = (1.0 - rel) * P.var(0)
    var_res = Rp.var(0)
    rel_res = 1.0 - noise / np.maximum(var_res, 1e-12)
    keep = rel_res >= rel_res_min
    info.loc[pc, "rel_res"] = rel_res
    pc_idx = np.where(pc)[0]
    info.loc[pc_idx[keep], "used"] = True
    Rp = Rp[:, keep]; Bz_p = Bz_p[:, keep]
    PtP = Rp.T @ Rp
    shrink = np.zeros(keep.sum())
    if eiv and keep.any():
        want = noise[keep] * n
        allowed = cap * np.diag(PtP)
        sub = np.minimum(want, allowed)
        shrink = sub / np.maximum(want, 1e-12)
        A = PtP - np.diag(sub)
        for _ in range(20):
            if np.linalg.eigvalsh(A).min() > 1e-6 * np.trace(A) / A.shape[0]:
                break
            sub = sub * 0.8; shrink = shrink * 0.8
            A = PtP - np.diag(sub)
        info["eiv_applied"] = True
    else:
        A = PtP
        info["eiv_applied"] = False
    info.loc[pc_idx[keep], "eiv_shrink"] = shrink
    if keep.any():
        Bp = np.linalg.solve(A + 1e-10 * np.trace(A) / A.shape[0] * np.eye(A.shape[0]), Rp.T @ Ry)
        R = Ry - Rp @ Bp
    else:
        Bp = np.zeros((0, Y.shape[1])); R = Ry
    B = np.zeros((X.shape[1], Y.shape[1]))
    B[pc_idx[keep]] = Bp
    B[~pc] = Bz_y - Bz_p @ Bp
    return R.astype(np.float32), B.astype(np.float32)


def n_dims_above_floor(sv, sv_floor):
    return int((sv > sv_floor.max()).sum())


def canonical_correlations(A, B):
    qa, _ = np.linalg.qr(A - A.mean(0)); qb, _ = np.linalg.qr(B - B.mean(0))
    return np.linalg.svd(qa.T @ qb, compute_uv=False)


def orient(scores, loadings, gene_names):
    """Flip each PC so its top-|loading| gene has a positive loading; returns scores, loadings, top gene names."""
    top = np.argmax(np.abs(loadings), axis=1)
    sign = np.sign(loadings[np.arange(len(top)), top]); sign[sign == 0] = 1
    return scores * sign[None, :], loadings * sign[:, None], [str(gene_names[t]) for t in top]


# ------------------------------------------------------------------ strata and stratified nulls
def strata_labels(nuc_labels, d, seg_code, n_bins=None, min_size=None):
    """nuclear cluster x depth quintile x segmentation method; strata below min_size get -1."""
    import neighborhood as NB
    n_bins = n_bins or C.N_BINS
    min_size = min_size or C.RESID_STRATA_MIN
    s = AS.cross_labels(AS.cross_labels(np.asarray(nuc_labels), NB.quantile_bins(np.log(np.maximum(d, 1)), n_bins)), np.asarray(seg_code))
    cnt = np.bincount(s[s >= 0])
    small = np.where(cnt < min_size)[0]
    s = s.copy(); s[np.isin(s, small)] = -1
    u, inv = np.unique(s[s >= 0], return_inverse=True)
    out = np.full(len(s), -1, dtype=np.int64); out[s >= 0] = inv
    return out


def stratified_permute(values, strata, rng):
    out = values.copy()
    for st in np.unique(strata[strata >= 0]):
        idx = np.where(strata == st)[0]
        out[idx] = values[rng.permutation(idx)]
    return out


def stratified_purity(labels, nn_ind, strata, n_perm=None, seed=None):
    """Per-cluster neighbour purity against permutation of labels within strata (cells with stratum -1 kept fixed)."""
    n_perm = n_perm or C.N_PERM_STRAT
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    L = labels.max() + 1

    def purity(lab):
        same = (lab[nn_ind] == lab[:, None]).mean(axis=1)
        return np.bincount(lab, weights=same, minlength=L) / np.maximum(np.bincount(lab, minlength=L), 1)

    obs = purity(labels)
    perms = np.stack([purity(stratified_permute(labels, strata, rng)) for _ in range(n_perm)])
    sizes = np.bincount(labels, minlength=L)
    return pd.DataFrame({"cluster": np.arange(L), "size": sizes, "purity": obs, "purity_perm": perms.mean(0),
                         "purity_excess": obs - perms.mean(0), "purity_p": (perms >= obs[None, :]).mean(0)})


def stratified_effects(labels, features, strata, tile_ids, n_boot=None, n_perm=None, seed=None):
    """Strata-weighted standardised mean difference (cluster vs rest within stratum) per feature, tile-bootstrap
    95% CI and stratified-permutation p. Returns long DataFrame."""
    n_boot = n_boot or C.N_BOOT
    n_perm = n_perm or C.N_PERM_STRAT
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels); strata = np.asarray(strata)
    ok = strata >= 0
    L, S = labels.max() + 1, strata[ok].max() + 1
    F = features.values.astype(float)
    F = np.where(np.isfinite(F), F, np.nanmedian(F, axis=0))
    # within-stratum SD per feature
    sd = np.zeros(F.shape[1])
    for j in range(F.shape[1]):
        m = np.bincount(strata[ok], weights=F[ok, j], minlength=S) / np.maximum(np.bincount(strata[ok], minlength=S), 1)
        sd[j] = np.std(F[ok, j] - m[strata[ok]]) + 1e-12

    def d_stat(lab, w):
        joint = strata[ok] * L + lab[ok]
        N = np.bincount(joint, weights=w[ok], minlength=S * L).reshape(S, L)
        Ns = N.sum(1, keepdims=True)
        out = np.zeros((L, F.shape[1]))
        for j in range(F.shape[1]):
            S1 = np.bincount(joint, weights=w[ok] * F[ok, j], minlength=S * L).reshape(S, L)
            S1s = S1.sum(1, keepdims=True)
            mean_in = S1 / np.maximum(N, 1e-12)
            mean_out = (S1s - S1) / np.maximum(Ns - N, 1e-12)
            valid = (N > 0) & ((Ns - N) > 0)
            diff = np.where(valid, mean_in - mean_out, 0.0)
            wgt = np.where(valid, N, 0.0)
            out[:, j] = (wgt * diff).sum(0) / np.maximum(wgt.sum(0), 1e-12) / sd[j]
        return out

    ones = np.ones(len(labels))
    obs = d_stat(labels, ones)
    lo, hi = CP.bootstrap_tiles(lambda w: d_stat(labels, w).ravel(), np.asarray(tile_ids), n_boot=n_boot, seed=seed)
    perms = np.stack([d_stat(stratified_permute(labels, strata, rng), ones) for _ in range(n_perm)])
    p = (np.abs(perms) >= np.abs(obs)[None]).mean(0)
    rows = []
    for c in range(L):
        for j, f in enumerate(features.columns):
            rows.append({"cluster": c, "feature": f, "d": obs[c, j], "d_lo": lo.reshape(L, -1)[c, j],
                         "d_hi": hi.reshape(L, -1)[c, j], "p_perm": p[c, j]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ continuous score tests
def morans_i(z, nn_ind):
    z = np.asarray(z, dtype=float); zc = z - z.mean()
    lag = zc[nn_ind].mean(axis=1)
    den = (zc ** 2).sum()
    return float((zc * lag).sum() / den) if den > 0 else np.nan


def morans_test(z, nn_ind, strata, n_perm=None, seed=None):
    n_perm = n_perm or C.N_PERM_STRAT
    seed = C.SEED if seed is None else seed
    rng = np.random.default_rng(seed)
    z = np.asarray(z, dtype=float)
    # centre within strata so the statistic reflects within-subtype coherence
    zc = z.copy()
    for st in np.unique(strata[strata >= 0]):
        idx = strata == st
        zc[idx] = z[idx] - z[idx].mean()
    zc[strata < 0] = 0.0
    obs = morans_i(zc, nn_ind)
    perms = np.array([morans_i(stratified_permute(zc, strata, rng), nn_ind) for _ in range(n_perm)])
    return {"I": obs, "I_perm_mean": float(perms.mean()), "I_perm_sd": float(perms.std()),
            "I_perm_p99": float(np.percentile(perms, 99)), "p_perm": float((perms >= obs).mean())}


def dose_response(score, feature, strata, tile_ids, n_bins=None, n_boot=None, seed=None):
    """Mean stratum-centred score (SD units) per feature quintile; effect = top - bottom quintile with tile CI."""
    import neighborhood as NB
    n_bins = n_bins or C.N_BINS
    n_boot = n_boot or C.N_BOOT
    score = np.asarray(score, dtype=float); strata = np.asarray(strata)
    ok = strata >= 0
    zc = score.copy()
    for st in np.unique(strata[ok]):
        idx = strata == st
        zc[idx] = score[idx] - score[idx].mean()
    zc = zc / (zc[ok].std() + 1e-12)
    b = NB.quantile_bins(np.asarray(feature, dtype=float), n_bins)
    use = ok & (b >= 0)

    def means(w):
        num = np.bincount(b[use], weights=w[use] * zc[use], minlength=n_bins)
        den = np.bincount(b[use], weights=w[use], minlength=n_bins)
        return num / np.maximum(den, 1e-12)

    m = means(np.ones(len(zc)))
    lo, hi = CP.bootstrap_tiles(lambda w: np.append(means(w), means(w)[-1] - means(w)[0]), np.asarray(tile_ids), n_boot=n_boot, seed=seed)
    return {"bin_means": m, "bin_lo": lo[:n_bins], "bin_hi": hi[:n_bins],
            "effect": float(m[-1] - m[0]), "effect_lo": float(lo[-1]), "effect_hi": float(hi[-1])}


# ------------------------------------------------------------------ eligibility and the panel rule
def eligibility(summary, n_cells, sample_seg_share, floor_purity_max):
    """Apply the pre-registered cluster rule; returns the summary with flags, rank and 'eligible'.

    The interface criterion may be met by either of the two interface features (tumor fraction = boundary vs core,
    immune kernel density = immune-rich); the better one ranks the cluster and gives its label.
    """
    s = summary.copy()
    s["ok_size"] = s["size"] >= max(C.RESID_MIN_CLUSTER, C.RESID_MIN_CLUSTER_FRAC * n_cells)
    s["ok_repro"] = s["repro_f1"] >= C.RESID_REPRO_F1
    s["ok_depth"] = s["depth_ratio"] >= C.RESID_DEPTH_RATIO_MIN
    s["ok_seg"] = ~((s["seg_share"] > C.RESID_SEG_SHARE_MAX) & (sample_seg_share < 0.8))
    s["morphology_driven"] = (s["d_log_cell_area"].abs() >= C.RESID_MORPH_SD_MAX) | (s["d_nuc_cell_ratio"].abs() >= C.RESID_MORPH_SD_MAX)
    s["ok_purity"] = (s["purity_excess"] >= C.RESID_PURITY_EXCESS) & (s["purity_excess"] > floor_purity_max)
    best_stat = np.zeros(len(s)); best_feat = np.array([""] * len(s), dtype=object); ok_eff = np.zeros(len(s), dtype=bool)
    for pf in C.RESID_INTERFACE_FEATURES:
        ci_excl = (s[f"d_{pf}_lo"] > 0) | (s[f"d_{pf}_hi"] < 0)
        same_sign = np.sign(s[f"d_{pf}"]) == np.sign(s[f"d_{pf}_h2"])
        ok = (s[f"d_{pf}"].abs() >= C.RESID_EFFECT_SD) & ci_excl & same_sign
        lb = np.minimum(np.abs(s[f"d_{pf}_lo"]), np.abs(s[f"d_{pf}_hi"]))
        stat = np.where(ok, np.minimum(lb, s[f"d_{pf}_h2"].abs()), 0.0)
        better = stat > best_stat
        best_stat = np.where(better, stat, best_stat); best_feat = np.where(better, pf, best_feat)
        ok_eff |= ok.values
    s["ok_effect"] = ok_eff
    s["effect_feature"] = best_feat
    s["rank_stat"] = best_stat
    s["ok_de"] = s["n_cyto_only_de"] >= C.RESID_MIN_CYTO_ONLY_DE
    s["eligible"] = s[["ok_size", "ok_repro", "ok_depth", "ok_seg", "ok_purity", "ok_effect", "ok_de"]].all(axis=1) & ~s["morphology_driven"]
    s["label"] = np.where(s["effect_feature"] == "immune_kernel", np.where(s["d_immune_kernel"] > 0, "immune-rich", "immune-poor"),
                          np.where(s["d_tumor_frac"] < 0, "boundary", "core"))
    s["rank"] = np.nan
    el = s[s["eligible"]].sort_values("rank_stat", ascending=False)
    s.loc[el.index, "rank"] = np.arange(1, len(el) + 1)
    return s


def zoom_window(xy, mask, size_um=None):
    """The size x size window (um) holding the most masked cells; grid search on a size/4 lattice."""
    size_um = size_um or C.RESID_ZOOM_UM
    pts = xy[mask]
    if len(pts) == 0:
        return None
    step = size_um / 4
    x0, y0 = np.floor(xy[:, 0].min()), np.floor(xy[:, 1].min())
    ix = ((pts[:, 0] - x0) // step).astype(int); iy = ((pts[:, 1] - y0) // step).astype(int)
    grid = np.zeros((ix.max() + 5, iy.max() + 5))
    np.add.at(grid, (ix, iy), 1)
    # sum over 4 x 4 lattice blocks
    cs = grid.cumsum(0).cumsum(1)
    cs = np.pad(cs, ((1, 0), (1, 0)))
    k = 4
    tot = cs[k:, k:] - cs[:-k, k:] - cs[k:, :-k] + cs[:-k, :-k]
    i, j = np.unravel_index(np.argmax(tot), tot.shape)
    return (float(x0 + i * step), float(y0 + j * step), float(size_um))
