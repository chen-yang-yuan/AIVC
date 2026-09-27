"""Step functions of the Fig. 1 exploration. The notebook and 1_heavy.py both call these.

Light steps (notebook, all cells): prepare, step0_qc, step1_compartment, step2_concordance, step5_controls,
step8_pathways_light. Heavy steps (HGCC, all kept cells; locally on a subsample in QUICK mode): nbhd, step3_clustering,
step4_de, step7_association, step8_pathways_assoc. Every step caches its tables under output/2_exploration/<ds>/
(or <ds>_quick/ in QUICK mode) and is skipped when the files exist unless force=True.
"""

import json
import time

import numpy as np
import pandas as pd
from scipy import sparse

import config as C
import io_utils as IO
import compartment as CP
import clustering as CL
import de as DE
import neighborhood as NB
import association as AS
import controls as CT

AXES_MAIN = ["subtype_nuc", "niche", "tumor_frac_bin"]


def _t(t0):
    return f"{time.time() - t0:6.1f}s"


# ====================================================================== preparation
def prepare(ds, quick=False, n_cells=None, seed=None, verbose=True):
    """Load one dataset and build every shared object (depths, masks, rarefied matrices, tiles, offsets)."""
    t0 = time.time()
    seed = C.SEED if seed is None else seed
    n_cells = n_cells or C.QUICK_N_CELLS
    ctx = {"ds": ds, "quick": quick, "seed": seed}
    ctx["out"] = IO.out_dir(ds + ("_quick" if quick else ""))
    ctx["fig"] = IO.out_dir(ds + ("_quick" if quick else ""), "fig")
    ctx["n_boot"] = C.QUICK_N_BOOT if quick else C.N_BOOT
    ctx["B_null"] = C.QUICK_B_NULL if quick else C.B_NULL
    ctx["n_sets"] = C.QUICK_N_MATCHED_SETS if quick else C.N_MATCHED_SETS

    panel = IO.load_panel()
    obs, cids = IO.load_tumor_obs(ds)
    Xn, Xc, genes = IO.load_compartments(ds, panel)
    dep = CP.cell_depths(Xn, Xc)
    has_nuc = obs["nucleus_area"].notna().values
    keep = has_nuc & (dep["d_match"].values >= C.D_MIN)
    if C.REQUIRE_ONE_NUCLEUS and "nucleus_count" in obs:
        keep &= obs["nucleus_count"].values == 1
    sel = np.where(keep)[0]
    if quick and len(sel) > n_cells:
        sel = np.sort(np.random.default_rng(seed).choice(sel, n_cells, replace=False))
    in_heavy = np.zeros(len(obs), dtype=bool); in_heavy[sel] = True

    cell_area = obs["cell_area"].values.astype(float)
    nuc_area = obs["nucleus_area"].values.astype(float)
    ratio = np.where(has_nuc, nuc_area / np.maximum(cell_area, 1e-6), np.nan)
    cells = pd.DataFrame({"cell_id": cids, "global_x": obs["global_x"].values, "global_y": obs["global_y"].values,
                          "nuc_depth": dep["nuc_depth"], "cyto_depth": dep["cyto_depth"],
                          "in_cell_depth": dep["in_cell_depth"], "nuc_frac": dep["nuc_frac"], "d_match": dep["d_match"],
                          "has_nucleus": has_nuc, "nucleus_count": obs.get("nucleus_count", pd.Series(1, index=obs.index)).values,
                          "seg_method": obs["segmentation_method"].values if "segmentation_method" in obs else "unknown",
                          "cell_area": cell_area, "nucleus_area": nuc_area, "nuc_cell_ratio": ratio,
                          "keep": keep, "in_heavy": in_heavy})
    cells["ratio_bin"] = NB.quantile_bins(ratio)
    cells["area_bin"] = NB.quantile_bins(cell_area)
    cells["seg_code"] = pd.Categorical(cells["seg_method"]).codes.astype(np.int64)

    # matched-depth matrices on the selected cells
    d = dep["d_match"].values[sel].astype(np.int64)
    nm, cm, tm, tf = CP.rarefy_matched(Xn[sel], Xc[sel], d, seed)
    mats = {"nuc_m": nm, "cyto_m": cm, "total_m": tm, "total_full": tf}
    mats["nuc_h1"], mats["nuc_h2"] = CP.thin_half(nm, seed + 11)
    mats["cyto_h1"], mats["cyto_h2"] = CP.thin_half(cm, seed + 13)
    xy = cells[["global_x", "global_y"]].values[sel]
    tile_um = CP.choose_tile_um(xy)
    tile_id = CP.spatial_tiles(xy, tile_um)
    cells["tile_id"] = -1
    cells.loc[sel, "tile_id"] = tile_id
    cells["tile_um"] = tile_um

    pair = CP.Pair(Xn, Xc)
    ctx.update(panel=panel, obs=obs, cids=cids, Xn=Xn, Xc=Xc, genes=genes, cells=cells, sel=sel, d=d, mats=mats,
               tile_id=tile_id, tile_um=tile_um, pair=pair, offset=pair.loo_offset(),
               mean_total=pair.gene_mean_total(), xy_sel=xy)
    ctx["gene_mask_cluster"] = CL.shared_gene_mask(nm, cm)
    if verbose:
        print(f"  [{C.short(ds)}] prepared: {len(obs)} tumor cells, keep {keep.mean():.3f}, "
              f"heavy set {len(sel)}, tile {tile_um} um ({tile_id.max() + 1} tiles), "
              f"cluster genes {ctx['gene_mask_cluster'].sum()}  {_t(t0)}")
    return ctx


# ====================================================================== light steps
def step0_qc(ctx, force=False):
    out = ctx["out"]
    cells = ctx["cells"]
    genes = list(ctx["genes"])
    ctrl_cyto = [g for g in C.POS_CTRL_CYTO + C.POS_CTRL_CYTO_CANDIDATES if g in genes]
    top = list(np.asarray(genes)[np.argsort(-ctx["mean_total"])[:C.N_TOP_ABUNDANT_CTRL]])
    qc = {"dataset": ctx["ds"], "cancer_type": C.CANCER_TYPE.get(ctx["ds"]), "quick": ctx["quick"],
          "n_tumor_cells": int(len(cells)), "n_keep": int(cells["keep"].sum()), "frac_keep": float(cells["keep"].mean()),
          "n_heavy": int(cells["in_heavy"].sum()), "frac_no_nucleus": float(1 - cells["has_nucleus"].mean()),
          "median_nuc_depth": float(cells["nuc_depth"].median()), "median_cyto_depth": float(cells["cyto_depth"].median()),
          "median_d_match": float(cells["d_match"].median()),
          "nuclear_share_of_reads": float(cells["nuc_depth"].sum() / cells["in_cell_depth"].sum()),
          "median_nuc_frac_per_cell": float(cells["nuc_frac"].median()),
          "seg_method_counts": cells["seg_method"].value_counts().to_dict(),
          "tile_um": int(ctx["tile_um"]), "n_tiles": int(ctx["tile_id"].max() + 1),
          "genes_mean_ge_0.5": int((ctx["mean_total"] >= 0.5).sum()), "genes_mean_ge_0.2": int((ctx["mean_total"] >= 0.2).sum()),
          "genes_mean_ge_0.05": int((ctx["mean_total"] >= 0.05).sum()),
          "pos_ctrl_nuclear_present": [g for g in C.POS_CTRL_NUCLEAR if g in genes],
          "pos_ctrl_cyto_present": ctrl_cyto, "top_abundant_genes": top}
    IO.save(out + "qc_summary.json", qc)
    IO.save(out + "cells.parquet", cells)
    ctx["qc"] = qc
    return qc


def step1_compartment(ctx, force=False):
    """Per-gene localisation (log-OR), reliability, classes and strata checks. All tumor cells."""
    out, t0 = ctx["out"], time.time()

    def compute():
        pair, off, genes = ctx["pair"], ctx["offset"], ctx["genes"]
        lo = CP.per_gene_logor(pair, off, genes)
        rel = CP.per_gene_reliability(pair, off, lo["beta"].values, B=ctx["B_null"], seed=ctx["seed"])
        df = pd.concat([lo, rel], axis=1)
        df["class"] = CP.classify_genes(df).values
        # strata checks: localisation per segmentation method and per area-ratio bin
        cells = ctx["cells"]
        for name, lab in (("seg", cells["seg_code"].values), ("ratio", cells["ratio_bin"].values)):
            beta, _, _, N = AS.binom_axis(pair, off, lab)
            for l in range(beta.shape[0]):
                df[f"beta_{name}{l}"] = beta[l]
            b = np.where(N >= 50, beta, np.nan)
            df[f"beta_{name}_range"] = np.nanmax(b, 0) - np.nanmin(b, 0)
        df["is_pos_ctrl_nuc"] = df["gene"].isin(C.POS_CTRL_NUCLEAR)
        df["is_pos_ctrl_cyto"] = df["gene"].isin(ctx["qc"]["pos_ctrl_cyto_present"] + ctx["qc"]["top_abundant_genes"])
        return df

    df = IO.cached(out + "genes_compartment.parquet", compute, force)
    ctx["genes_compartment"] = df
    print(f"    step1 compartment: {df['class'].value_counts().to_dict()}  {_t(t0)}")
    return df


def step2_concordance(ctx, force=False):
    """Split-half concordance at cell level (selected cells, matched depth) and tile level."""
    out, t0 = ctx["out"], time.time()
    genes, mats = ctx["genes"], ctx["mats"]
    mean_sel = np.asarray((ctx["Xn"][ctx["sel"]] + ctx["Xc"][ctx["sel"]]).mean(0)).ravel()

    def cell_level():
        return CP.concordance(mats["nuc_m"], mats["cyto_m"], genes, seed=ctx["seed"], gene_mask=mean_sel >= C.MEAN_MIN_CELL)

    def tile_level():
        T = CP.tile_indicator(ctx["tile_id"])
        return CP.concordance(CP.pool(mats["nuc_m"], T), CP.pool(mats["cyto_m"], T), genes, seed=ctx["seed"],
                              gene_mask=mean_sel >= C.MEAN_MIN_TILE)

    cc = IO.cached(out + "genes_concordance_cell.parquet", cell_level, force)
    ct = IO.cached(out + "genes_concordance_tile.parquet", tile_level, force)
    ctx["conc_cell"], ctx["conc_tile"] = cc, ct
    print(f"    step2 concordance: cell {cc['r_true'].notna().sum()} genes, tile {ct['r_true'].notna().sum()} genes; "
          f"median r_true {cc['r_true'].median():.2f} / {ct['r_true'].median():.2f}  {_t(t0)}")
    return cc, ct


def step5_controls(ctx, force=False):
    """Abundance-bin percentiles for the per-gene statistics and the SG-set vs matched-control comparison."""
    out, t0 = ctx["out"], time.time()
    gc = CT.add_bin_percentiles(ctx["genes_compartment"], ["beta", "rel_param", "beta_seg_range", "beta_ratio_range", "phi"])
    IO.save(out + "genes_compartment.parquet", gc)
    ctx["genes_compartment"] = gc
    for key, fn in (("conc_cell", "genes_concordance_cell.parquet"), ("conc_tile", "genes_concordance_tile.parquet")):
        df = ctx[key]
        df = CT.add_bin_percentiles(df, ["r_true", "rel_n", "rel_c"], mean_col="mean_per_unit")
        df["r_true_low_pct"] = 100 - df["r_true_pct"]
        IO.save(out + fn, df)
        ctx[key] = df

    def sets():
        panel, genes = list(ctx["panel"]), ctx["genes"]
        sg = IO.load_sg_genes(panel)
        ctrl = CT.matched_control_genes(sg, panel, ctx["mean_total"], seed=ctx["seed"])
        pos = {g: i for i, g in enumerate(genes)}
        rows = []
        for name, gl in (("SG_markers", sg), ("matched_control", ctrl)):
            idx = np.array([pos[g] for g in gl])
            M = sparse.csr_matrix((np.ones(len(idx)), (idx, np.zeros(len(idx), dtype=int))), shape=(len(genes), 1))
            pr = CP.pair_from_sets(ctx["Xn"], ctx["Xc"], M, ctx["pair"].M, ctx["pair"].Q)
            off = pr.loo_offset()
            lo = CP.per_gene_logor(pr, off, [name])
            rel = CP.per_gene_reliability(pr, off, lo["beta"].values, B=ctx["B_null"], seed=ctx["seed"])
            r = pd.concat([lo, rel], axis=1).iloc[0].to_dict()
            r["set"] = name; r["n_genes"] = len(gl)
            r["log_mean_expr"] = float(np.log1p(ctx["mean_total"][idx]).mean())
            rows.append(r)
        return pd.DataFrame(rows)

    ctx["control_sets"] = IO.cached(out + "controls_sets.parquet", sets, force)
    print(f"    step5 controls: percentiles added; SG vs matched control beta "
          f"{ctx['control_sets']['beta'].round(3).tolist()}, reliability {ctx['control_sets']['rel_param'].round(3).tolist()}  {_t(t0)}")
    return ctx["control_sets"]


def _pathway_objects(ctx):
    if "pw_M" not in ctx:
        sg = IO.load_sg_genes(list(ctx["panel"]))
        M, info = IO.load_pathways(ctx["genes"], extra_sets={"SG_markers": ("SG_markers.xlsx", sg)})
        ctx["pw_M"], ctx["pw_info"] = M, info
    return ctx["pw_M"], ctx["pw_info"]


def step8_pathways_light(ctx, force=False):
    """Pathway-level localisation, reliability and tile concordance, with matched random-set nulls."""
    out, t0 = ctx["out"], time.time()
    M, info = _pathway_objects(ctx)
    pair, Xn, Xc = ctx["pair"], ctx["Xn"], ctx["Xc"]

    def compute():
        pr = CP.pair_from_sets(Xn, Xc, M, pair.M, pair.Q)
        off = pr.loo_offset()
        lo = CP.per_gene_logor(pr, off, info["pathway"].values)
        rel = CP.per_gene_reliability(pr, off, lo["beta"].values, B=ctx["B_null"], seed=ctx["seed"])
        df = pd.concat([info.reset_index(drop=True), lo.drop(columns=["gene"]), rel], axis=1)
        # matched random-set null per pathway
        bins = CT.expression_bins(ctx["mean_total"])
        Mc = sparse.csc_matrix(M)
        null_beta, null_rel, pct_beta, pct_rel = [], [], [], []
        for j in range(M.shape[1]):
            members = Mc.indices[Mc.indptr[j]:Mc.indptr[j + 1]]
            sets = CT.matched_random_sets(members, bins, n_sets=ctx["n_sets"], seed=ctx["seed"] + j)
            rows = np.concatenate(sets); cols = np.repeat(np.arange(len(sets)), [len(s) for s in sets])
            Mn = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(M.shape[0], len(sets)))
            pn = CP.pair_from_sets(Xn, Xc, Mn, pair.M, pair.Q)
            on = pn.loo_offset()
            ln = CP.per_gene_logor(pn, on, np.arange(len(sets)))
            rn = CP.per_gene_reliability(pn, on, ln["beta"].values, B=1, seed=ctx["seed"])
            null_beta.append(np.nanmean(np.abs(ln["beta"]))); null_rel.append(np.nanmean(rn["rel_param"]))
            pct_beta.append(CT.set_null_percentile(abs(df.loc[j, "beta"]), np.abs(ln["beta"].values)))
            pct_rel.append(CT.set_null_percentile(df.loc[j, "rel_param"], rn["rel_param"].values))
        df["null_abs_beta_mean"] = null_beta; df["null_rel_mean"] = null_rel
        df["abs_beta_null_pct"] = pct_beta; df["rel_null_pct"] = pct_rel
        return df

    def tile_conc():
        T = CP.tile_indicator(ctx["tile_id"])
        an = CP.aggregate(ctx["mats"]["nuc_m"], M); ac = CP.aggregate(ctx["mats"]["cyto_m"], M)
        return CP.concordance(CP.pool(an, T), CP.pool(ac, T), info["pathway"].values, seed=ctx["seed"])

    ctx["pathways_compartment"] = IO.cached(out + "pathways_compartment.parquet", compute, force)
    ctx["pathways_concordance"] = IO.cached(out + "pathways_concordance.parquet", tile_conc, force)
    print(f"    step8 pathways (light): {len(info)} sets; reliability > null 95th pct in "
          f"{int((ctx['pathways_compartment']['rel_null_pct'] >= 95).sum())}  {_t(t0)}")
    return ctx["pathways_compartment"], ctx["pathways_concordance"]


# ====================================================================== heavy steps
def nbhd(ctx, force=False):
    """Neighbourhood features, niches and bins for the selected cells."""
    out, t0 = ctx["out"], time.time()

    def compute():
        all_obs = IO.load_all_obs(ctx["ds"])
        feat = NB.nbhd_features(all_obs, ctx["cids"][ctx["sel"]])
        comp = feat[[c for c in feat.columns if c.startswith("comp_")]]
        feat["niche"] = NB.niche_labels(comp, seed=ctx["seed"])
        feat["crowding_bin"] = NB.quantile_bins(feat["n_tumor_50"].values)
        feat["tumor_frac_bin"] = NB.quantile_bins(feat["tumor_frac"].values)
        feat["immune_frac_bin"] = NB.quantile_bins(feat["immune_kernel"].values)
        feat["dist_nontumor_bin"] = NB.quantile_bins(feat["dist_nontumor"].values)
        feat.insert(0, "cell_id", ctx["cids"][ctx["sel"]])
        return feat

    ctx["nbhd"] = IO.cached(out + "nbhd_features.parquet", compute, force)
    print(f"    nbhd features: {ctx['nbhd'].shape[1]} columns  {_t(t0)}")
    return ctx["nbhd"]


def step3_clustering(ctx, force=False, verbose=True):
    out, t0 = ctx["out"], time.time()
    mats, gm = ctx["mats"], ctx["gene_mask_cluster"]

    def compute():
        labels, pcs = {}, {}
        for name in C.CLUSTER_MATRICES:
            t1 = time.time()
            lab, p = CL.cluster_matrix(mats[name], gm, seed=ctx["seed"])
            for r, l in lab.items():
                labels[f"{name}_r{r}"] = l
            pcs[name] = p
            if verbose:
                print(f"      {name:10s} " + ", ".join(f"res {r}: {l.max() + 1}" for r, l in lab.items()) + f"  {_t(t1)}")
        df = pd.DataFrame(labels)
        df.insert(0, "cell_id", ctx["cids"][ctx["sel"]])
        IO.save(out + "pcs.npz", pcs)
        return df

    cl = IO.cached(out + "clusters.parquet", compute, force)
    ctx["clusters"] = cl
    r = C.PRIMARY_RES
    lab = {m: cl[f"{m}_r{r}"].values for m in C.CLUSTER_MATRICES}

    def summary():
        pw = CL.all_pairwise(lab)
        for res in C.RESOLUTIONS:
            if res == r:
                continue
            lab2 = {m: cl[f"{m}_r{res}"].values for m in C.CLUSTER_MATRICES}
            p2 = CL.all_pairwise(lab2); p2["resolution"] = res
            pw["resolution"] = r
            pw = pd.concat([pw, p2], ignore_index=True)
        pw["resolution"] = pw["resolution"].fillna(r)
        return pw

    ctx["cluster_pairwise"] = IO.cached(out + "clustering_pairwise.parquet", summary, force)
    ctx["cyto_only"] = IO.cached(out + "cyto_only_clusters.parquet", lambda: CL.cyto_only_clusters(lab), force)
    ctx["recovery"] = IO.cached(out + "cluster_recovery_total_by_nuc.parquet",
                                lambda: CL.cluster_recovery(lab["total_full"], lab["nuc_m"]), force)
    IO.save(out + "jaccard_cyto_vs_nuc.csv", CL.jaccard_matrix(lab["cyto_m"], lab["nuc_m"]).reset_index())
    IO.save(out + "jaccard_nuc_vs_total_full.csv", CL.jaccard_matrix(lab["nuc_m"], lab["total_full"]).reset_index())

    def spatial():
        nn = CL.knn_indices(ctx["xy_sel"])
        rows = []
        for m in ("cyto_m", "nuc_m", "total_full"):
            s = CL.spatial_purity(lab[m], nn, seed=ctx["seed"]); s["matrix"] = m
            s = s.merge(CL.cluster_depth(lab[m], ctx["d"]), on=["cluster", "size"])
            rows.append(s)
        return pd.concat(rows, ignore_index=True)

    ctx["cluster_spatial"] = IO.cached(out + "cluster_spatial.parquet", spatial, force)

    def features():
        f = ctx["nbhd"]
        cols = ["tumor_frac", "immune_kernel", "stromal_kernel", "log1p_crowding", "log1p_dist_nontumor"]
        feats = pd.concat([f[cols].reset_index(drop=True),
                           ctx["cells"].loc[ctx["sel"], ["nuc_cell_ratio", "cell_area"]].reset_index(drop=True),
                           pd.DataFrame({"d_match": ctx["d"]})], axis=1)
        rows = []
        for m in ("cyto_m", "nuc_m"):
            e = CL.cluster_feature_enrichment(lab[m], feats, seed=ctx["seed"]); e["matrix"] = m
            rows.append(e)
        return pd.concat(rows, ignore_index=True)

    if "nbhd" in ctx:
        ctx["cluster_features"] = IO.cached(out + "cluster_features.parquet", features, force)
    ceil = ctx["cluster_pairwise"].query("resolution == @r").set_index(["a", "b"])["ari"]
    print(f"    step3 clustering: ARI nuc/cyto {ceil.get(('nuc_m', 'cyto_m'), np.nan):.2f}, "
          f"nuc/total_full {ceil.get(('nuc_m', 'total_full'), np.nan):.2f}, ceilings nuc-halves "
          f"{ceil.get(('nuc_h1', 'nuc_h2'), np.nan):.2f}, cyto-halves {ceil.get(('cyto_h1', 'cyto_h2'), np.nan):.2f}, "
          f"cross-halves {ceil.get(('nuc_h1', 'cyto_h1'), np.nan):.2f}; cyto-only clusters "
          f"{int(ctx['cyto_only']['cyto_only'].sum())}  {_t(t0)}")
    return cl


def step4_de(ctx, force=False):
    out, t0 = ctx["out"], time.time()
    mats, genes, cl = ctx["mats"], ctx["genes"], ctx["clusters"]
    r = C.PRIMARY_RES
    mean_sel = np.asarray((ctx["Xn"][ctx["sel"]] + ctx["Xc"][ctx["sel"]]).mean(0)).ravel()
    gmask = mean_sel >= C.MEAN_MIN_TILE
    res = {}
    for part, lab in (("nucpart", cl[f"nuc_h1_r{r}"].values), ("cytopart", cl[f"cyto_h1_r{r}"].values)):
        def paired(lab=lab):
            df = DE.paired_logfc(mats["nuc_h2"], mats["cyto_h2"], lab, genes, ctx["tile_id"], gene_mask=gmask,
                                 n_boot=ctx["n_boot"], seed=ctx["seed"])
            return DE.call_de(df)
        res[part] = IO.cached(out + f"de_paired_{part}.parquet", paired, force)

        def topk(lab=lab):
            tn = DE.wilcoxon_topk(mats["nuc_h2"], genes, lab, gene_mask=ctx["gene_mask_cluster"])
            tc = DE.wilcoxon_topk(mats["cyto_h2"], genes, lab, gene_mask=ctx["gene_mask_cluster"])
            return {"nuclear": tn, "cytoplasmic": tc}
        tk = IO.cached(out + f"de_topk_{part}.json", topk, force)
        ov = DE.topk_overlap(tk["nuclear"], tk["cytoplasmic"]); ov["partition"] = part
        res[part + "_overlap"] = ov
    ov = pd.concat([res["nucpart_overlap"], res["cytopart_overlap"]], ignore_index=True)
    IO.save(out + "de_topk_overlap.parquet", ov)
    # cluster centroids of the cytoplasmic partition for cross-sample recurrence
    cen = DE.cluster_centroids(mats["cyto_m"], cl[f"cyto_m_r{r}"].values)
    IO.save(out + "centroids_cyto_m.npz", {"centroids": cen, "genes": np.asarray(genes)})
    ctx["de_nucpart"], ctx["de_cytopart"], ctx["de_overlap"] = res["nucpart"], res["cytopart"], ov
    d = res["nucpart"]
    print(f"    step4 DE (nuclear partition): cyto-only {int(d['cyto_only'].sum())}, nuc-only {int(d['nuc_only'].sum())}, "
          f"shared {int(d['shared'].sum())}; top-{C.TOPK_DE} marker Jaccard median "
          f"{res['nucpart_overlap']['jaccard'].median():.2f}  {_t(t0)}")
    return res


def _axes(ctx):
    cl, f, cells = ctx["clusters"], ctx["nbhd"], ctx["cells"].loc[ctx["sel"]]
    r = C.PRIMARY_RES
    axes = {"subtype_nuc": cl[f"nuc_h1_r{r}"].values, "subtype_cyto": cl[f"cyto_h1_r{r}"].values,
            "subtype_total": cl[f"total_full_r{r}"].values,
            "niche": f["niche"].values, "crowding_bin": f["crowding_bin"].values,
            "tumor_frac_bin": f["tumor_frac_bin"].values, "immune_frac_bin": f["immune_frac_bin"].values,
            "dist_nontumor_bin": f["dist_nontumor_bin"].values,
            "morph_ratio_bin": cells["ratio_bin"].values, "morph_area_bin": cells["area_bin"].values,
            "seg_method": cells["seg_code"].values}
    return axes


def step7_association(ctx, force=False):
    out, t0 = ctx["out"], time.time()
    sel = ctx["sel"]
    Xn, Xc = ctx["Xn"][sel], ctx["Xc"][sel]
    pair = CP.Pair(Xn, Xc)
    off = pair.loo_offset()
    gidx = np.where(ctx["mean_total"] >= C.MEAN_MIN_ASSOC)[0]
    axes = _axes(ctx)

    def compute():
        rows = []
        for name, lab in axes.items():
            t1 = time.time()
            rows.append(AS.run_axis(pair, off, Xn, Xc, lab, gidx, ctx["genes"], name, seed=ctx["seed"]))
            print(f"      axis {name:18s} {int(lab.max() + 1):3d} levels  {_t(t1)}")
        # neighbourhood beyond subtype: cross-classification increments
        for name in ("niche", "tumor_frac_bin", "immune_frac_bin", "morph_ratio_bin"):
            t1 = time.time()
            cross = AS.cross_labels(axes["subtype_nuc"], axes[name])
            rows.append(AS.run_axis(pair, off, Xn, Xc, cross, gidx, ctx["genes"], f"subtype_nuc_x_{name}", seed=ctx["seed"]))
            print(f"      axis subtype_nuc_x_{name:14s} {int(cross.max() + 1):3d} levels  {_t(t1)}")
        long = pd.concat(rows, ignore_index=True)
        return long

    long = IO.cached(out + "assoc_genes_long.parquet", compute, force)
    wide = AS.wide_summary(long)
    for name in ("niche", "tumor_frac_bin", "immune_frac_bin", "morph_ratio_bin"):
        for m in ("cyto", "nuc", "frac"):
            a, b = f"subtype_nuc_x_{name}__{m}", f"subtype_nuc__{m}"
            if a in wide and b in wide:
                wide[f"{name}_beyond_subtype__{m}"] = wide[a] - wide[b]
    IO.save(out + "assoc_genes_wide.parquet", wide)
    ctx["assoc_long"], ctx["assoc_wide"] = long, wide
    top = wide.set_index("gene")[[c for c in wide.columns if c.endswith("__frac") and "_x_" not in c]].max(axis=1)
    print(f"    step7 association: {len(gidx)} genes x {long['axis'].nunique()} axes; top fraction-model genes "
          f"{top.sort_values(ascending=False).head(5).round(3).to_dict()}  {_t(t0)}")
    return long, wide


def step8_pathways_assoc(ctx, force=False):
    out, t0 = ctx["out"], time.time()
    M, info = _pathway_objects(ctx)
    sel = ctx["sel"]
    Xn, Xc = ctx["Xn"][sel], ctx["Xc"][sel]
    depth_M = np.asarray(Xc.sum(axis=1)).ravel().astype(float); depth_Q = np.asarray(Xn.sum(axis=1)).ravel().astype(float)
    An, Ac = CP.aggregate(Xn, M), CP.aggregate(Xc, M)
    pr = CP.Pair(An, Ac, M=depth_M, Q=depth_Q)
    off = pr.loo_offset()
    axes = _axes(ctx)
    names = info["pathway"].values
    all_idx = np.arange(len(names))

    def compute():
        rows = [AS.run_axis(pr, off, An, Ac, lab, all_idx, names, name, seed=ctx["seed"]) for name, lab in axes.items()]
        for name in ("niche", "tumor_frac_bin"):
            cross = AS.cross_labels(axes["subtype_nuc"], axes[name])
            rows.append(AS.run_axis(pr, off, An, Ac, cross, all_idx, names, f"subtype_nuc_x_{name}", seed=ctx["seed"]))
        long = pd.concat(rows, ignore_index=True)
        # matched random-set null for the main axes
        bins = CT.expression_bins(ctx["mean_total"])
        Mc = sparse.csc_matrix(M)
        null_rows = []
        for j in range(M.shape[1]):
            members = Mc.indices[Mc.indptr[j]:Mc.indptr[j + 1]]
            sets = CT.matched_random_sets(members, bins, n_sets=ctx["n_sets"], seed=ctx["seed"] + j)
            rr = np.concatenate(sets); cc = np.repeat(np.arange(len(sets)), [len(s) for s in sets])
            Mn = sparse.csr_matrix((np.ones(len(rr)), (rr, cc)), shape=(M.shape[0], len(sets)))
            Bn, Bc = CP.aggregate(Xn, Mn), CP.aggregate(Xc, Mn)
            pn = CP.Pair(Bn, Bc, M=depth_M, Q=depth_Q); on = pn.loo_offset()
            for name in AXES_MAIN:
                ln = AS.run_axis(pn, on, Bn, Bc, axes[name], np.arange(len(sets)), np.arange(len(sets)), name, n_perm=0)
                for m in ("cyto", "nuc", "frac"):
                    obs_v = long.query("axis == @name and model == @m and gene == @names[@j]")["dev_expl"]
                    nullv = ln.query("model == @m")["dev_expl"].values
                    null_rows.append({"pathway": names[j], "axis": name, "model": m,
                                      "dev_expl_null_mean": float(np.nanmean(nullv)),
                                      "dev_expl_null_pct": CT.set_null_percentile(float(obs_v.iloc[0]) if len(obs_v) else np.nan, nullv)})
        null = pd.DataFrame(null_rows)
        long = long.rename(columns={"gene": "pathway"}).merge(null, on=["pathway", "axis", "model"], how="left")
        return long

    ctx["pathways_assoc"] = IO.cached(out + "pathways_assoc_long.parquet", compute, force)
    sig = ctx["pathways_assoc"].query("model == 'frac' and dev_expl_null_pct >= 95")
    print(f"    step8 pathways (assoc): {len(names)} sets; fraction-model above null 95th pct: "
          f"{sig.groupby('axis').size().to_dict()}  {_t(t0)}")
    return ctx["pathways_assoc"]


def save_cells_heavy(ctx):
    """cell_id + heavy-step per-cell columns (labels, neighbourhood features, axes) for the selected cells."""
    df = pd.DataFrame({"cell_id": ctx["cids"][ctx["sel"]], "tile_id": ctx["tile_id"], "d_match": ctx["d"]})
    if "clusters" in ctx:
        df = df.merge(ctx["clusters"], on="cell_id", how="left")
    if "resid_clusters" in ctx:
        df = df.merge(ctx["resid_clusters"], on="cell_id", how="left")
    if "resid_scores" in ctx:
        df = df.merge(ctx["resid_scores"][["cell_id", "resid_pc1_h1", "resid_pc2_h1", "resid_pc3_h1"]], on="cell_id", how="left")
    if "nbhd" in ctx:
        keep = [c for c in ctx["nbhd"].columns if not c.startswith("comp_")] + \
               [f"comp_{t}" for t in ("Malignant cell", "Fibroblast (CAF)", "Endothelial cell", "CD8+ T cell", "Myeloid cell")]
        df = df.merge(ctx["nbhd"][[c for c in keep if c in ctx["nbhd"]]], on="cell_id", how="left")
    IO.save(ctx["out"] + "cells_heavy.parquet", df)
    return df


def run_heavy(ctx, steps=("3", "3b", "4", "7", "8"), force=False):
    nbhd(ctx, force)
    if "3" in steps:
        step3_clustering(ctx, force)
    elif IO.exists(ctx["out"] + "clusters.parquet"):
        ctx["clusters"] = IO.load(ctx["out"] + "clusters.parquet")
    if "3b" in steps:
        step3b_residual(ctx, force)
    if "4" in steps:
        step4_de(ctx, force)
    if "7" in steps:
        step7_association(ctx, force)
    if "8" in steps:
        step8_pathways_assoc(ctx, force)
    save_cells_heavy(ctx)


def load_heavy(ctx):
    """Load heavy-step outputs into ctx if present (after `make out-pull`). Returns the list of missing files."""
    out = ctx["out"]
    files = {"nbhd": "nbhd_features.parquet", "clusters": "clusters.parquet",
             "cluster_pairwise": "clustering_pairwise.parquet", "cyto_only": "cyto_only_clusters.parquet",
             "recovery": "cluster_recovery_total_by_nuc.parquet", "cluster_spatial": "cluster_spatial.parquet",
             "cluster_features": "cluster_features.parquet", "de_nucpart": "de_paired_nucpart.parquet",
             "de_cytopart": "de_paired_cytopart.parquet", "de_overlap": "de_topk_overlap.parquet",
             "assoc_long": "assoc_genes_long.parquet", "assoc_wide": "assoc_genes_wide.parquet",
             "pathways_assoc": "pathways_assoc_long.parquet", "cells_heavy": "cells_heavy.parquet",
             "resid_diag": "resid_diagnostics.json", "resid_summary": "resid_cluster_summary.parquet",
             "resid_moran": "resid_moran.parquet", "resid_scores": "resid_scores.parquet",
             "resid_clusters": "resid_clusters.parquet", "cluster_spatial_strat": "cluster_spatial_stratified.parquet",
             "de_residpart": "de_paired_residpart.parquet"}
    missing = []
    for key, fn in files.items():
        if IO.exists(out + fn):
            ctx[key] = IO.load(out + fn)
        else:
            missing.append(fn)
    return missing


# ====================================================================== step 3b: cytoplasm beyond nucleus
def step3b_residual(ctx, force=False, verbose=True):
    """Residualise cytoplasmic expression on the nuclear design; test residual scores and clusters for spatial
    coherence and boundary/immune association beyond subtype; apply the pre-registered panel rule."""
    import residual as RS
    out, t0 = ctx["out"], time.time()
    if IO.exists(out + "resid_diagnostics.json") and not force:
        ctx["resid_diag"] = IO.load(out + "resid_diagnostics.json")
        for key, fn in (("resid_summary", "resid_cluster_summary.parquet"), ("resid_moran", "resid_moran.parquet"),
                        ("resid_scores", "resid_scores.parquet"), ("resid_clusters", "resid_clusters.parquet"),
                        ("cluster_spatial_strat", "cluster_spatial_stratified.parquet"), ("de_residpart", "de_paired_residpart.parquet")):
            if IO.exists(out + fn):
                ctx[key] = IO.load(out + fn)
        print(f"    [cache] step3b: panel_exists={ctx['resid_diag'].get('panel_exists')}  {_t(t0)}")
        return ctx["resid_diag"]

    mats, gm, genes, cl = ctx["mats"], ctx["gene_mask_cluster"], np.asarray(ctx["genes"]), ctx["clusters"]
    r = C.PRIMARY_RES
    cells = ctx["cells"].loc[ctx["sel"]].reset_index(drop=True)
    d = ctx["d"].astype(float)
    n_perm = C.QUICK_N_PERM_STRAT if ctx["quick"] else C.N_PERM_STRAT
    n_boot = ctx["n_boot"]
    seed = ctx["seed"]
    nn = CL.knn_indices(ctx["xy_sel"])
    feats = ctx["nbhd"].copy()
    feats["log_cell_area"] = np.log(np.maximum(cells["cell_area"].values, 1e-6))
    feats["nuc_cell_ratio"] = cells["nuc_cell_ratio"].values
    F = feats[C.RESID_FEATURES]
    strata = RS.strata_labels(cl[f"nuc_h1_r{r}"].values, d, cells["seg_code"].values)
    gene_names = genes[gm]
    target = float(np.median(d)) / 2.0

    # ---- dense log-normalised halves and the nuclear design
    Yn1 = RS.lognorm_dense(mats["nuc_h1"], gm, target); Yn2 = RS.lognorm_dense(mats["nuc_h2"], gm, target)
    Yc1 = RS.lognorm_dense(mats["cyto_h1"], gm, target); Yc2 = RS.lognorm_dense(mats["cyto_h2"], gm, target)
    sn1, Ln, mn, _ = RS.fit_pca(Yn1, C.RESID_N_PCS_NUC, seed)
    sn2 = RS.project(Yn2, Ln, mn)
    rel_n = RS.pc_reliability(sn1, sn2)
    X, info = RS.design_matrix(sn1, rel_n, cl[f"nuc_h1_r{r}"].values, np.log(np.maximum(d, 1)),
                               np.log(np.maximum(cells["cell_area"].values, 1e-6)), cells["nuc_cell_ratio"].values,
                               cells["seg_code"].values)
    Rc1, Bc1 = RS.residualize(Yc1, X, info); Rc2, _ = RS.residualize(Yc2, X, info); Rn2, _ = RS.residualize(Yn2, X, info)
    del Yn2
    # residual PCA of the cytoplasmic half 1; project half 2; floor from nuclear half 2
    s1, L1, m1, sv1 = RS.fit_pca(Rc1, C.RESID_N_PCS, seed)
    s1, L1, top_genes = RS.orient(s1, L1, gene_names)
    s2 = RS.project(Rc2, L1, m1)
    sf, Lf, mf, svf = RS.fit_pca(Rn2, C.RESID_N_PCS, seed)
    rel_resid = RS.pc_reliability(s1, s2)
    n_dims = RS.n_dims_above_floor(sv1, svf)
    cc = RS.canonical_correlations(s1[:, :10], sn2[:, :10])
    del Rc1, Rc2, Rn2
    # converse: nuclear half 1 residualised on the cytoplasmic design
    sc1, Lc, mc, _ = RS.fit_pca(Yc1, C.RESID_N_PCS_NUC, seed)
    rel_c = RS.pc_reliability(sc1, RS.project(Yc2, Lc, mc))
    Xc, infoc = RS.design_matrix(sc1, rel_c, cl[f"cyto_h1_r{r}"].values, np.log(np.maximum(d, 1)),
                                 np.log(np.maximum(cells["cell_area"].values, 1e-6)), cells["nuc_cell_ratio"].values,
                                 cells["seg_code"].values)
    Rn1c, _ = RS.residualize(Yn1, Xc, infoc)
    sconv, _, _, svconv = RS.fit_pca(Rn1c, C.RESID_N_PCS, seed)
    del Yn1, Yc1, Yc2, Rn1c
    if verbose:
        print(f"      nuclear PCs kept {int(info['is_pc'].sum())}/{len(rel_n)} (rel PC1-3 {np.round(rel_n[:3], 2)}), "
              f"used after partialling {int(info['used'][info['is_pc']].sum())}, EIV shrink {np.nanmean(info['eiv_shrink']):.2f}; residual dims above floor {n_dims}; "
              f"resid score rel PC1-3 {np.round(rel_resid[:3], 2)}; max canonical corr with nuc_h2 {cc.max():.2f}")

    # ---- continuous score tests (PCs 1-3)
    moran_rows = []
    for k in range(C.RESID_N_SCORE_PCS):
        mt = RS.morans_test(s1[:, k], nn, strata, n_perm, seed + k)
        mt_f = RS.morans_test(sf[:, k], nn, strata, n_perm, seed + k)
        noise = RS.morans_i(s1[:, k] - s2[:, k], nn)
        row = {"pc": k + 1, "top_gene": top_genes[k], "reliability": float(2 * rel_resid[k] / (1 + rel_resid[k])),
               "reliability_half": float(rel_resid[k]), "singular_value": float(sv1[k]), "floor_singular_value_max": float(svf.max()),
               "I": mt["I"], "I_perm_p99": mt["I_perm_p99"], "I_p": mt["p_perm"], "I_floor": mt_f["I"], "I_noise_control": noise}
        for f in ("tumor_frac", "immune_kernel", "log1p_dist_nontumor"):
            for half, sc_ in (("h1", s1[:, k]), ("h2", s2[:, k])):
                dr = RS.dose_response(sc_, F[f].values, strata, ctx["tile_id"], n_boot=n_boot, seed=seed + k)
                row[f"dr_{f}_{half}"] = dr["effect"]; row[f"dr_{f}_{half}_lo"] = dr["effect_lo"]; row[f"dr_{f}_{half}_hi"] = dr["effect_hi"]
                if half == "h1":
                    for q in range(len(dr["bin_means"])):
                        row[f"dr_{f}_q{q + 1}"] = dr["bin_means"][q]
        moran_rows.append(row)
    moran = pd.DataFrame(moran_rows)
    pa = np.zeros(len(moran), dtype=bool)
    moran["A_feature"] = ""
    for pf in C.RESID_INTERFACE_FEATURES:
        ok = ((moran[f"dr_{pf}_h1"].abs() >= C.RESID_EFFECT_SD) & (moran[f"dr_{pf}_h2"].abs() >= C.RESID_EFFECT_SD)
              & ((moran[f"dr_{pf}_h1_lo"] > 0) | (moran[f"dr_{pf}_h1_hi"] < 0))
              & ((moran[f"dr_{pf}_h2_lo"] > 0) | (moran[f"dr_{pf}_h2_hi"] < 0))
              & (np.sign(moran[f"dr_{pf}_h1"]) == np.sign(moran[f"dr_{pf}_h2"])))
        moran.loc[ok & ~pa, "A_feature"] = pf
        pa |= ok.values
    moran["passes_A"] = ((moran["reliability"] >= C.RESID_SCORE_REL_MIN) & (moran["I"] >= C.RESID_MORAN_MIN)
                         & (moran["I"] > moran["I_perm_p99"]) & (moran["I"] > moran["I_floor"]) & pa)
    moran["passes_A"] &= n_dims >= 1
    IO.save(out + "resid_moran.parquet", moran)

    # ---- residual clusters (display + gate), floor and converse partitions, nested stats for cyto_m
    lab1 = CL.cluster_from_pcs(s1, seed=seed)
    lab2 = CL.cluster_from_pcs(s2, seed=seed)
    labf = CL.cluster_from_pcs(sf, seed=seed)
    labconv = CL.cluster_from_pcs(sconv, seed=seed)
    clusters = pd.DataFrame({"cell_id": cells["cell_id"].values, "strata": strata})
    for nm, lab in (("resid_cyto_h1", lab1), ("resid_cyto_h2", lab2), ("resid_floor_nuc_h2", labf), ("resid_converse_nuc_h1", labconv)):
        for res, l in lab.items():
            clusters[f"{nm}_r{res}"] = l
    IO.save(out + "resid_clusters.parquet", clusters)
    scores = pd.DataFrame({"cell_id": cells["cell_id"].values})
    for k in range(5):
        scores[f"resid_pc{k + 1}_h1"] = s1[:, k]; scores[f"resid_pc{k + 1}_h2"] = s2[:, k]; scores[f"floor_pc{k + 1}"] = sf[:, k]
    IO.save(out + "resid_scores.parquet", scores)
    IO.save(out + "resid_loadings.npz", {"genes": gene_names, "loadings_cyto_h1": L1, "sv_cyto_h1": sv1, "sv_floor": svf,
                                          "sv_converse": svconv, "nuc_rel": rel_n, "resid_rel_half": rel_resid,
                                          "design_columns": info["column"].values.astype(str), "B_cyto_h1": Bc1})

    def cluster_stats(lab, name):
        p = RS.stratified_purity(lab, nn, strata, n_perm, seed)
        e = RS.stratified_effects(lab, F, strata, ctx["tile_id"], n_boot=n_boot, n_perm=n_perm, seed=seed)
        w = e.pivot(index="cluster", columns="feature", values="d"); w.columns = [f"d_{c}" for c in w.columns]
        lo = e.pivot(index="cluster", columns="feature", values="d_lo"); lo.columns = [f"d_{c}_lo" for c in lo.columns]
        hi = e.pivot(index="cluster", columns="feature", values="d_hi"); hi.columns = [f"d_{c}_hi" for c in hi.columns]
        pp = e.pivot(index="cluster", columns="feature", values="p_perm"); pp.columns = [f"d_{c}_p" for c in pp.columns]
        s = p.set_index("cluster").join([w, lo, hi, pp]).reset_index()
        dep = pd.DataFrame({"cluster": lab, "d": d, "seg": cells["seg_code"].values}).groupby("cluster")
        s["median_depth"] = dep["d"].median().reindex(s["cluster"]).values
        s["depth_ratio"] = s["median_depth"] / float(np.median(d))
        s["seg_share"] = dep["seg"].agg(lambda x: x.value_counts(normalize=True).iloc[0]).reindex(s["cluster"]).values
        s.insert(0, "partition", name)
        return s

    l1, l2 = lab1[r], lab2[r]
    summ = cluster_stats(l1, "resid_cyto_h1")
    summ_h2 = cluster_stats(l2, "resid_cyto_h2")
    rec = CL.cluster_recovery(l1, l2)
    summ["repro_f1"] = rec["best_f1"].values; summ["matched_h2"] = rec["best_other"].values
    for pf in C.RESID_INTERFACE_FEATURES:
        summ[f"d_{pf}_h2"] = summ_h2.set_index("cluster")[f"d_{pf}"].reindex(summ["matched_h2"]).values
    floor_s = cluster_stats(labf[r], "resid_floor_nuc_h2")
    conv_s = cluster_stats(labconv[r], "resid_converse_nuc_h1")
    nested = cluster_stats(cl[f"cyto_m_r{r}"].values, "cyto_m")
    IO.save(out + "cluster_spatial_stratified.parquet", pd.concat([nested, floor_s, conv_s, summ_h2], ignore_index=True))
    # DE link on the independent halves, within strata
    mean_sel = np.asarray((ctx["Xn"][ctx["sel"]] + ctx["Xc"][ctx["sel"]]).mean(0)).ravel()
    de = DE.stratified_paired_logfc(mats["nuc_h2"], mats["cyto_h2"], l1, strata, genes, ctx["tile_id"],
                                    gene_mask=mean_sel >= C.MEAN_MIN_TILE, n_boot=n_boot, seed=seed)
    de = DE.call_de(de)
    IO.save(out + "de_paired_residpart.parquet", de)
    summ["n_cyto_only_de"] = de.groupby("cluster")["cyto_only"].sum().reindex(summ["cluster"]).fillna(0).astype(int).values
    summ["n_nuc_only_de"] = de.groupby("cluster")["nuc_only"].sum().reindex(summ["cluster"]).fillna(0).astype(int).values
    sample_seg_share = float(cells["seg_code"].value_counts(normalize=True).iloc[0])
    summ = RS.eligibility(summ, len(cells), sample_seg_share, float(floor_s["purity_excess"].max()))
    summ["n_dims_above_floor"] = n_dims
    IO.save(out + "resid_cluster_summary.parquet", summ)

    # ---- the panel rule
    passes_A = bool(moran["passes_A"].any())
    passes_B = bool(summ["eligible"].any())
    chosen = int(summ.loc[summ["rank"] == 1, "cluster"].iloc[0]) if passes_B else None
    zoom = RS.zoom_window(ctx["xy_sel"], l1 == chosen) if chosen is not None else None
    diag = {"dataset": ctx["ds"], "n_cells": int(len(cells)), "n_strata": int(strata.max() + 1),
            "frac_cells_in_strata": float((strata >= 0).mean()),
            "nuclear_pcs_kept": int(info["is_pc"].sum()), "nuclear_pc_reliability": [float(x) for x in rel_n],
            "eiv_applied": bool(info["eiv_applied"].iloc[0]), "eiv_mean_shrink": float(np.nanmean(info["eiv_shrink"])),
            "nuclear_pcs_used_after_partialling": int(info["used"][info["is_pc"]].sum()),
            "resid_score_reliability_half_max": float(np.max(rel_resid)),
            "n_dims_above_floor": n_dims, "sv_cyto_resid_top5": [float(x) for x in sv1[:5]], "sv_floor_max": float(svf.max()),
            "sv_converse_top5": [float(x) for x in svconv[:5]],
            "max_canonical_corr_resid_vs_nuc_h2": float(cc.max()),
            "ari_resid_vs_nuc_h2": float(CL.compare_partitions(l1, cl[f"nuc_h2_r{r}"].values)["ari"]),
            "ari_resid_h1_vs_h2": float(CL.compare_partitions(l1, l2)["ari"]),
            "resid_score_reliability": [float(x) for x in moran["reliability"]],
            "passes_A_continuous": passes_A, "passes_B_cluster": passes_B,
            "panel_exists": passes_A and passes_B,
            "verdict": ("panel" if passes_A and passes_B else "continuous axis only" if passes_A else
                        "cluster only (score fails)" if passes_B else "no cytoplasmic structure beyond nucleus"),
            "chosen_cluster": chosen, "chosen_label": (str(summ.loc[summ["rank"] == 1, "label"].iloc[0]) if chosen is not None else None),
            "chosen_feature": (str(summ.loc[summ["rank"] == 1, "effect_feature"].iloc[0]) if chosen is not None else None),
            "A_feature": [str(x) for x in moran["A_feature"]],
            "chosen_rank_stat": (float(summ.loc[summ["rank"] == 1, "rank_stat"].iloc[0]) if chosen is not None else None),
            "zoom_window": zoom, "n_eligible_clusters": int(summ["eligible"].sum()),
            "n_clusters_resid": int(l1.max() + 1), "n_clusters_floor": int(labf[r].max() + 1)}
    IO.save(out + "resid_diagnostics.json", diag)
    ctx.update(resid_diag=diag, resid_summary=summ, resid_moran=moran, resid_scores=scores, resid_clusters=clusters,
               cluster_spatial_strat=pd.concat([nested, floor_s, conv_s, summ_h2], ignore_index=True), de_residpart=de)
    print(f"    step3b residual: verdict '{diag['verdict']}'; dims above floor {n_dims}; PC1 I {moran['I'].iloc[0]:.3f} "
          f"(floor {moran['I_floor'].iloc[0]:.3f}, p99 {moran['I_perm_p99'].iloc[0]:.3f}); eligible clusters "
          f"{diag['n_eligible_clusters']}/{diag['n_clusters_resid']}; ARI resid h1/h2 {diag['ari_resid_h1_vs_h2']:.2f}  {_t(t0)}")
    return diag


# ====================================================================== step 0b: annotation QC (report only)
def transfer_model(force=False):
    """The non-epithelial label classifier trained on the 10x-labelled datasets, cached across datasets."""
    import annotation_qc as AQ
    xo = IO.out_dir("cross_sample")
    p = xo + "annotation_transfer_model.pkl"
    if IO.exists(p) and not force:
        import pickle
        return pickle.load(open(p, "rb"))
    t0 = time.time()
    clf = AQ.fit_transfer()
    import pickle
    pickle.dump(clf, open(p, "wb"))
    ref = IO.cached(xo + "annotation_transfer_lodo.parquet", AQ.lodo_reference, force)
    print(f"    transfer classifier fitted on {C.TRANSFER_TRAIN}; LODO accuracy median {ref['accuracy'].median():.2f}  {_t(t0)}")
    return clf


def step0b_annotation_qc(ctx, force=False, clf=None):
    """Marker-set scores, class profiles, tumor-set purity flags, malignant-cluster profiles, label transfer."""
    import annotation_qc as AQ
    out, t0, ds = ctx["out"], time.time(), ctx["ds"]
    if IO.exists(out + "annotation_qc_flags.json") and not force:
        ctx["annot_flags"] = IO.load(out + "annotation_qc_flags.json")
        ctx["annot_classes"] = IO.load(out + "annotation_qc_class_scores.parquet")
        ctx["annot_clusters"] = IO.load(out + "annotation_qc_tumor_clusters.parquet") if IO.exists(out + "annotation_qc_tumor_clusters.parquet") else pd.DataFrame()
        ctx["annot_transfer"] = IO.load(out + "annotation_transfer.parquet")
        ctx["annot_provenance"] = IO.load(out + "annotation_provenance.json")
        if IO.exists(out + "annotation_qc_cells.parquet"):
            ctx["annot_cells"] = IO.load(out + "annotation_qc_cells.parquet")
        print(f"    [cache] step0b annotation QC  {_t(t0)}")
        return ctx["annot_flags"]
    Y, obs, genes = AQ.load_all_cells(ds)
    scores = AQ.marker_scores(Y, genes, ds)
    prov = AQ.provenance(ds, obs)
    classes = AQ.class_profile(scores, obs)
    flags, summary = AQ.purity_flags(scores, obs, ds)
    clusters = AQ.tumor_cluster_profiles(scores, obs, ds)
    clf = clf or transfer_model()
    transfer, pred, pmax = AQ.apply_transfer(clf, Y, obs, ds)
    summary["marker_sets"] = scores.attrs["sets"]
    mal = obs["cell_type_merged"].values == C.TUMOR_TYPE
    nonepi = (pmax >= C.TRANSFER_CONF) & (pred != C.TRANSFER_EPITHELIAL_CLASS)
    summary["frac_malignant_confident_nonepithelial"] = float((nonepi & mal).sum() / max(mal.sum(), 1))
    summary["frac_malignant_confident_epithelial"] = float(((pmax >= C.TRANSFER_CONF) & (pred == C.TRANSFER_EPITHELIAL_CLASS) & mal).sum() / max(mal.sum(), 1))
    per_cell = pd.concat([obs[["cell_id", "cell_type_merged"]].reset_index(drop=True), scores, flags], axis=1)
    per_cell["transfer_pred"] = pred; per_cell["transfer_pmax"] = pmax
    IO.save(out + "annotation_qc_cells.parquet", per_cell)
    IO.save(out + "annotation_qc_class_scores.parquet", classes)
    IO.save(out + "annotation_qc_flags.json", summary)
    IO.save(out + "annotation_provenance.json", prov)
    if len(clusters):
        IO.save(out + "annotation_qc_tumor_clusters.parquet", clusters)
    IO.save(out + "annotation_transfer.parquet", transfer)
    ctx.update(annot_flags=summary, annot_classes=classes, annot_clusters=clusters, annot_transfer=transfer,
               annot_provenance=prov, annot_cells=per_cell)
    print(f"    step0b annotation QC: {prov['source']}; malignant flags lineage-neg {summary['frac_lineage_negative']:.3f}, "
          f"immune-like {summary['frac_immune_like']:.3f}, confident non-epithelial by transfer "
          f"{summary['frac_malignant_confident_nonepithelial']:.3f}; benign-like malignant clusters "
          f"{int(clusters['benign_like'].sum()) if len(clusters) else 'n/a'}  {_t(t0)}")
    return summary
