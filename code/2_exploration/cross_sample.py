"""Per-dataset candidate calls (step 9) and the cross-sample section (B6): recurrence of per-gene patterns,
pairwise agreement of localisation between samples, recurrence of cytoplasmic clusters, final candidate table."""

import numpy as np
import pandas as pd

import config as C
import io_utils as IO


# ------------------------------------------------------------------ per-dataset candidates
def candidates(ctx):
    """One row per gene with boolean criteria and a score; requires the light steps and (if present) heavy results."""
    gc = ctx["genes_compartment"].set_index("gene")
    ct = ctx.get("conc_tile")
    cc = ctx.get("conc_cell")
    df = pd.DataFrame(index=gc.index)
    df["mean_total"] = gc["mean_total"]
    df["beta"] = gc["beta"]; df["class"] = gc["class"]
    df["rel_param"] = gc["rel_param"]; df["rel_param_pct"] = gc.get("rel_param_pct")
    df["reliable"] = (gc["rel_param"] >= C.REL_CANDIDATE) & (gc["n_units_rel"] >= C.MIN_UNITS)
    df["ctrl_rel_top"] = gc.get("rel_param_pct", pd.Series(np.nan, index=gc.index)) >= C.CTRL_PERCENTILE
    df["ctrl_seg_ok"] = ~(gc.get("beta_seg_range_pct", pd.Series(np.nan, index=gc.index)) >= C.CTRL_PERCENTILE)
    df["localized"] = gc["class"].isin(["nuclear-retained", "cytoplasm-enriched"])
    for name, tab in (("tile", ct), ("cell", cc)):
        if tab is None:
            continue
        t = tab.set_index("gene")
        df[f"r_true_{name}"] = t["r_true"].reindex(df.index)
        df[f"divergent_{name}"] = (t["r_true"] <= C.R_TRUE_DIVERGENT).reindex(df.index).fillna(False) & \
                                  (t.get("r_true_low_pct", pd.Series(np.nan, index=t.index)) >= C.CTRL_PERCENTILE).reindex(df.index).fillna(False)
    df["divergent"] = df.get("divergent_tile", False) | df.get("divergent_cell", False)
    if "de_nucpart" in ctx:
        d = ctx["de_nucpart"]
        df["cyto_only_de"] = d.groupby("gene")["cyto_only"].any().reindex(df.index).fillna(False)
        df["nuc_only_de"] = d.groupby("gene")["nuc_only"].any().reindex(df.index).fillna(False)
        df["max_abs_delta"] = d.groupby("gene")["delta"].apply(lambda s: s.abs().max()).reindex(df.index)
    else:
        df["cyto_only_de"] = False; df["nuc_only_de"] = False
    sens_cols = []
    if "assoc_wide" in ctx:
        w = ctx["assoc_wide"].set_index("gene")
        for a in C.CANDIDATE_AXES:
            f, cmn = f"{a}__frac", f"{a}__cyto_minus_nuc"
            if f in w:
                df[f"dev_frac_{a}"] = w[f].reindex(df.index)
                df[f"dev_cmn_{a}"] = w[cmn].reindex(df.index) if cmn in w else np.nan
                df[f"sens_{a}"] = ((w[f] >= C.ASSOC_MIN_EXCESS) | (w.get(cmn, pd.Series(np.nan, index=w.index)) >= C.ASSOC_MIN_EXCESS)).reindex(df.index).fillna(False)
                sens_cols.append(f"sens_{a}")
        for a in ("niche", "tumor_frac_bin", "immune_frac_bin"):
            col = f"{a}_beyond_subtype__frac"
            if col in w:
                df[f"beyond_subtype_{a}"] = w[col].reindex(df.index)
    df["sens_any"] = df[sens_cols].any(axis=1) if sens_cols else False
    df["n_sens_axes"] = df[sens_cols].sum(axis=1) if sens_cols else 0
    # annotation
    M, info = ctx.get("pw_M"), ctx.get("pw_info")
    if M is not None:
        Mc = M.tocsr()
        ann = []
        for i in range(M.shape[0]):
            js = Mc.indices[Mc.indptr[i]:Mc.indptr[i + 1]]
            ann.append(";".join(info["pathway"].values[js]))
        df["pathways"] = pd.Series(ann, index=ctx["genes"]).reindex(df.index)
        df["in_SG"] = df["pathways"].str.contains("SG_markers").fillna(False)
        df["stress_annotated"] = df["pathways"].str.contains("STRESS|UNFOLDED|HYPOXIA|HEAT|OXIDATIVE|APOPTOSIS|SG_markers", case=False).fillna(False)
    df["score"] = (df["reliable"].astype(int) + df["divergent"].astype(int) + df["cyto_only_de"].astype(int)
                   + df["sens_any"].astype(int) + df["ctrl_rel_top"].astype(int))
    df["candidate"] = df["reliable"] & df["ctrl_seg_ok"] & (df["divergent"] | df["cyto_only_de"] | df["sens_any"])
    df = df.reset_index().rename(columns={"index": "gene"})
    df.insert(1, "dataset", C.short(ctx["ds"]))
    return df.sort_values(["candidate", "score", "rel_param"], ascending=[False, False, False])


# ------------------------------------------------------------------ cross-sample
def merge_tables(dss, fn, suffix=""):
    rows = []
    for ds in dss:
        p = IO.out_dir(ds + suffix) + fn
        if IO.exists(p):
            t = IO.load(p); t["dataset"] = C.short(ds); rows.append(t)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def recurrence(cand_all):
    """Per gene: in how many samples each criterion holds; recurrent = candidate in >= MIN_SAMPLES_RECURRENT."""
    flags = ["candidate", "reliable", "divergent", "cyto_only_de", "sens_any", "localized"]
    flags += [c for c in cand_all.columns if c.startswith("sens_") and c not in flags]
    g = cand_all.groupby("gene")
    rec = g[flags].sum().astype(int)
    rec.columns = ["n_" + c for c in rec.columns]
    rec["n_samples_tested"] = g.size()
    rec["n_samples_reliable"] = rec["n_reliable"]
    rec["mean_beta"] = g["beta"].mean(); rec["sd_beta"] = g["beta"].std()
    rec["mean_rel"] = g["rel_param"].mean()
    rec["classes"] = g["class"].apply(lambda s: ";".join(sorted(set(s))))
    rec["datasets_candidate"] = cand_all[cand_all["candidate"]].groupby("gene")["dataset"].apply(lambda s: ";".join(sorted(s)))
    for c in ("pathways", "in_SG", "stress_annotated"):
        if c in cand_all:
            rec[c] = g[c].first()
    rec["recurrent"] = rec["n_candidate"] >= C.MIN_SAMPLES_RECURRENT
    return rec.reset_index().sort_values(["n_candidate", "n_sens_any", "mean_rel"], ascending=False)


def gene_by_sample(cand_all, col):
    return cand_all.pivot_table(index="gene", columns="dataset", values=col)


def cluster_recurrence(dss, suffix="", min_r=0.8):
    """Correlate cytoplasmic-partition cluster centroids between samples; best partner per cluster."""
    cen = {}
    for ds in dss:
        p = IO.out_dir(ds + suffix) + "centroids_cyto_m.npz"
        if IO.exists(p):
            z = IO.load(p); cen[C.short(ds)] = z["centroids"]
    rows = []
    names = list(cen)
    for i, a in enumerate(names):
        for b in names:
            if a == b:
                continue
            A, B = cen[a], cen[b]
            # centre genes within each sample to remove sample-wide shifts
            A = A - A.mean(0, keepdims=True); B = B - B.mean(0, keepdims=True)
            A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
            num = A @ B.T
            den = np.outer(np.sqrt((A ** 2).sum(axis=1)), np.sqrt((B ** 2).sum(axis=1)))
            R = num / np.maximum(den, 1e-12)
            for ca in range(R.shape[0]):
                j = int(np.argmax(R[ca]))
                rows.append({"sample": a, "cluster": ca, "other_sample": b, "best_cluster": j, "r": float(R[ca, j])})
    df = pd.DataFrame(rows)
    if len(df):
        df["matched"] = df["r"] >= min_r
    return df


def decision_summary(dss, suffix=""):
    """The numbers the design decisions rest on, one row per sample."""
    rows = []
    for ds in dss:
        o = IO.out_dir(ds + suffix)
        r = {"dataset": C.short(ds)}
        if IO.exists(o + "qc_summary.json"):
            q = IO.load(o + "qc_summary.json")
            r.update(n_cells=q["n_tumor_cells"], frac_keep=round(q["frac_keep"], 3), median_d_match=q["median_d_match"],
                     nuclear_share=round(q["nuclear_share_of_reads"], 3), genes_ge_0_5=q["genes_mean_ge_0.5"])
        if IO.exists(o + "genes_compartment.parquet"):
            g = IO.load(o + "genes_compartment.parquet")
            r.update(n_nuclear_retained=int((g["class"] == "nuclear-retained").sum()),
                     n_cyto_enriched=int((g["class"] == "cytoplasm-enriched").sum()),
                     n_reliable=int(((g["rel_param"] >= C.REL_CANDIDATE) & (g["n_units_rel"] >= C.MIN_UNITS)).sum()),
                     median_rel=round(float(g["rel_param"].median()), 3))
        if IO.exists(o + "controls_sets.parquet"):
            s = IO.load(o + "controls_sets.parquet").set_index("set")
            r.update(SG_rel=round(float(s.loc["SG_markers", "rel_param"]), 3), ctrl_rel=round(float(s.loc["matched_control", "rel_param"]), 3))
        if IO.exists(o + "clustering_pairwise.parquet"):
            p = IO.load(o + "clustering_pairwise.parquet"); p = p[p["resolution"] == C.PRIMARY_RES].set_index(["a", "b"])["ari"]
            r.update(ari_nuc_cyto=round(p.get(("nuc_m", "cyto_m"), np.nan), 3), ari_nuc_total=round(p.get(("nuc_m", "total_full"), np.nan), 3),
                     ari_nuc_halves=round(p.get(("nuc_h1", "nuc_h2"), np.nan), 3), ari_cyto_halves=round(p.get(("cyto_h1", "cyto_h2"), np.nan), 3),
                     ari_cross_halves=round(p.get(("nuc_h1", "cyto_h1"), np.nan), 3))
        if IO.exists(o + "cyto_only_clusters.parquet"):
            r["n_cyto_only_clusters"] = int(IO.load(o + "cyto_only_clusters.parquet")["cyto_only"].sum())
        if IO.exists(o + "cluster_recovery_total_by_nuc.parquet"):
            r["total_recovered_by_nuc_minF1"] = round(float(IO.load(o + "cluster_recovery_total_by_nuc.parquet")["best_f1"].min()), 3)
        if IO.exists(o + "de_paired_nucpart.parquet"):
            d = IO.load(o + "de_paired_nucpart.parquet")
            r.update(n_cyto_only_de=int(d.groupby("gene")["cyto_only"].any().sum()), n_nuc_only_de=int(d.groupby("gene")["nuc_only"].any().sum()))
        if IO.exists(o + "assoc_genes_wide.parquet"):
            w = IO.load(o + "assoc_genes_wide.parquet")
            for a in ("subtype_nuc", "niche", "tumor_frac_bin"):
                if f"{a}__frac" in w:
                    r[f"n_sens_{a}"] = int((w[f"{a}__frac"] >= C.ASSOC_MIN_EXCESS).sum())
            if "niche_beyond_subtype__frac" in w:
                r["n_niche_beyond_subtype"] = int((w["niche_beyond_subtype__frac"] >= C.ASSOC_MIN_EXCESS).sum())
        if IO.exists(o + "candidates.parquet"):
            r["n_candidates"] = int(IO.load(o + "candidates.parquet")["candidate"].sum())
        if IO.exists(o + "resid_diagnostics.json"):
            dg = IO.load(o + "resid_diagnostics.json")
            r.update(resid_verdict=dg["verdict"], resid_dims_above_floor=dg["n_dims_above_floor"],
                     resid_eligible_clusters=dg["n_eligible_clusters"], resid_chosen_cluster=dg["chosen_cluster"],
                     resid_chosen_label=dg["chosen_label"], resid_rank_stat=dg["chosen_rank_stat"])
            if IO.exists(o + "resid_moran.parquet"):
                mo = IO.load(o + "resid_moran.parquet")
                r.update(resid_pc1_I=round(float(mo["I"].iloc[0]), 3), resid_pc1_I_floor=round(float(mo["I_floor"].iloc[0]), 3),
                         resid_pc1_rel=round(float(mo["reliability"].iloc[0]), 3),
                         resid_pc1_dr_tumor_frac=round(float(mo[f"dr_{C.RESID_PRIMARY_FEATURE}_h1"].iloc[0]), 3))
        rows.append(r)
    return pd.DataFrame(rows)


def residual_axis_recurrence(dss, suffix=""):
    """Correlation of the cytoplasmic residual PC1 loadings between samples over their shared genes."""
    load = {}
    for ds in dss:
        p = IO.out_dir(ds + suffix) + "resid_loadings.npz"
        if IO.exists(p):
            z = IO.load(p); load[C.short(ds)] = pd.Series(z["loadings_cyto_h1"][0], index=z["genes"].astype(str))
    names = list(load)
    if len(names) < 2:
        return pd.DataFrame()
    M = pd.DataFrame(np.nan, index=names, columns=names)
    for a in names:
        for b in names:
            common = load[a].index.intersection(load[b].index)
            M.loc[a, b] = float(np.corrcoef(load[a].loc[common], load[b].loc[common])[0, 1]) if len(common) > 10 else np.nan
    return M.reset_index().rename(columns={"index": "sample"})


def provenance_table(dss, suffix=""):
    rows = []
    for ds in dss:
        o = IO.out_dir(ds + suffix)
        if IO.exists(o + "annotation_provenance.json"):
            r = IO.load(o + "annotation_provenance.json")
            if IO.exists(o + "annotation_qc_flags.json"):
                f = IO.load(o + "annotation_qc_flags.json")
                r.update({k: v for k, v in f.items() if k.startswith("frac_")})
            if IO.exists(o + "annotation_qc_tumor_clusters.parquet"):
                t = IO.load(o + "annotation_qc_tumor_clusters.parquet")
                r["benign_like_clusters"] = int(t["benign_like"].sum()); r["benign_like_cells_frac"] = float(t.loc[t["benign_like"], "n_cells"].sum() / t["n_cells"].sum())
            rows.append(r)
    return pd.DataFrame(rows)
