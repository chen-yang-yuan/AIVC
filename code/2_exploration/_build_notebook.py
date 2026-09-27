"""Builds fig1_exploration.ipynb (run once; the notebook is then edited/run in Jupyter)."""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# Fig. 1 exploration: nuclear vs cytoplasmic expression in tumor cells

Six Xenium 5K FFPE samples (BC, OC, CC, LC, Prostate, Skin), tumor cells only. Narrative: (A) the compartments differ,
(B) the difference is biology rather than artifact, (C) the difference relates to the model's inputs. Design and
status: `plans/fig1_exploration.md`. Execution order and file inventory: `README.md` in this folder.

**Where things run.** Light steps (0, 1, 2, 5, 8-light, 9, cross-sample) run here on all cells. Heavy steps (3
clustering, 4 DE, 7 association, 8 pathway association) run on HGCC via `1_heavy.py` / `1_heavy.sh` and their tables are
read here after `make out-pull`. With `QUICK = True` the heavy steps also run here on a 20k-cell subsample (outputs go
to `output/2_exploration/<ds>_quick/`), which is how to develop and check the pipeline locally.""")

code("""import warnings; warnings.filterwarnings("ignore")
import json, time, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import config as C, io_utils as IO, compartment as CP, clustering as CL, de as DE, neighborhood as NB
import association as AS, controls as CT, plotting as PL, pipeline as P, cross_sample as XS

pd.set_option("display.width", 200); pd.set_option("display.max_columns", 60)

# ------------------------------------------------------------------ run settings
RUN_DATASETS = C.DATASETS            # narrow while iterating, e.g. ["Xenium_5K_LC"]
QUICK = False                        # True: 20k-cell subsample, heavy steps run locally, outputs in <ds>_quick/
RUN_HEAVY_LOCALLY = QUICK            # run steps 3/4/7/8 here (only sensible with QUICK); else read HGCC results
FORCE = dict(annotation=False, compartment=False, concordance=False, controls=False, pathways=False,
             clustering=False, residual=False, de=False, association=False, pathways_assoc=False, figures=True)
SUFFIX = "_quick" if QUICK else \"\"
print("datasets:", [C.short(d) for d in RUN_DATASETS], "| quick:", QUICK)""")

md("""## Per-dataset loop

Each iteration prepares one dataset (depths, cell filters, matched-depth matrices, half-splits, tiles, offsets), runs
the light steps, loads or runs the heavy steps, draws the figures and writes `candidates.parquet`.""")

code("""CTX = {}
for ds in RUN_DATASETS:
    t0 = time.time()
    print(f"========== {ds} ==========")
    ctx = P.prepare(ds, quick=QUICK)
    P.step0_qc(ctx)
    # ---- Step 0b: annotation QC report (all cells; report only, the tumor definition is unchanged)
    P.step0b_annotation_qc(ctx, FORCE["annotation"])

    # ---- Part A1: global description and per-gene localisation (all cells)
    P.step1_compartment(ctx, FORCE["compartment"])
    # ---- Part A2: split-half concordance (cell level on kept cells, tile level)
    P.step2_concordance(ctx, FORCE["concordance"])
    # ---- Part B5: abundance-bin percentiles, strata, SG vs matched control
    P.step5_controls(ctx, FORCE["controls"])
    # ---- Part C8 (light): pathway-level localisation, reliability, concordance with matched random-set nulls
    P.step8_pathways_light(ctx, FORCE["pathways"])

    # ---- heavy steps: A3 clustering, A4 DE, C7 association, C8 pathway association
    if RUN_HEAVY_LOCALLY:
        P.nbhd(ctx)
        P.step3_clustering(ctx, FORCE["clustering"])
        P.step3b_residual(ctx, FORCE["residual"])
        P.step4_de(ctx, FORCE["de"])
        P.step7_association(ctx, FORCE["association"])
        P.step8_pathways_assoc(ctx, FORCE["pathways_assoc"])
        P.save_cells_heavy(ctx)
    else:
        missing = P.load_heavy(ctx)
        if missing:
            print(f"    [heavy] missing {len(missing)} files (run 1_heavy.sh on HGCC, then `make out-pull`): {missing[:4]} ...")

    # ---- Part C9: candidate genes for this dataset
    cand = XS.candidates(ctx)
    IO.save(ctx["out"] + "candidates.parquet", cand)
    ctx["candidates"] = cand
    print(f"    step9 candidates: {int(cand['candidate'].sum())} genes "
          f"(reliable {int(cand['reliable'].sum())}, divergent {int(cand['divergent'].sum())}, "
          f"cyto-only DE {int(cand['cyto_only_de'].sum())}, axis-sensitive {int(cand['sens_any'].sum())})")
    CTX[ds] = ctx
    print(f"    done in {(time.time() - t0) / 60:.1f} min")""")

md("""## Per-dataset figures

Panel order follows the tracker: read fractions and positive controls (A1), concordance (A2), three-way clustering
with spatial maps (A3), DE comparison (A4), association with the input axes (C7), pathway summary (C8).""")

code("""def figures(ctx):
    ds, fig, cells = ctx["ds"], ctx["fig"], ctx["cells"]
    g = ctx["genes_compartment"]
    # A0: annotation QC
    PL.marker_score_heatmap(ctx["annot_classes"], fig + "A0_marker_scores.jpeg", ds)
    PL.tumor_cluster_heatmap(ctx.get("annot_clusters"), fig + "A0_tumor_cluster_profiles.jpeg", ds)
    if "annot_cells" in ctx:
        ac = ctx["annot_cells"]; mal = ac["cell_type_merged"] == "Malignant cell"
        flagged = (ac.loc[mal, ["lineage_negative", "immune_like", "endothelial_like"]].any(axis=1)).values
        xy_all = cells[["global_x", "global_y"]].values
        PL.spatial_map(xy_all, np.where(flagged, "flagged", "ok"), fig + "A0_spatial_flags.jpeg", categorical=True,
                       title=f"{C.short(ds)}: tumor cells flagged by marker QC ({flagged.mean():.1%})")
    # A1
    PL.read_fractions(cells, fig + "A1_read_fractions.jpeg", ds)
    ctrl_nuc = ctx["qc"]["pos_ctrl_nuclear_present"]
    ctrl_cyto = ctx["qc"]["pos_ctrl_cyto_present"] + ctx["qc"]["top_abundant_genes"]
    PL.logor_vs_abundance(g, fig + "A1_logor_vs_abundance.jpeg", ds, label_genes=ctrl_nuc + ctrl_cyto[:6])
    PL.positive_controls(g, fig + "A1_positive_controls.jpeg", ds, ctrl_nuc, ctrl_cyto)
    PL.spatial_map(cells[["global_x", "global_y"]].values, cells["nuc_frac"].values, fig + "A1_spatial_nuc_frac.jpeg",
                   title=f"{C.short(ds)}: per-cell nuclear fraction", vmin=0.1, vmax=0.9)
    # A2
    top_div = ctx["conc_tile"].dropna(subset=["r_true"]).nsmallest(8, "r_true")["gene"].tolist()
    PL.concordance_scatter(ctx["conc_cell"], ctx["conc_tile"], fig + "A2_concordance_scatter.jpeg", ds, highlight=top_div)
    # A3
    if "clusters" in ctx:
        pw = ctx["cluster_pairwise"]; pw = pw[pw["resolution"] == C.PRIMARY_RES]
        PL.ari_heatmap(pw, fig + "A3_ari_heatmap.jpeg", ds)
        PL.jaccard_heatmap(IO.load(ctx["out"] + "jaccard_cyto_vs_nuc.csv").set_index("a"), fig + "A3_jaccard_cyto_vs_nuc.jpeg", ds,
                           "nuclear cluster", "cytoplasmic cluster")
        xy = ctx["xy_sel"]; r = C.PRIMARY_RES
        for m in ("nuc_m", "cyto_m", "total_full"):
            PL.spatial_map(xy, ctx["clusters"][f"{m}_r{r}"].values, fig + f"A3_spatial_{m}.jpeg", categorical=True,
                           title=f"{C.short(ds)}: {m} partition (res {r})")
        co = ctx["cyto_only"]
        novel = co.loc[co["novel_vs_nuc"], "cyto_cluster"].tolist()
        if novel:
            PL.spatial_map(xy, ctx["clusters"][f"cyto_m_r{r}"].values, fig + "A3_spatial_cyto_only.jpeg", categorical=True,
                           highlight=novel, title=f"{C.short(ds)}: cytoplasm-specific clusters {novel}")
    # A3b: cytoplasm beyond nucleus
    if "resid_diag" in ctx:
        dg, xy = ctx["resid_diag"], ctx["xy_sel"]
        sc_ = ctx["resid_scores"]
        v = sc_["resid_pc1_h1"].values; lim = np.nanpercentile(np.abs(v), 99)
        PL.spatial_map(xy, np.clip(v, -lim, lim), fig + "A3b_map_resid_pc1.jpeg", cmap="RdBu_r", vmin=-lim, vmax=lim,
                       title=f"{C.short(ds)}: cytoplasmic residual PC1 ({ctx['resid_moran']['top_gene'].iloc[0]}), verdict: {dg['verdict']}")
        lab = ctx["resid_clusters"][f"resid_cyto_h1_r{C.PRIMARY_RES}"].values
        PL.spatial_map(xy, lab, fig + "A3b_map_resid_clusters.jpeg", categorical=True, title=f"{C.short(ds)}: cytoplasmic residual clusters")
        PL.residual_dotplot(ctx["resid_summary"], fig + "A3b_dotplot.jpeg", ds,
                            floor=ctx["cluster_spatial_strat"].query("partition == 'resid_floor_nuc_h2'"))
        PL.dose_response_plot(ctx["resid_moran"], fig + "A3b_dose_response.jpeg", ds)
        if dg.get("chosen_cluster") is not None:
            allobs = IO.load_all_obs(ds)
            PL.highlight_map(allobs[["global_x", "global_y"]].values, allobs["cell_type_merged"].values, xy,
                             lab == dg["chosen_cluster"], fig + "A3b_highlight.jpeg",
                             title=f"{C.short(ds)}: residual cluster {dg['chosen_cluster']} ({dg['chosen_label']})", zoom=dg.get("zoom_window"))
    # A4
    if "de_nucpart" in ctx:
        PL.paired_logfc(ctx["de_nucpart"], fig + "A4_paired_logfc_nucpart.jpeg", ds, "nuclear-partition")
        PL.paired_logfc(ctx["de_cytopart"], fig + "A4_paired_logfc_cytopart.jpeg", ds, "cytoplasmic-partition")
        ov = ctx["de_overlap"]
        PL.topk_overlap_bars(ov[ov["partition"] == "nucpart"], fig + "A4_topk_overlap.jpeg", ds)
    # C7
    if "assoc_wide" in ctx:
        axes = [a for a in C.CANDIDATE_AXES + ["subtype_cyto", "crowding_bin", "dist_nontumor_bin", "morph_area_bin", "seg_method"]]
        PL.axis_effect_heatmap(ctx["assoc_wide"], fig + "C7_axis_effects_frac.jpeg", ds, axes, model="frac")
        PL.axis_effect_heatmap(ctx["assoc_wide"], fig + "C7_axis_effects_cyto_minus_nuc.jpeg", ds, axes, model="cyto_minus_nuc")
        cand = ctx["candidates"]
        top = cand[cand["candidate"]].head(4)["gene"].tolist()
        pos = {g: i for i, g in enumerate(ctx["genes"])}
        for gname in top:
            j = pos[gname]
            k = np.asarray(ctx["Xc"][ctx["sel"]][:, j].todense()).ravel(); n = k + np.asarray(ctx["Xn"][ctx["sel"]][:, j].todense()).ravel()
            frac = np.where(n >= 3, (k + 0.5) / (n + 1), np.nan)
            PL.spatial_map(ctx["xy_sel"], frac, fig + f"C7_spatial_cytofrac_{gname}.jpeg", title=f"{C.short(ds)}: {gname} cytoplasmic fraction (n>=3)", vmin=0, vmax=1)
    # C8
    PL.pathway_summary(ctx["pathways_compartment"], fig + "C8_pathway_reliability.jpeg", ds, col="rel_param")
    if "pathways_assoc" in ctx:
        pa = ctx["pathways_assoc"]; sub = pa[(pa["model"] == "frac") & (pa["axis"].isin(C.CANDIDATE_AXES))]
        piv = sub.pivot_table(index="pathway", columns="axis", values="dev_expl_excess")
        piv = piv.loc[piv.max(axis=1).sort_values(ascending=False).head(25).index]
        import seaborn as sns
        fig_, ax = plt.subplots(figsize=(5, 7)); sns.heatmap(piv, cmap="Reds", ax=ax, yticklabels=[p[:45] for p in piv.index]); ax.tick_params(labelsize=6)
        ax.set_title(f"{C.short(ds)}: pathway cytoplasmic fraction, deviance explained beyond permutation", fontsize=8)
        PL.savefig(fig + "C8_pathway_axis_effects.jpeg")

if FORCE["figures"]:
    for ds, ctx in CTX.items():
        figures(ctx); print("figures written:", ctx["fig"])""")

md("""## Inspect one dataset

Quick looks at the per-dataset tables (change `ds`).""")

code("""ds = RUN_DATASETS[0]; ctx = CTX[ds]
print(json.dumps({k: v for k, v in ctx["qc"].items() if not isinstance(v, (list, dict))}, indent=1))
print(json.dumps(ctx["annot_provenance"], indent=1)); print(json.dumps({k: v for k, v in ctx["annot_flags"].items() if k != "marker_sets"}, indent=1))
display(ctx["annot_classes"][["cell_type_merged", "n_cells"] + [c for c in ctx["annot_classes"].columns if c.endswith("_mean")]].round(2))
if len(ctx.get("annot_clusters", [])): display(ctx["annot_clusters"].round(2))
display(ctx["annot_transfer"].round(3))
g = ctx["genes_compartment"]
display(g[g["is_pos_ctrl_nuc"] | g["is_pos_ctrl_cyto"]][["gene", "mean_total", "cyto_frac_pooled", "beta", "ci_lo", "ci_hi", "class", "rel_param", "beta_seg_range"]].round(3))
display(g.sort_values("rel_param", ascending=False).head(15)[["gene", "mean_total", "beta", "class", "rel_param", "rel_param_pct", "beta_seg_range_pct"]].round(3))
display(ctx["conc_tile"].dropna(subset=["r_true"]).nsmallest(15, "r_true")[["gene", "mean_per_unit", "r_nc", "rel_n", "rel_c", "r_true", "r_true_low_pct"]].round(3))
display(ctx["control_sets"][["set", "n_genes", "log_mean_expr", "beta", "rel_param"]].round(3))""")

code("""if "clusters" in ctx:
    pw = ctx["cluster_pairwise"]; display(pw[pw["resolution"] == C.PRIMARY_RES].round(3))
    display(ctx["cyto_only"].round(3)); display(ctx["recovery"].round(3))
    display(ctx["cluster_spatial"].round(3))
    cf = ctx.get("cluster_features")
    if cf is not None:
        display(cf[cf["matrix"] == "cyto_m"].pivot_table(index="cluster", columns="feature", values="z").round(1))
if "resid_diag" in ctx:
    print(json.dumps({k: v for k, v in ctx["resid_diag"].items() if k != "nuclear_pc_reliability"}, indent=1, default=str))
    display(ctx["resid_moran"][["pc", "top_gene", "reliability", "I", "I_perm_p99", "I_floor", "I_noise_control", "dr_tumor_frac_h1", "dr_tumor_frac_h2", "dr_immune_kernel_h1", "passes_A"]].round(3))
    display(ctx["resid_summary"][["cluster", "size", "purity_excess", "d_tumor_frac", "d_tumor_frac_lo", "d_tumor_frac_hi", "d_immune_kernel", "d_log_cell_area", "depth_ratio", "repro_f1", "n_cyto_only_de", "eligible", "morphology_driven", "rank", "label"]].round(3))
    st = ctx["cluster_spatial_strat"]; display(st[st["partition"] == "cyto_m"][["cluster", "size", "purity", "purity_excess", "d_tumor_frac", "d_immune_kernel", "depth_ratio"]].round(3))
if "de_nucpart" in ctx:
    d = ctx["de_nucpart"]; display(d[d["cyto_only"]].sort_values("delta", key=np.abs, ascending=False).head(20).round(2))
    display(ctx["de_overlap"])
if "assoc_wide" in ctx:
    w = ctx["assoc_wide"].set_index("gene")
    cols = [c for c in w.columns if c.endswith("__frac") and "_x_" not in c]
    display(w[cols].max(axis=1).sort_values(ascending=False).head(15).to_frame("max dev_expl_excess (frac)").join(w[cols]).round(4))
    beyond = [c for c in w.columns if "beyond_subtype__frac" in c]
    display(w[beyond].max(axis=1).sort_values(ascending=False).head(15).to_frame("max beyond subtype (frac)").join(w[beyond]).round(4))
display(ctx["candidates"].head(30).round(3))""")

md("""## Cross-sample section (B6, C9)

Recurrence of per-gene patterns and cytoplasmic clusters across the six samples, agreement of localisation between
samples, and the final candidate table.""")

code("""XO = IO.out_dir("cross_sample" + SUFFIX); XF = IO.out_dir("cross_sample" + SUFFIX, "fig")
# ---- annotation provenance and QC across samples (step 0b)
prov = XS.provenance_table(RUN_DATASETS, SUFFIX); IO.save(XO + "annotation_provenance.csv", prov); display(prov.round(3))
if IO.exists(IO.out_dir("cross_sample") + "annotation_transfer_lodo.parquet"):
    lodo = IO.load(IO.out_dir("cross_sample") + "annotation_transfer_lodo.parquet")
    print("leave-one-dataset-out accuracy of the non-epithelial classifier within BC/OC/CC:"); display(lodo.pivot_table(index="class", columns="held_out", values="accuracy").round(2))
tr = XS.merge_tables(RUN_DATASETS, "annotation_transfer.parquet", SUFFIX)
if len(tr):
    print("label transfer onto the manually annotated datasets (agreement for non-epithelial classes; confident non-epithelial fraction for tumor/unresolved classes):")
    display(tr[tr["dataset"].isin([C.short(d) for d in C.MANUAL_ANNOTATION])][["dataset", "label", "n", "agreement", "frac_confident", "top_prediction", "top_prediction_frac", "frac_confident_immune", "frac_confident_stromal"]].round(3))
cand_all = XS.merge_tables(RUN_DATASETS, "candidates.parquet", SUFFIX)
rec = XS.recurrence(cand_all); IO.save(XO + "gene_recurrence.parquet", rec)
IO.save(XO + "candidates_final.csv", rec[rec["n_candidate"] >= 1])
print("genes candidate in >= 1 / 2 / 3 samples:", [(rec["n_candidate"] >= k).sum() for k in (1, 2, 3)])
display(rec.head(40))
PL.recurrence_bar(rec[rec["n_candidate"] >= 1], XF + "B6_recurrence.jpeg", col="n_candidate", title="candidate genes: samples per gene")

beta = XS.gene_by_sample(cand_all, "beta"); IO.save(XO + "beta_by_sample.csv", beta.reset_index())
ok = beta.notna().sum(axis=1) == beta.shape[1]
PL.pairwise_scatter(beta[ok], XF + "B6_logor_pairwise.jpeg", "per-gene log odds ratio (cytoplasm vs offset)", lim=(-4, 3))
rt = XS.gene_by_sample(cand_all, "r_true_tile"); IO.save(XO + "r_true_tile_by_sample.csv", rt.reset_index())
PL.pairwise_scatter(rt.dropna(), XF + "B6_rtrue_pairwise.jpeg", "disattenuated nuclear-cytoplasmic r (tiles)", lim=(0, 1.2))
print("between-sample correlation of beta:"); display(beta.corr().round(2))

ax_rec = XS.residual_axis_recurrence(RUN_DATASETS, SUFFIX)
if len(ax_rec):
    IO.save(XO + "residual_axis_recurrence.csv", ax_rec); print("residual PC1 loading correlation between samples:"); display(ax_rec.round(2))
cl_rec = XS.cluster_recurrence(RUN_DATASETS, SUFFIX)
if len(cl_rec):
    IO.save(XO + "cluster_recurrence.parquet", cl_rec)
    display(cl_rec.groupby(["sample", "cluster"])["matched"].sum().unstack(fill_value=0).T if "matched" in cl_rec else cl_rec)""")

md("""## Decision summary

The numbers the target-form, identity-source and target-set decisions rest on. Interpretation goes into
`plans/fig1_exploration.md` after the full run.""")

code("""summary = XS.decision_summary(RUN_DATASETS, SUFFIX)
IO.save(XO + "decision_summary.csv", summary)
display(summary.T)
print(\"\"\"
Reading guide
- nuclear_share, median_d_match: depth available for matched comparisons (CC is the shallow sample).
- SG_rel vs ctrl_rel: localisation reliability of the SG set against an expression-matched control set.
- ari_nuc_cyto vs ari_cross_halves / ari_nuc_halves / ari_cyto_halves: the compartments carry different partitions
  only if cross-compartment agreement is well below the within-compartment half-split ceiling.
- ari_nuc_total, total_recovered_by_nuc_minF1: whether nuclear-only clustering recovers the subtypes seen in total
  expression (identity-source decision).
- n_cyto_only_de, n_sens_*, n_niche_beyond_subtype, n_candidates: size of the cytoplasm-specific signal per sample.
- resid_verdict, resid_dims_above_floor, resid_pc1_I, resid_chosen_cluster: the "cytoplasm beyond nucleus" panel rule
  (step 3b); "panel" means a reproducible, spatially coherent residual cluster with a boundary/immune effect exists.
\"\"\")""")

nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"name": "python3", "display_name": "Python 3 (AIVC-env)", "language": "python"},
                  "language_info": {"name": "python"}}
nbf.write(nb, "fig1_exploration.ipynb")
print("written fig1_exploration.ipynb with", len(cells), "cells")
