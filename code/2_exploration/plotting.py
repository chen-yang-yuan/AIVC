"""Figure helpers following the repo conventions (axes stripped for spatial maps, 300 dpi, jpeg/png)."""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import seaborn as sns

import config as C

color_cts = clr.LinearSegmentedColormap.from_list(
    "magma", ["#000003", "#3B0F6F", "#8C2980", "#F66E5B", "#FD9F6C", "#FBFCBF"], N=256)
CLASS_COLORS = {"nuclear-retained": "#3B4CC0", "cytoplasm-enriched": "#B40426", "balanced": "#999999",
                "ambiguous": "#DDAA33", "low-coverage": "#DDDDDD"}
DPI = 300


def savefig(path):
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close()


def slide_figsize(xy, short_edge=6.0):
    ext = np.nanmax(xy, 0) - np.nanmin(xy, 0)
    scale = short_edge / max(ext.min(), 1.0)
    return (max(ext[0] * scale, 2.0), max(ext[1] * scale, 2.0))


def spatial_map(xy, values, path, title="", categorical=False, s=0.5, cmap=None, vmin=None, vmax=None,
                highlight=None, legend=True):
    """Scatter of tumor cells in slide coordinates coloured by a value or a categorical label."""
    fig, ax = plt.subplots(figsize=slide_figsize(xy))
    if categorical:
        values = pd.Series(values).astype(str).values
        cats = pd.unique(values)
        palette = sns.color_palette("tab20", n_colors=max(len(cats), 1))
        for i, c in enumerate(cats):
            m = values == c
            col = palette[i % 20]
            if highlight is not None and c not in [str(h) for h in highlight]:
                col = "#DDDDDD"
            ax.scatter(xy[m, 0], xy[m, 1], s=s, c=[col], label=f"{c} (n={m.sum()})", linewidths=0, rasterized=True)
        if legend:
            ax.legend(markerscale=12, fontsize=6, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    else:
        sca = ax.scatter(xy[:, 0], xy[:, 1], s=s, c=values, cmap=cmap or color_cts, vmin=vmin, vmax=vmax,
                         linewidths=0, rasterized=True)
        cb = plt.colorbar(sca, ax=ax, fraction=0.03, pad=0.01)
        cb.ax.tick_params(labelsize=6)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    savefig(path)


def read_fractions(cells, path, ds):
    """Per-cell nuclear fraction: histogram, by segmentation method, by nucleus/cell area ratio bin."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2))
    ax = axes[0]
    ax.hist(cells["nuc_frac"].dropna(), bins=60, color="#555555")
    ax.axvline(np.nanmedian(cells["nuc_frac"]), color="red", lw=1)
    ax.set_xlabel("nuclear fraction of in-cell reads"); ax.set_ylabel("tumor cells")
    ax.set_title(f"{C.short(ds)}: median {np.nanmedian(cells['nuc_frac']):.2f}", fontsize=9)
    ax = axes[1]
    order = cells["seg_method"].value_counts().index
    sns.boxplot(data=cells, x="seg_method", y="nuc_frac", order=order, ax=ax, fliersize=0, color="#BBBBBB")
    ax.set_xticklabels([str(o).replace("Segmented by ", "")[:22] for o in order], rotation=20, fontsize=7)
    ax.set_xlabel(""); ax.set_ylabel("nuclear fraction")
    ax = axes[2]
    sns.boxplot(data=cells, x="ratio_bin", y="nuc_frac", ax=ax, fliersize=0, color="#BBBBBB")
    ax.set_xlabel("nucleus/cell area ratio quintile"); ax.set_ylabel("nuclear fraction")
    plt.tight_layout()
    savefig(path)


def logor_vs_abundance(gdf, path, ds, label_genes=None):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    x = np.log10(gdf["mean_total"].clip(lower=1e-4))
    for cls, col in CLASS_COLORS.items():
        m = gdf["class"] == cls
        ax.scatter(x[m], gdf.loc[m, "beta"], s=5, c=col, label=f"{cls} ({m.sum()})", linewidths=0, alpha=0.7)
    ax.axhline(0, color="k", lw=0.5); ax.axhline(C.LOG2, color="k", lw=0.5, ls="--"); ax.axhline(-C.LOG2, color="k", lw=0.5, ls="--")
    if label_genes:
        sub = gdf[gdf["gene"].isin(label_genes)]
        for _, r in sub.iterrows():
            ax.annotate(r["gene"], (np.log10(max(r["mean_total"], 1e-4)), r["beta"]), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
            ax.scatter([np.log10(max(r["mean_total"], 1e-4))], [r["beta"]], s=18, facecolors="none", edgecolors="k", linewidths=0.6)
    ax.set_xlabel("log10 mean in-cell count per tumor cell"); ax.set_ylabel("log odds ratio (cytoplasm vs cell offset)")
    ax.set_title(f"{C.short(ds)}: per-gene localisation", fontsize=9)
    ax.legend(fontsize=6, frameon=False)
    savefig(path)


def positive_controls(gdf, path, ds, nuc_ctrl, cyto_ctrl):
    sub = gdf[gdf["gene"].isin(list(nuc_ctrl) + list(cyto_ctrl))].copy()
    sub["kind"] = np.where(sub["gene"].isin(nuc_ctrl), "nuclear-retained control", "cytoplasmic control")
    sub = sub.sort_values(["kind", "beta"])
    fig, ax = plt.subplots(figsize=(max(4, 0.35 * len(sub) + 1), 3.5))
    cols = np.where(sub["kind"].str.startswith("nuclear"), CLASS_COLORS["nuclear-retained"], CLASS_COLORS["cytoplasm-enriched"])
    ax.bar(np.arange(len(sub)), sub["beta"], yerr=1.96 * sub["se"], color=cols, capsize=2)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(np.arange(len(sub))); ax.set_xticklabels(sub["gene"], rotation=60, fontsize=7)
    ax.set_ylabel("log odds ratio"); ax.set_title(f"{C.short(ds)}: positive controls", fontsize=9)
    savefig(path)


def concordance_scatter(cell_df, tile_df, path, ds, highlight=None):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, df, name in zip(axes, [cell_df, tile_df], ["cell level", "tile level"]):
        if df is None or len(df) == 0:
            continue
        ok = df["r_true"].notna()
        sca = ax.scatter(np.log10(df.loc[ok, "mean_per_unit"].clip(lower=1e-3)), df.loc[ok, "r_true"], s=6,
                         c=np.sqrt(df.loc[ok, "rel_n"].clip(0, 1) * df.loc[ok, "rel_c"].clip(0, 1)), cmap="viridis",
                         vmin=0, vmax=1, linewidths=0)
        ax.axhline(1, color="k", lw=0.5, ls="--")
        if highlight:
            sub = df[df["gene"].isin(highlight) & ok]
            for _, r in sub.iterrows():
                ax.annotate(r["gene"], (np.log10(max(r["mean_per_unit"], 1e-3)), r["r_true"]), fontsize=6)
        ax.set_xlabel("log10 mean count per unit"); ax.set_ylabel("disattenuated nuclear-cytoplasmic r")
        ax.set_title(f"{C.short(ds)}: {name} ({ok.sum()} genes)", fontsize=9)
        plt.colorbar(sca, ax=ax, label="sqrt(rel_n * rel_c)")
    plt.tight_layout()
    savefig(path)


def ari_heatmap(pairwise, path, ds, names=None, metric="ari"):
    names = names or C.CLUSTER_MATRICES
    M = pd.DataFrame(np.nan, index=names, columns=names)
    for _, r in pairwise.iterrows():
        M.loc[r["a"], r["b"]] = r[metric]; M.loc[r["b"], r["a"]] = r[metric]
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(M.astype(float), annot=True, fmt=".2f", cmap="viridis", vmin=0, vmax=1, ax=ax, annot_kws={"size": 7})
    ax.set_title(f"{C.short(ds)}: {metric.upper()} between partitions", fontsize=9)
    savefig(path)


def jaccard_heatmap(J, path, ds, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(0.5 * J.shape[1] + 2, 0.4 * J.shape[0] + 1.5))
    sns.heatmap(J, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1, ax=ax, annot_kws={"size": 6})
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(f"{C.short(ds)}: cluster Jaccard", fontsize=9)
    savefig(path)


def paired_logfc(de, path, ds, partition_name):
    fig, ax = plt.subplots(figsize=(5, 5))
    col = np.where(de["cyto_only"], CLASS_COLORS["cytoplasm-enriched"],
                   np.where(de["nuc_only"], CLASS_COLORS["nuclear-retained"], "#BBBBBB"))
    ax.scatter(de["lfc_nuc"], de["lfc_cyto"], s=4, c=col, linewidths=0, alpha=0.6)
    lim = np.nanpercentile(np.abs(de[["lfc_nuc", "lfc_cyto"]].values), 99.5)
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2 FC, nuclear half"); ax.set_ylabel("log2 FC, cytoplasmic half")
    ax.set_title(f"{C.short(ds)}: DE across {partition_name} clusters "
                 f"(cyto-only {int(de['cyto_only'].sum())}, nuc-only {int(de['nuc_only'].sum())})", fontsize=8)
    savefig(path)


def topk_overlap_bars(ov, path, ds):
    fig, ax = plt.subplots(figsize=(max(3, 0.5 * len(ov) + 1), 3))
    x = np.arange(len(ov))
    ax.bar(x, ov["n_shared"], color="#777777", label="shared")
    ax.bar(x, ov["n_a_only"], bottom=ov["n_shared"], color=CLASS_COLORS["nuclear-retained"], label="nuclear only")
    ax.bar(x, ov["n_b_only"], bottom=ov["n_shared"] + ov["n_a_only"], color=CLASS_COLORS["cytoplasm-enriched"], label="cytoplasmic only")
    ax.set_xticks(x); ax.set_xticklabels(ov["cluster"]); ax.set_xlabel("cluster"); ax.set_ylabel(f"top-{C.TOPK_DE} markers")
    ax.legend(fontsize=6, frameon=False); ax.set_title(f"{C.short(ds)}: marker overlap", fontsize=9)
    savefig(path)


def axis_effect_heatmap(wide, path, ds, axes, n_top=40, model="cyto_minus_nuc"):
    cols = [f"{a}__{model}" for a in axes if f"{a}__{model}" in wide]
    if not cols:
        return
    sub = wide.set_index("gene")[cols]
    top = sub.abs().max(axis=1).sort_values(ascending=False).head(n_top).index
    fig, ax = plt.subplots(figsize=(0.6 * len(cols) + 2, 0.22 * len(top) + 1.5))
    v = np.nanmax(np.abs(sub.loc[top].values)) if len(top) else 1
    sns.heatmap(sub.loc[top], cmap="RdBu_r", center=0, vmin=-v, vmax=v, ax=ax, annot=False,
                xticklabels=[c.replace(f"__{model}", "") for c in cols], yticklabels=True)
    ax.tick_params(axis="y", labelsize=6); ax.tick_params(axis="x", labelsize=7, rotation=45)
    ax.set_title(f"{C.short(ds)}: deviance explained ({model})", fontsize=9)
    savefig(path)


def pathway_summary(pw, path, ds, n_top=25, col="rel_param"):
    sub = pw.dropna(subset=[col]).sort_values(col, ascending=False).head(n_top)
    fig, ax = plt.subplots(figsize=(6, 0.25 * len(sub) + 1.5))
    ax.barh(np.arange(len(sub)), sub[col], color="#555555")
    ax.set_yticks(np.arange(len(sub))); ax.set_yticklabels([p[:45] for p in sub["pathway"]], fontsize=6)
    ax.invert_yaxis(); ax.set_xlabel(col); ax.set_title(f"{C.short(ds)}: pathway-level {col}", fontsize=9)
    savefig(path)


def recurrence_bar(rec, path, col="n_samples", title=""):
    counts = rec[col].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.bar(counts.index.astype(str), counts.values, color="#555555")
    ax.set_xlabel("number of samples"); ax.set_ylabel("genes"); ax.set_title(title, fontsize=9)
    savefig(path)


def pairwise_scatter(mat, path, title, lim=None):
    """Lower-triangle scatter matrix of a gene x sample table."""
    cols = list(mat.columns)
    n = len(cols)
    if n < 2:
        return
    fig, axes = plt.subplots(n, n, figsize=(1.6 * n, 1.6 * n), squeeze=False)
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if j >= i:
                ax.axis("off"); continue
            ok = mat[[cols[i], cols[j]]].notna().all(axis=1)
            ax.scatter(mat.loc[ok, cols[j]], mat.loc[ok, cols[i]], s=2, c="#444444", linewidths=0, alpha=0.5)
            r = np.corrcoef(mat.loc[ok, cols[j]], mat.loc[ok, cols[i]])[0, 1] if ok.sum() > 2 else np.nan
            ax.set_title(f"r={r:.2f}", fontsize=6, pad=1)
            if lim:
                ax.set_xlim(lim); ax.set_ylim(lim)
            ax.tick_params(labelsize=5)
            if i == n - 1: ax.set_xlabel(cols[j], fontsize=7)
            if j == 0: ax.set_ylabel(cols[i], fontsize=7)
    fig.suptitle(title, fontsize=9)
    savefig(path)


# ------------------------------------------------------------------ step 3b figures
def highlight_map(xy_all, types_all, xy_sel, mask, path, title="", zoom=None, immune_types=None):
    """All cells of the slide in grey, immune cells in colour, tumor cells pale, the highlighted cells dark; with a
    zoom inset on the pre-registered window."""
    immune_types = immune_types or C.IMMUNE_TYPES
    types_all = np.asarray(types_all)
    is_imm = np.isin(types_all, immune_types)
    is_tum = types_all == C.TUMOR_TYPE
    fig, ax = plt.subplots(figsize=slide_figsize(xy_all, 7.0))

    def draw(a, s):
        a.scatter(xy_all[~is_tum & ~is_imm, 0], xy_all[~is_tum & ~is_imm, 1], s=s, c="#E6E6E6", linewidths=0, rasterized=True)
        a.scatter(xy_all[is_imm, 0], xy_all[is_imm, 1], s=s, c="#F4A261", linewidths=0, rasterized=True, label="immune cells")
        a.scatter(xy_sel[~mask, 0], xy_sel[~mask, 1], s=s, c="#B8C4D6", linewidths=0, rasterized=True, label="other tumor cells")
        a.scatter(xy_sel[mask, 0], xy_sel[mask, 1], s=s * 1.5, c="#1F3A93", linewidths=0, rasterized=True, label="highlighted cluster")

    draw(ax, 0.4)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.legend(markerscale=15, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title(title, fontsize=9)
    if zoom is not None:
        x0, y0, w = zoom
        ax.add_patch(matplotlib.patches.Rectangle((x0, y0), w, w, fill=False, ec="k", lw=0.8))
        ins = ax.inset_axes([1.02, 0.0, 0.5, 0.5])
        draw(ins, 6)
        ins.set_xlim(x0, x0 + w); ins.set_ylim(y0, y0 + w); ins.set_aspect("equal")
        ins.set_xticks([]); ins.set_yticks([]); ins.set_title(f"{int(w)} um window", fontsize=7)
    savefig(path)


def residual_dotplot(summary, path, ds, floor=None):
    """Per residual cluster: purity excess, d on the features with CI, reproducibility, depth ratio, DE count."""
    cols = [("purity_excess", None), ("d_tumor_frac", ("d_tumor_frac_lo", "d_tumor_frac_hi")),
            ("d_immune_kernel", ("d_immune_kernel_lo", "d_immune_kernel_hi")),
            ("d_log1p_dist_nontumor", ("d_log1p_dist_nontumor_lo", "d_log1p_dist_nontumor_hi")),
            ("repro_f1", None), ("depth_ratio", None), ("n_cyto_only_de", None)]
    s = summary.sort_values("cluster")
    n = len(s)
    fig, axes = plt.subplots(1, len(cols), figsize=(2.0 * len(cols), 0.3 * n + 1.5), sharey=True)
    y = np.arange(n)
    for ax, (c, ci) in zip(axes, cols):
        v = s[c].values.astype(float)
        col = np.where(s["eligible"].values, "#1F3A93", np.where(s["morphology_driven"].values, "#DDAA33", "#999999"))
        ax.scatter(v, y, s=np.clip(s["size"].values / s["size"].max() * 120, 10, 120), c=col, zorder=3)
        if ci:
            ax.hlines(y, s[ci[0]].values, s[ci[1]].values, color="#555555", lw=1, zorder=2)
        if floor is not None and c in floor:
            ax.axvline(floor[c].max(), color="#AAAAAA", ls="--", lw=0.8)
        if c.startswith("d_") or c == "purity_excess":
            ax.axvline(0, color="k", lw=0.5)
        ax.set_title(c.replace("d_", "d: "), fontsize=7); ax.tick_params(labelsize=6)
    axes[0].set_yticks(y); axes[0].set_yticklabels([f"{int(c)} ({lab})" + (" *" if r == 1 else "") for c, lab, r in zip(s["cluster"], s["label"], s["rank"].fillna(0))], fontsize=6)
    fig.suptitle(f"{C.short(ds)}: cytoplasmic residual clusters (blue eligible, yellow morphology-driven, dashed = floor max)", fontsize=8)
    plt.tight_layout()
    savefig(path)


def dose_response_plot(moran, path, ds, features=("tumor_frac", "immune_kernel", "log1p_dist_nontumor")):
    fig, axes = plt.subplots(1, len(features), figsize=(3.2 * len(features), 2.8), sharey=True)
    for ax, f in zip(np.atleast_1d(axes), features):
        for _, r in moran.iterrows():
            qs = [c for c in moran.columns if c.startswith(f"dr_{f}_q")]
            ax.plot(np.arange(1, len(qs) + 1), r[qs].values.astype(float), marker="o", ms=3,
                    label=f"PC{int(r['pc'])} ({r['top_gene']}), effect {r[f'dr_{f}_h1']:+.2f} [{r[f'dr_{f}_h1_lo']:+.2f}, {r[f'dr_{f}_h1_hi']:+.2f}]")
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel(f"{f} quintile"); ax.tick_params(labelsize=7)
        ax.legend(fontsize=5, frameon=False)
    np.atleast_1d(axes)[0].set_ylabel("residual score (SD, within strata)")
    fig.suptitle(f"{C.short(ds)}: cytoplasmic residual score vs neighbourhood, within nuclear subtype", fontsize=8)
    plt.tight_layout()
    savefig(path)


# ------------------------------------------------------------------ step 0b figures
def marker_score_heatmap(classes, path, ds):
    cols = [c for c in classes.columns if c.endswith("_mean")]
    M = classes.set_index("cell_type_merged")[cols]; M.columns = [c.replace("_mean", "") for c in cols]
    Z = (M - M.mean(0)) / (M.std(0) + 1e-9)
    fig, ax = plt.subplots(figsize=(0.6 * Z.shape[1] + 2.5, 0.35 * Z.shape[0] + 1.5))
    sns.heatmap(Z, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, annot=M.round(2), fmt="", annot_kws={"size": 6}, ax=ax)
    ax.set_title(f"{C.short(ds)}: marker-set scores per class (colour = z across classes, text = mean log1p CPM)", fontsize=8)
    ax.tick_params(labelsize=7)
    savefig(path)


def tumor_cluster_heatmap(clusters, path, ds):
    if clusters is None or len(clusters) == 0:
        return
    zc = [c for c in clusters.columns if c.startswith("z_")]
    Z = clusters.set_index("cluster")[zc]; Z.columns = [c[2:] for c in zc]
    labels = [f"{i} (n={n}){' benign-like' if b else ''}{' contaminated' if c else ''}"
              for i, n, b, c in zip(clusters["cluster"], clusters["n_cells"], clusters["benign_like"], clusters["contaminated"])]
    fig, ax = plt.subplots(figsize=(0.6 * Z.shape[1] + 2.5, 0.35 * Z.shape[0] + 1.5))
    sns.heatmap(Z, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, annot=True, fmt=".1f", annot_kws={"size": 6}, ax=ax, yticklabels=labels)
    ax.set_title(f"{C.short(ds)}: graph clusters called malignant, marker z-scores across clusters", fontsize=8)
    ax.tick_params(labelsize=7)
    savefig(path)
