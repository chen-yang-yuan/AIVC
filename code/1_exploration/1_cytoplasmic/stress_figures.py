"""Figure helpers of 6_stress_substates.ipynb: spatial maps, boxplots in physical units, stacked compositions,
dot plots, enrichment heatmaps, composition bars, forest plots, paired bars and Sankey-style alluvial flows between
labellings (`alluvial_figure`). Substate colours come from subtype_utils.LABEL_COLORS (fixed vocabulary shared by
every subtyping run); stress subtypes have their own fixed palette below."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # shared subtype_utils.py
from subtype_utils import LABEL_COLORS, LABEL_ORDER

SAMPLE_COLORS = {"BC": "#1f77b4", "OC": "#ff7f0e", "CC": "#2ca02c", "LC": "#d62728", "Prostate": "#9467bd"}
COMP_COLORS = {"cyto": "#d95f02", "nuc": "#1b9e77", "joint": "#7570b3"}

STRESS_COLORS = {
    "Unstressed": "#bdbdbd",
    "Vascular deprivation": "#1f78b4", "Crowding": "#6a3d9a", "Core isolation": "#a6cee3",
    "Immune attack": "#e31a1c", "Myeloid contact": "#fb9a99", "Stromal contact": "#b15928",
    "Hypoxia (expr)": "#33a02c", "Heat shock": "#ff7f00", "Inflammatory (TNFA/NFKB)": "#fdbf6f",
    "Interferon (IFNG)": "#cab2d6",
}
_EXTRA = ["#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"]


def stress_color(name):
    """Colour of a stress subtype (fixed palette; unknown names get a fallback colour)."""
    if name in STRESS_COLORS:
        return STRESS_COLORS[name]
    return _EXTRA[abs(hash(name)) % len(_EXTRA)]


def label_palette(categories):
    """Colour of every substate label (fixed) or stress subtype."""
    return {c: LABEL_COLORS.get(c, stress_color(c)) for c in categories}


def ordered(categories):
    cats = list(pd.unique(np.asarray(categories, dtype=object)))
    known = [c for c in LABEL_ORDER if c in cats]
    rest = sorted(c for c in cats if c not in LABEL_ORDER and c != "Unstressed")
    return known + rest + (["Unstressed"] if "Unstressed" in cats else [])


# ----------------------------------------------------------------------------------------------------------------
# Spatial panels
# ----------------------------------------------------------------------------------------------------------------
def _clean(ax, title):
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=8)


def categorical_panel(ax, xy, labels, title="", palette=None, s=0.3, legend=True, highlight=None, order=None):
    """Scatter coloured by category. With `highlight`, only that category is coloured (others light grey)."""
    labels = np.asarray(labels, dtype=object)
    cats = order or ordered(labels)
    pal = palette or label_palette(cats)
    if highlight is not None:
        ax.scatter(xy[:, 0], xy[:, 1], s=s, c="#e0e0e0", linewidths=0, rasterized=True)
        m = labels == highlight
        ax.scatter(xy[m, 0], xy[m, 1], s=max(s, 0.6), c=pal.get(highlight, "#d62728"), linewidths=0, rasterized=True)
        _clean(ax, f"{title}\n{highlight}: {m.sum():,} of {len(m):,} cells")
        return
    colors = np.array([pal[c] for c in cats], dtype=object)
    codes = pd.Categorical(labels, categories=cats).codes
    o = np.random.default_rng(0).permutation(len(codes))
    ax.scatter(xy[o, 0], xy[o, 1], s=s, c=colors[codes[o]], linewidths=0, rasterized=True)
    _clean(ax, title)
    if legend:
        counts = pd.Series(labels).value_counts()
        handles = [Line2D([], [], marker="o", ls="", ms=4, color=pal[c], label=f"{c} ({counts.get(c, 0):,})") for c in cats]
        ax.legend(handles=handles, fontsize=5, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)


def continuous_panel(ax, xy, values, title="", cmap="magma", s=0.3, pct=(2, 98), label=""):
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    lo, hi = np.nanpercentile(v, pct[0]), np.nanpercentile(v, pct[1])
    o = np.argsort(np.where(ok, v, -np.inf))
    sc_ = ax.scatter(xy[o, 0], xy[o, 1], s=s, c=np.clip(np.where(ok, v, lo)[o], lo, hi), cmap=cmap, vmin=lo, vmax=hi, linewidths=0, rasterized=True)
    _clean(ax, title)
    cb = plt.colorbar(sc_, ax=ax, fraction=0.03, pad=0.01)
    cb.ax.tick_params(labelsize=6)
    if label:
        cb.set_label(label, fontsize=6)


def spatial_pair_figure(stress, labels, sample_levels, sample_codes, pair, axis_title, path, s=0.3):
    """2 x n_samples: substate highlighted over all tumor cells (top) and the stress axis as a map (bottom)."""
    sub, axname = pair
    n = len(sample_levels)
    fig, axes = plt.subplots(2, n, figsize=(3.6 * n, 7.2))
    for si, sname in enumerate(sample_levels):
        m = sample_codes == si
        xy = stress.loc[m, ["global_x", "global_y"]].values
        categorical_panel(axes[0, si], xy, labels[m], title=sname, highlight=sub, s=s)
        continuous_panel(axes[1, si], xy, stress.loc[m, axname].values, title=f"{sname}: {axis_title}", s=s)
    fig.suptitle(f"{sub} substate (top) vs {axis_title} (bottom)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def spatial_labels_figure(stress, columns, sample_mask, title, path, s=0.3, palettes=None):
    """One row of categorical panels (columns = dict name -> label array over all cells) for one sample."""
    xy = stress.loc[sample_mask, ["global_x", "global_y"]].values
    fig, axes = plt.subplots(1, len(columns), figsize=(5.2 * len(columns), 5.5))
    axes = np.atleast_1d(axes)
    for ax, (name, lab) in zip(axes, columns.items()):
        categorical_panel(ax, xy, np.asarray(lab, dtype=object)[sample_mask], title=name, s=s,
                          palette=(palettes or {}).get(name))
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def spatial_axes_figure(stress, axes_titles, sample_mask, title, path, s=0.3, ncol=3):
    names = list(axes_titles)
    nrow = int(np.ceil(len(names) / ncol))
    xy = stress.loc[sample_mask, ["global_x", "global_y"]].values
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.3 * nrow))
    for k, ax in enumerate(np.atleast_1d(axes).ravel()):
        if k >= len(names):
            ax.axis("off"); continue
        continuous_panel(ax, xy, stress.loc[sample_mask, names[k]].values, title=axes_titles[names[k]], s=s)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def _box_stats(v):
    q1, med, q3 = np.nanpercentile(v, [25, 50, 75])
    return {"med": med, "q1": q1, "q3": q3, "whislo": np.nanpercentile(v, 5), "whishi": np.nanpercentile(v, 95), "fliers": []}


def boxplot_axes_figure(stress, labels, sample_codes, sample_levels, axes_titles, path, label_order, ncol=5):
    """Pooled box (5-95 % whiskers) of each raw axis per substate, one dot per sample median."""
    names = list(axes_titles)
    nrow = int(np.ceil(len(names) / ncol))
    pal = label_palette(label_order)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.4 * nrow))
    labels = np.asarray(labels, dtype=object)
    for k, ax in enumerate(np.atleast_1d(axes).ravel()):
        if k >= len(names):
            ax.axis("off"); continue
        v = stress[names[k]].values.astype(float)
        stats_ = [_box_stats(v[labels == l]) for l in label_order]
        bp = ax.bxp(stats_, positions=np.arange(len(label_order)), widths=0.6, showfliers=False, patch_artist=True)
        for patch, l in zip(bp["boxes"], label_order):
            patch.set_facecolor(pal[l]); patch.set_alpha(0.7); patch.set_edgecolor("k"); patch.set_linewidth(0.5)
        for si, sname in enumerate(sample_levels):
            ms = sample_codes == si
            meds = [np.nanmedian(v[ms & (labels == l)]) if (ms & (labels == l)).sum() >= 50 else np.nan for l in label_order]
            ax.plot(np.arange(len(label_order)) + (si - 2) * 0.08, meds, ls="", marker="o", ms=2.5, color=SAMPLE_COLORS[sname], label=sname if k == 0 else None)
        ax.set_xticks(np.arange(len(label_order))); ax.set_xticklabels(label_order, rotation=90, fontsize=6)
        ax.set_title(axes_titles[names[k]], fontsize=8); ax.tick_params(axis="y", labelsize=6); ax.grid(alpha=0.3, axis="y")
    handles = [Line2D([], [], marker="o", ls="", color=SAMPLE_COLORS[s], label=s) for s in sample_levels]
    fig.legend(handles=handles, loc="upper right", fontsize=7, title="sample median", title_fontsize=7, frameon=False)
    fig.suptitle("Raw stress axes per substate (box: pooled 25-75 %, whiskers 5-95 %; dots: per-sample medians)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def stacked_groups_figure(stress, labels, axes_names, group_names, path, label_order, ncol=6):
    """Share of the Low / Mid / High (Away / Close) cells inside each substate, one panel per axis."""
    nrow = int(np.ceil(len(axes_names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.0 * nrow), sharey=True)
    labels = np.asarray(labels, dtype=object)
    cols3 = ["#4575b4", "#d9d9d9", "#d73027"]
    for k, ax in enumerate(np.atleast_1d(axes).ravel()):
        if k >= len(axes_names):
            ax.axis("off"); continue
        axname = axes_names[k]
        grp = stress[f"grp_{axname}"].values
        order = group_names[axname]
        share = pd.crosstab(pd.Categorical(labels, categories=label_order), pd.Categorical(grp, categories=order), normalize="index").reindex(label_order)
        bottom = np.zeros(len(label_order))
        cols = cols3 if len(order) == 3 else [cols3[0], cols3[2]]
        for g, c in zip(order, cols):
            ax.bar(np.arange(len(label_order)), share[g].values, bottom=bottom, color=c, width=0.8, label=g)
            bottom += share[g].values
        ax.axhline(1 / len(order), color="k", lw=0.5, ls="--")
        ax.set_xticks(np.arange(len(label_order))); ax.set_xticklabels(label_order, rotation=90, fontsize=5)
        ax.set_title(axname, fontsize=8); ax.legend(fontsize=5, loc="upper right", frameon=False); ax.tick_params(axis="y", labelsize=6)
    fig.suptitle("Share of each stress group inside every substate (dashed line: share expected under independence)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def dotplot_figure(med, share_high, path, title, vmax=1.0):
    """Rows = axes, columns = substates; colour = median z, size = share of the substate in the High group."""
    fig, ax = plt.subplots(figsize=(0.55 * med.shape[1] + 3, 0.28 * med.shape[0] + 2))
    rows, cols = list(med.index), list(med.columns)
    X, Y = np.meshgrid(np.arange(len(cols)), np.arange(len(rows)))
    sizes = 250 * np.nan_to_num(share_high.values) ** 1.0
    sc_ = ax.scatter(X.ravel(), Y.ravel(), s=sizes.ravel(), c=med.values.ravel(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, edgecolors="k", linewidths=0.3)
    ax.set_xticks(np.arange(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(rows))); ax.set_yticklabels(rows, fontsize=7); ax.invert_yaxis()
    ax.set_xlim(-0.6, len(cols) - 0.4)
    cb = fig.colorbar(sc_, ax=ax, fraction=0.03, pad=0.02); cb.set_label("median pooled z", fontsize=7); cb.ax.tick_params(labelsize=6)
    for frac in (0.2, 0.33, 0.5):
        ax.scatter([], [], s=250 * frac, c="w", edgecolors="k", label=f"{int(frac * 100)} % High")
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.12, 1.0), frameon=False, title="dot size", title_fontsize=6)
    ax.set_title(title, fontsize=9)
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------
# Label vs label
# ----------------------------------------------------------------------------------------------------------------
def enrichment_heatmaps(items, path, title, vmax=2.0):
    """items: list of (panel title, contingency DataFrame, long DataFrame with log2_enrichment). Colour = log2
    (observed / expected), text = cell count."""
    n = len(items)
    widths = [max(3.0, 0.5 * it[1].shape[1] + 1.5) for it in items]
    fig, axes = plt.subplots(1, n, figsize=(sum(widths) + 1, 0.42 * max(it[1].shape[0] for it in items) + 2.5), gridspec_kw={"width_ratios": widths, "wspace": 0.08})
    axes = np.atleast_1d(axes)
    for k, (ax, (t, cont, long)) in enumerate(zip(axes, items)):
        rows, cols = list(cont.index), list(cont.columns)
        E = long.pivot(index="substate", columns="stress_subtype", values="log2_enrichment").reindex(index=rows, columns=cols).values
        im = ax.imshow(E, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(cols)):
                ax.text(j, i, f"{int(cont.values[i, j]):,}", ha="center", va="center", fontsize=5, color="k")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontsize=6)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows if k == 0 else [], fontsize=6)
        ax.set_title(t, fontsize=8)
    fig.colorbar(im, ax=axes, shrink=0.6, label="log2 (observed / expected)")
    fig.suptitle(title, fontsize=10)
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def composition_figure(items, path, title):
    """items: list of (panel title, substate labels, stress labels). Two rows: stress subtypes within each
    substate (top), substates within each stress subtype (bottom)."""
    n = len(items)
    fig, axes = plt.subplots(2, n, figsize=(5.5 * n, 8))
    axes = np.atleast_2d(axes)
    for k, (t, a, b) in enumerate(items):
        a = np.asarray(a, dtype=object); b = np.asarray(b, dtype=object)
        for r, (x, y, xl) in enumerate([(a, b, "substate"), (b, a, "stress subtype")]):
            xo, yo = ordered(x), ordered(y)
            pal = label_palette(yo)
            share = pd.crosstab(pd.Categorical(x, categories=xo), pd.Categorical(y, categories=yo), normalize="index").reindex(xo)
            ax = axes[r, k]
            left = np.zeros(len(xo))
            for c in yo:
                ax.barh(np.arange(len(xo)), share[c].values, left=left, color=pal[c], height=0.8, label=c)
                left += share[c].values
            ax.set_yticks(np.arange(len(xo))); ax.set_yticklabels([f"{l} ({(x == l).sum():,})" for l in xo], fontsize=6); ax.invert_yaxis()
            ax.set_xlim(0, 1); ax.set_xlabel(f"share of the {xl}'s cells", fontsize=7); ax.tick_params(axis="x", labelsize=6)
            ax.set_title(f"{t}: {'stress subtypes within substates' if r == 0 else 'substates within stress subtypes'}", fontsize=8)
            ax.legend(handles=[Patch(color=pal[c], label=c) for c in yo], fontsize=5, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def forest_figure(tab, sample_levels, path, title, value="log2_or_pooled", lo="ci_low", hi="ci_high", per_sample_prefix="log2_or_"):
    """Pooled estimate with CI (black) and per-sample estimates (coloured dots) for the rows of `tab`."""
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(tab) + 1.5))
    y = np.arange(len(tab))
    ax.errorbar(tab[value], y, xerr=[tab[value] - tab[lo], tab[hi] - tab[value]], fmt="s", color="k", ms=4, capsize=2, lw=1, label="pooled (block-bootstrap 95 % CI)")
    for si, s in enumerate(sample_levels):
        ax.plot(tab[f"{per_sample_prefix}{s}"], y + (si - 2) * 0.12, ls="", marker="o", ms=3, color=SAMPLE_COLORS[s], label=s)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(tab["name"], fontsize=7); ax.invert_yaxis(); ax.grid(alpha=0.3, axis="x")
    ax.set_xlabel("log2 odds ratio (High / Close group inside vs outside the substate)", fontsize=8)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_title(title, fontsize=9)
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def heatmap(ax, M, row_labels, col_labels, title, vmax=None, cmap="RdBu_r", row_breaks=()):
    v = np.nanmax(np.abs(M)) if vmax is None else vmax
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=-v, vmax=v)
    ax.set_xticks(range(len(col_labels))); ax.set_xticklabels(col_labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(row_labels))); ax.set_yticklabels(row_labels, fontsize=7)
    for rb in row_breaks:
        ax.axhline(rb - 0.5, color="k", lw=0.8)
    ax.set_title(title, fontsize=10)
    return im


def paired_bars(ax, values_a, values_b, names, title, label_a="cytoplasmic", label_b="nuclear", values_c=None, label_c="both", xline=0.0, row_breaks=()):
    y = np.arange(len(names))
    w = 0.25 if values_c is not None else 0.4
    ax.barh(y - w, values_a, w, label=label_a, color=COMP_COLORS["cyto"])
    ax.barh(y if values_c is not None else y + w, values_b, w, label=label_b, color=COMP_COLORS["nuc"])
    if values_c is not None:
        ax.barh(y + w, values_c, w, label=label_c, color=COMP_COLORS["joint"])
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7); ax.invert_yaxis(); ax.set_title(title, fontsize=9); ax.grid(alpha=0.3)
    for rb in row_breaks:
        ax.axhline(rb - 0.5, color="k", lw=0.6)
    if xline is not None:
        ax.axvline(xline, color="grey", lw=0.8)


# ----------------------------------------------------------------------------------------------------------------
# Sankey / alluvial flows between labellings of the same cells
# ----------------------------------------------------------------------------------------------------------------
def alluvial_figure(columns, path, title, color_col=1, highlight=None, min_share=0.0, gap=0.015, width=0.12, figsize=None, fontsize=7):
    """Sankey-style figure: `columns` = list of (column title, labels array) for the same cells (2 or 3 columns).
    Node height = cell count; ribbon width between adjacent columns = cells carrying that pair of labels.
    Ribbons take the colour of the node in column `color_col`; those in `highlight` (set of (left column index,
    left label, right label)) are drawn opaque, the rest faint. `min_share` skips ribbons thinner than that
    share of the smaller of the two nodes."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch, Rectangle
    highlight = highlight or set()
    ncol = len(columns)
    labs = [np.asarray(c[1], dtype=object) for c in columns]
    n = len(labs[0])
    cats = [ordered(l) for l in labs]
    pal = [label_palette(c) for c in cats]
    counts = [pd.Series(l).value_counts() for l in labs]
    scale = [1.0 - gap * (len(c) - 1) for c in cats]
    pos = []
    for k in range(ncol):
        y, d = 1.0, {}
        for c in cats[k]:
            h = scale[k] * counts[k].get(c, 0) / n
            d[c] = (y - h, y); y -= h + gap
        pos.append(d)
    xs = np.linspace(0.0, 1.0, ncol)
    fig, ax = plt.subplots(figsize=figsize or (4.8 * ncol + 1.5, 9))
    for k in range(ncol - 1):
        ct = pd.crosstab(labs[k], labs[k + 1]).reindex(index=cats[k], columns=cats[k + 1], fill_value=0)
        off_l = {c: pos[k][c][1] for c in cats[k]}
        off_r = {c: pos[k + 1][c][1] for c in cats[k + 1]}
        x0, x1 = xs[k] + width / 2, xs[k + 1] - width / 2
        xm = (x0 + x1) / 2
        patches = []
        for a in cats[k]:
            for b in cats[k + 1]:
                v = ct.loc[a, b]
                if v == 0:
                    continue
                hl_, hr_ = scale[k] * v / n, scale[k + 1] * v / n
                yl1, yr1 = off_l[a], off_r[b]
                yl0, yr0 = yl1 - hl_, yr1 - hr_
                off_l[a] = yl0; off_r[b] = yr0
                if v < min_share * min(counts[k].get(a, 0), counts[k + 1].get(b, 0)):
                    continue
                verts = [(x0, yl1), (xm, yl1), (xm, yr1), (x1, yr1), (x1, yr0), (xm, yr0), (xm, yl0), (x0, yl0), (x0, yl1)]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
                col = pal[k][a] if color_col == k else (pal[k + 1][b] if color_col == k + 1 else pal[k][a])
                hi = (k, a, b) in highlight
                patches.append((hi, PathPatch(Path(verts, codes), facecolor=col, edgecolor="none", alpha=0.85 if hi else 0.12, zorder=2 if hi else 1)))
        for _, p in sorted(patches, key=lambda t: t[0]):
            ax.add_patch(p)
    for k in range(ncol):
        for c in cats[k]:
            y0, y1 = pos[k][c]
            ax.add_patch(Rectangle((xs[k] - width / 2, y0), width, y1 - y0, facecolor=pal[k][c], edgecolor="k", lw=0.4, zorder=3))
            txt = f"{c} ({counts[k].get(c, 0):,})"
            if k == 0:
                ax.text(xs[k] - width / 2 - 0.01, (y0 + y1) / 2, txt, ha="right", va="center", fontsize=fontsize)
            elif k == ncol - 1:
                ax.text(xs[k] + width / 2 + 0.01, (y0 + y1) / 2, txt, ha="left", va="center", fontsize=fontsize)
            elif y1 - y0 >= 0.012:
                ax.text(xs[k], (y0 + y1) / 2, txt, ha="center", va="center", fontsize=fontsize - 1, zorder=4,
                        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.0))
        ax.text(xs[k], 1.025, columns[k][0], ha="center", va="bottom", fontsize=fontsize + 2, fontweight="bold")
    ax.set_xlim(-0.42, 1.42); ax.set_ylim(-0.01, 1.07); ax.axis("off")
    ax.set_title(title, fontsize=9, pad=18)
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)
