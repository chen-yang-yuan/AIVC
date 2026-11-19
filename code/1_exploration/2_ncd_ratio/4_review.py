"""NCD ratio analysis, step 4, review: compares the 24 embedding tasks of 4_embedding.py and audits which gene
list and which per-sample correction to carry on (plans/ncd_ratio_exploration.md, section 5.4).

Light and local: reads only the small files of output/1_exploration/ncd_ratio/4_embedding/<list>/<correction>/<run>/
(cells.csv.gz, pca_variance.csv, checks.csv, cluster_gene_means.csv.gz), pulled from HGCC with `make ncd-pull`.

One row per (list, correction) in comparison.csv, at the reference resolution REF_RES; by_resolution.csv has the
resolution-dependent numbers for every resolution of RES_GRID.

  Criterion              Measure
  sample effect          share of a cell's neighbours from its own sample divided by the sample's share of the cells
                         (mean over samples; 1 = perfect mixing); Cramer's V of cluster and sample
  depth and geometry     largest |Spearman| of PC1 to PC10 with log depth and with the nuclear share of counts; share
                         of the variance of log depth and of the nuclear share between clusters (eta squared)
  more than counting     variance explained by the first N_PCS PCs, data and null; number of leading PCs whose
                         variance exceeds the PC of the same rank of the null by more than 10 %; clusters of the null
  repeatable             adjusted Rand index (ARI) of the clusters of the two transcript halves (same cells,
                         independent transcripts); ARI of each half with the data
  agreement of lists     ARI of the clusters of the data between the three lists, same correction (list_agreement.csv)
  spatial pattern        share of a cell's K_SPACE nearest tumor cells in the tissue that are in its cluster, divided
                         by the share expected from the cluster sizes in the sample (mean over samples); the same for
                         the null
  gene content           marker genes of a cluster = its N_MARKERS genes with the highest mean corrected score;
                         clusters with a curated program or hallmark set enriched among the markers (hypergeometric
                         against the genes of the list, p < P_PROGRAM); surface and secreted genes among all markers
                         (annotation of 3_audit/gene_audit.csv)

Figures: umap_<list>_<correction>.png (data and null; by sample, cluster, depth, nuclear share), overview.png.

Run from this folder (AIVC-env, minutes):  python 4_review.py [--quick]
Output: output/1_exploration/ncd_ratio/4_embedding/review/ (emptied first; --quick reads and writes 4_embedding_quick/)
"""
import argparse
import os
import shutil
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import hypergeom
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # shared subtype_utils.py
import ncd_ratio_utils as nr
import subtype_utils as su

REF_RES = 1.0
K_SPACE = 10
N_MARKERS = 30
P_PROGRAM = 1e-3
N_PLOT = 150000   # cells drawn in a UMAP panel
LIST_NAMES = {"primary": "V-abs 500", "larger": "V-abs 1000", "union": "union 630"}


def load(root, lst, corr, run):
    d = f"{root}{lst}/{corr}/{run}/"
    assert os.path.exists(d + "summary.txt"), f"task {lst}/{corr}/{run} is not complete in {root}"
    return pd.read_csv(d + "cells.csv.gz"), pd.read_csv(d + "pca_variance.csv"), pd.read_csv(d + "checks.csv"), d


def cramers_v(a, b):
    t = pd.crosstab(a, b).values.astype(float)
    n = t.sum()
    e = t.sum(axis=1, keepdims=True) * t.sum(axis=0, keepdims=True) / n
    chi2 = ((t - e) ** 2 / np.maximum(e, 1e-12)).sum()
    return float(np.sqrt(chi2 / (n * max(min(t.shape) - 1, 1))))


def eta2(values, labels):
    """Share of the variance of values that lies between the groups of labels."""
    v = pd.Series(values).astype(float)
    ok = v.notna().values
    v, g = v[ok], pd.Series(labels)[ok].values
    m = v.groupby(g).agg(["mean", "size"])
    return float((m["size"] * (m["mean"] - v.mean()) ** 2).sum() / ((v - v.mean()) ** 2).sum())


def spatial_ratio(cells, labels):
    """Observed / expected share of the K_SPACE nearest tumor cells of the same sample that share the cell's cluster."""
    out = []
    for s, idx in cells.groupby("sample").indices.items():
        xy = cells[["global_x", "global_y"]].values[idx]
        lab = np.asarray(labels)[idx]
        _, nb = cKDTree(xy).query(xy, k=K_SPACE + 1)
        obs = (lab[nb[:, 1:]] == lab[:, None]).mean()
        p = np.bincount(lab) / len(lab)
        out.append(obs / (p ** 2).sum())
    return float(np.mean(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="read and write 4_embedding_quick/")
    args = ap.parse_args()
    t0 = time.time()
    root = nr.step_dir("4_embedding", args.quick)
    out = root + "review/"
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)

    groups = {f"curated: {k}": v for k, v in su.CURATED_PROGRAMS.items()}
    groups.update({f"hallmark: {k.replace('HALLMARK_', '').replace('_', ' ').capitalize()}": v
                   for k, v in su.read_gmt(nr.UTILS_PATH + "hallmark_pathways_filtered.gmt").items()})
    ga = pd.read_csv(nr.step_dir("3_audit") + "gene_audit.csv").set_index("gene")

    rows, rows_res, labels_ref, marker_rows = [], [], {}, []
    for lst in nr.LISTS:
        genes = nr.load_target_genes(lst)
        for corr in nr.CORRECTIONS:
            cd, vd, kd, ddir = load(root, lst, corr, "data")
            cn, vn, kn, _ = load(root, lst, corr, "null")
            ca, _, _, _ = load(root, lst, corr, "half_a")
            cb, _, _, _ = load(root, lst, corr, "half_b")
            for other in (cn, ca, cb):
                assert (other["cell_id"].values == cd["cell_id"].values).all() and (other["sample"].values == cd["sample"].values).all()
            key = f"cluster_r{REF_RES}"
            labels_ref[lst, corr] = cd[key].values
            npc = nr.N_PCS
            sp = lambda k, name: k[(k["check"] == "spearman_pc_covariate") & (k["sample"] == "pooled") & (k["covariate"] == name)].set_index("pc")["value"]
            nn = lambda k: k[k["check"].isin(["same_sample_neighbour_share", "sample_share_of_cells"])].pivot(index="sample", columns="check", values="value")
            mix_d, mix_n = nn(kd), nn(kn)
            depth_rho, lean_rho = sp(kd, "log_depth"), sp(kd, "nuclear_share")
            vdr, vnr = vd["variance_ratio"].values, vn["variance_ratio"].values
            m = min(len(vdr), len(vnr), npc)
            above = vdr[:m] > 1.1 * vnr[:m]
            lead = int(np.argmin(above)) if not above.all() else m   # leading PCs above the null

            # marker genes and their content, reference resolution
            gm = pd.read_csv(ddir + "cluster_gene_means.csv.gz")
            gm = gm[np.isclose(gm["resolution"], REF_RES)]
            all_markers, enriched = set(), 0
            for k, g in gm.groupby("cluster"):
                top = g.sort_values("mean_corrected", ascending=False).head(N_MARKERS)["gene"].values
                all_markers.update(top)
                best = ("", 1.0, 0)
                for name, gg in groups.items():
                    inlist = np.isin(genes, gg).sum()
                    hit = np.isin(top, gg).sum()
                    if inlist >= 5 and hit >= 3:
                        p = hypergeom.sf(hit - 1, len(genes), inlist, len(top))
                        if p < best[1]:
                            best = (name, p, hit)
                enriched += best[1] < P_PROGRAM
                marker_rows.append({"list": lst, "correction": corr, "cluster": k, "cells": int((cd[key] == k).sum()),
                                    "mean_r_of_cluster": float(g["mean_r"].mean()), "best_program": best[0], "p": best[1], "genes_hit": best[2],
                                    "markers": " ".join(top[:N_MARKERS])})
            am = [g for g in all_markers if g in ga.index]
            rows.append({
                "list": lst, "correction": corr, "genes": len(genes), "cells": len(cd),
                "neighbours_from_own_sample_over_chance": float((mix_d["same_sample_neighbour_share"] / mix_d["sample_share_of_cells"]).mean()),
                "same_for_the_null": float((mix_n["same_sample_neighbour_share"] / mix_n["sample_share_of_cells"]).mean()),
                "cramers_v_cluster_sample": cramers_v(cd[key], cd["sample"]),
                "max_abs_rho_pc_depth": float(depth_rho.abs().max()), "pc_of_max_depth": int(depth_rho.abs().idxmax()),
                "max_abs_rho_pc_nuclear_share": float(lean_rho.abs().max()), "pc_of_max_nuclear_share": int(lean_rho.abs().idxmax()),
                "eta2_depth_by_cluster": eta2(cd["log_depth"], cd[key]), "eta2_nuclear_share_by_cluster": eta2(cd["nuclear_share"], cd[key]),
                f"variance_first_{npc}_pcs_data": float(vdr[:npc].sum()), f"variance_first_{npc}_pcs_null": float(vnr[:npc].sum()),
                "variance_pc1_data": float(vdr[0]), "variance_pc1_null": float(vnr[0]), "leading_pcs_above_null": lead,
                "clusters_data": int(cd[key].nunique()), "clusters_null": int(cn[key].nunique()),
                "ari_half_a_half_b": adjusted_rand_score(ca[key], cb[key]),
                "ari_half_with_data": float(np.mean([adjusted_rand_score(ca[key], cd[key]), adjusted_rand_score(cb[key], cd[key])])),
                "ari_null_with_data": adjusted_rand_score(cn[key], cd[key]),
                "spatial_same_cluster_over_chance_data": spatial_ratio(cd, cd[key].values),
                "spatial_same_cluster_over_chance_null": spatial_ratio(cn, cn[key].values),
                "clusters_with_program_markers": int(enriched), "marker_genes": len(all_markers),
                "surface_markers": int(ga.loc[am, "surface"].sum()), "secreted_markers": int(ga.loc[am, "secreted"].sum())})
            for res in nr.RES_GRID:
                k = f"cluster_r{res}"
                rows_res.append({"list": lst, "correction": corr, "resolution": res, "clusters_data": int(cd[k].nunique()), "clusters_null": int(cn[k].nunique()),
                                 "smallest_cluster_data": int(cd[k].value_counts().min()),
                                 "ari_half_a_half_b": adjusted_rand_score(ca[k], cb[k]), "ari_null_with_data": adjusted_rand_score(cn[k], cd[k]),
                                 "cramers_v_cluster_sample": cramers_v(cd[k], cd["sample"]),
                                 "eta2_depth_by_cluster": eta2(cd["log_depth"], cd[k]), "eta2_nuclear_share_by_cluster": eta2(cd["nuclear_share"], cd[k])})

            # UMAP figure: data (top) and null (bottom)
            rng = np.random.default_rng(nr.SEED)
            pick = rng.permutation(len(cd))[:N_PLOT]
            fig, ax = plt.subplots(2, 4, figsize=(22, 10.5))
            for i, (c_, nm) in enumerate([(cd, "data"), (cn, "null")]):
                x, y = c_["umap_1"].values[pick], c_["umap_2"].values[pick]
                size = max(0.3, 30000 / len(pick))
                for j, s in enumerate(sorted(cd["sample"].unique())):
                    mm = c_["sample"].values[pick] == s
                    ax[i, 0].scatter(x[mm], y[mm], s=size, color=plt.cm.tab10(j), label=s, rasterized=True, linewidths=0)
                ax[i, 0].legend(frameon=False, markerscale=10 / np.sqrt(size), fontsize=8)
                lab = c_[key].values[pick]
                ax[i, 1].scatter(x, y, s=size, c=lab, cmap="tab20", rasterized=True, linewidths=0)
                for k_ in np.unique(lab):
                    ax[i, 1].text(np.median(x[lab == k_]), np.median(y[lab == k_]), str(k_), fontsize=9, ha="center", weight="bold")
                for col, name, cmap in [(2, "log_depth", "viridis"), (3, "nuclear_share", "coolwarm")]:
                    v = c_[name].values[pick]
                    lo, hi = np.quantile(v, [0.02, 0.98])
                    h = ax[i, col].scatter(x, y, s=size, c=v, cmap=cmap, vmin=lo, vmax=hi, rasterized=True, linewidths=0)
                    fig.colorbar(h, ax=ax[i, col], shrink=0.7)
                for col, ttl in enumerate(["sample", f"clusters (resolution {REF_RES})", "log10 depth of the cell", "nuclear share of the cell's transcripts"]):
                    ax[i, col].set_title(f"{nm}: {ttl}"); ax[i, col].set_xticks([]); ax[i, col].set_yticks([])
            fig.suptitle(f"{LIST_NAMES[lst]}, correction {corr}", fontsize=14)
            fig.tight_layout()
            fig.savefig(out + f"umap_{lst}_{corr}.png", dpi=110)
            plt.close(fig)
            print(f"[{lst}/{corr}] done ({time.time() - t0:.0f} s)", flush=True)

    comp, by_res, markers = pd.DataFrame(rows), pd.DataFrame(rows_res), pd.DataFrame(marker_rows)
    comp.to_csv(out + "comparison.csv", index=False)
    by_res.to_csv(out + "by_resolution.csv", index=False)
    markers.to_csv(out + "cluster_markers.csv", index=False)
    keys = list(labels_ref)
    agree = pd.DataFrame([[adjusted_rand_score(labels_ref[a], labels_ref[b]) for b in keys] for a in keys],
                         index=["/".join(k) for k in keys], columns=["/".join(k) for k in keys])
    agree.to_csv(out + "list_agreement.csv")

    # overview figure
    names = [f"{LIST_NAMES[l]}\n{c}" for l, c in zip(comp["list"], comp["correction"])]
    xx = np.arange(len(comp)); w = 0.38
    panels = [("Sample effect: neighbours from the own sample / chance", "neighbours_from_own_sample_over_chance", "same_for_the_null"),
              (f"Variance explained by the first {nr.N_PCS} PCs", f"variance_first_{nr.N_PCS}_pcs_data", f"variance_first_{nr.N_PCS}_pcs_null"),
              (f"Clusters at resolution {REF_RES}", "clusters_data", "clusters_null"),
              ("Spatial pattern: same-cluster tissue neighbours / chance", "spatial_same_cluster_over_chance_data", "spatial_same_cluster_over_chance_null"),
              ("Repeatable: ARI of the clusters, half A vs half B (grey: null vs data)", "ari_half_a_half_b", "ari_null_with_data"),
              ("Largest |Spearman| of PC1-10 with depth (grey: with nuclear share)", "max_abs_rho_pc_depth", "max_abs_rho_pc_nuclear_share"),
              ("Variance between clusters: depth (grey: nuclear share)", "eta2_depth_by_cluster", "eta2_nuclear_share_by_cluster")]
    fig, ax = plt.subplots(2, 4, figsize=(26, 11))
    for a, (ttl, c1, c2) in zip(ax.ravel(), panels):
        a.bar(xx - w / 2, comp[c1], w, color="#1b9e77", label="data" if "null" in c2 else None)
        a.bar(xx + w / 2, comp[c2], w, color="#bbbbbb", label="null" if "null" in c2 and not c2.startswith("ari") else None)
        a.set_xticks(xx); a.set_xticklabels(names, fontsize=8); a.set_title(ttl, fontsize=10)
        if a.get_legend_handles_labels()[0]:
            a.legend(frameon=False)
    a = ax.ravel()[-1]
    im = a.imshow(agree.values, cmap="Blues", vmin=0, vmax=1)
    a.set_xticks(xx); a.set_xticklabels(names, fontsize=7, rotation=45); a.set_yticks(xx); a.set_yticklabels(names, fontsize=7)
    for i in range(len(keys)):
        for j in range(len(keys)):
            a.text(j, i, f"{agree.values[i, j]:.2f}", ha="center", va="center", fontsize=8, color="w" if agree.values[i, j] > 0.5 else "k")
    a.set_title(f"ARI of the clusters of the data between tasks (resolution {REF_RES})", fontsize=10)
    fig.tight_layout()
    fig.savefig(out + "overview.png", dpi=130)
    plt.close(fig)

    pd.set_option("display.width", 320); pd.set_option("display.max_columns", 60); pd.set_option("display.max_colwidth", 60)
    num = lambda v: f"{v:.0f}" if float(v).is_integer() else f"{v:.3f}"
    lines = [f"NCD ratio analysis, step 4, review of the embeddings{' (QUICK RUN)' if args.quick else ''}", "",
             f"Reference resolution {REF_RES}; {K_SPACE} tissue neighbours; markers = {N_MARKERS} genes per cluster; program enrichment p < {P_PROGRAM}", "",
             "Comparison (one column per list and correction)", comp.set_index(["list", "correction"]).T.to_string(float_format=num), "",
             "By resolution", by_res.to_string(index=False, float_format=num), "",
             "Agreement of the clusters of the data between tasks (ARI)", agree.to_string(float_format=lambda v: f"{v:.2f}"), "",
             "Clusters and their marker genes (first 12 markers)",
             markers.assign(markers=markers["markers"].str.split().str[:12].str.join(" ")).to_string(index=False, float_format=lambda v: f"{v:.3g}"), ""]
    with open(out + "summary.txt", "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    print(f"done in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
