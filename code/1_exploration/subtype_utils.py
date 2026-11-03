"""Shared functions of the tumor-cell subtyping notebooks (3_clustering.ipynb per sample, 4_pooled.ipynb pooled).

The notebooks cluster tumor cells on the depth-corrected cytoplasmic expression of the cytoplasm-enriched gene set
(embedding from 2_embedding.py / 4_pooled_embedding.py on HGCC), annotate every cluster with a biological program
derived from its enriched genes, quantify per-subtype characteristics and draw the figures. Everything that must be
identical between the per-sample and the pooled analysis lives here: the parameters, the program dictionary and
label vocabulary, program scoring, the resolution sweep and clustering, cluster diagnostics, differential
expression, the annotation rule, the characteristics tables and the figures.

Conventions
- `adata_tumor` carries the cell metadata (coordinates, areas, cell types, per-cell scores and spatial context) and
  the cluster / label columns `subtype_louvain_cyto` and `subtype_label`; `adata_tumor_cyto` is the embedded object
  (X = regressed + scaled expression, layers["lognorm"] = log-normalized, kNN graph, PCA / UMAP) with the cluster in
  `subtype_louvain`. Per sample they are two objects with identical cell order; in the pooled notebook they are the
  same object.
- A cluster id is specific to one clustering; the annotated label (vocabulary LABEL_ORDER, colours LABEL_COLORS) is
  the unit of comparison across samples and between the per-sample and the pooled analysis.
- `paths` is a dict with the output subfolders "clustering", "DE", "annotation", "spatial", "characteristics" (and
  the sample's "processed" data folder for the characteristics).
"""
import re

import matplotlib.colors as clr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import hypergeom, spearmanr
from sklearn.neighbors import NearestNeighbors

# ----------------------------------------------------------------------------------------------------------------
# Analysis parameters (shared by all samples and by the pooled analysis)
# ----------------------------------------------------------------------------------------------------------------
RES_GRID = [0.5, 0.8, 1.0, 1.2, 1.5]  # resolutions tried in the sweep
LOGFC_THR = 0.5         # DE: minimum log2 fold change of a marker (cluster vs rest)
PCT_THR = 0.1           # DE: minimum fraction of the cluster's cells expressing the marker
P_OVERLAP = 1e-3        # annotation: hypergeometric p of the DE list vs a program below which the program labels the cluster
Z_SCORE = 0.5           # annotation: mean z of a program score above which the program labels the cluster (no overlap hit)
MIN_TARGET_COUNTS = 20  # annotation: median target-gene transcripts below which a cluster is "Low depth"
BORDER_DIST = 25        # annotation: um; an EMT/stroma cluster whose median distance to fibroblasts is below this is "Stroma border"
RADIUS = 30             # characteristics: um; neighbourhood radius for the cell-type composition around each tumor cell
MIN_PROGRAM_GENES = 5   # programs with fewer genes in the gene set are dropped

# Reference cell groups for the spatial context (names of obs["cell_type_merged"]; verified by the notebooks)
spatial_refs = {"immune": ["CD4+ T cell", "CD8+ T cell", "T cell", "B cell", "Dendritic cell", "Myeloid cell", "Mast cell"],
                "tcell": ["CD4+ T cell", "CD8+ T cell", "T cell"],
                "cd8_tcell": ["CD8+ T cell"],
                "myeloid": ["Myeloid cell", "Dendritic cell"],
                "bcell": ["B cell"],
                "fibroblast": ["Fibroblast (CAF)"],
                "endothelial": ["Endothelial cell"]}

# Cell-type vocabulary of obs["cell_type_merged"] across the samples (a name outside it is reported as a warning)
known_types = {"Malignant cell", "Epithelial cell (non-malignant)", "Fibroblast (CAF)", "Smooth muscle cell", "Pericyte",
               "Endothelial cell", "Lymphatic endothelial cell", "Mesothelial cell", "Adipocyte",
               "CD4+ T cell", "CD8+ T cell", "T cell", "B cell", "Dendritic cell", "Myeloid cell", "Mast cell",
               "Mixed", "Unknown"}

# ----------------------------------------------------------------------------------------------------------------
# Colours and figure defaults
# ----------------------------------------------------------------------------------------------------------------
cmap = plt.get_cmap("tab20")
color_dct = [clr.rgb2hex(cmap(i)) for i in range(20)]  # cluster ids
color_cts = clr.LinearSegmentedColormap.from_list("magma", ["#000003", "#3B0F6F", "#8C2980", "#F66E5B", "#FD9F6C", "#FBFCBF"], N=256)  # per-cell values


def set_figsize(figsize):
    """scanpy figure defaults at the given size, with the grid that sc.set_figure_params turns on switched off again."""
    sc.set_figure_params(figsize=figsize)
    plt.rcParams["axes.grid"] = False


# ----------------------------------------------------------------------------------------------------------------
# Gene programs and label vocabulary
# ----------------------------------------------------------------------------------------------------------------
def read_gmt(path):
    """GMT file -> {set name: [genes]}."""
    gmt = {}
    for line in open(path):
        fields = line.rstrip("\n").split("\t")
        gmt[fields[0]] = [g for g in fields[2:] if g]
    return gmt


# Curated programs. The filtered GMTs do not isolate them: the hallmark interferon-gamma and complement sets share
# CFB, HIF1A, CDKN1A, CTSL and LGMN with other programs, and no hallmark set carries the HSPA8 / HSPB8 / SQSTM1
# heat-shock program, the lysosomal program (CTSL, LGMN, LIPA, GRN) or sterol synthesis (HMGCR, SQLE, DHCR24).
CURATED_PROGRAMS = {
    "Interferon": ["IFIT1", "IFIT2", "IFIT3", "IFI35", "IFIH1", "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20",
                   "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFI6", "IFI27", "IFI44", "IFI44L", "IFITM1", "IFITM3",
                   "RSAD2", "DDX58", "DDX60", "EIF2AK2", "BST2", "XAF1", "GBP1", "GBP2", "GBP4", "GBP5", "CXCL9", "CXCL10",
                   "CXCL11", "PARP9", "PARP14", "SAMD9", "SAMD9L", "HERC5", "HERC6", "USP18", "TRIM22", "APOL1", "APOL6",
                   "PLSCR1", "EPSTI1", "CMPK2", "LGALS3BP"],
    "Complement / acute phase": ["SERPINA3", "CFB", "C4A", "C4B", "UBD", "C3", "C1R", "C1S", "C2", "C5", "C7", "CFH", "CFI",
                                 "CFD", "SAA1", "SAA2", "ORM1", "ORM2", "CP", "HP", "FGA", "FGB", "FGG", "CRP", "LCN2",
                                 "SERPINA1", "SERPING1", "CD55", "CD59", "CD46", "LBP"],
    "Stress response": ["HSPA1A", "HSPA1B", "HSPA8", "HSPB8", "HSPB1", "HSPH1", "HSP90AA1", "HSP90AB1", "HSPD1", "HSPE1",
                        "DNAJA1", "DNAJB1", "DNAJB4", "BAG3", "SQSTM1", "CRYAB", "HSPA4", "HSPA6", "CHORDC1", "AHSA1",
                        "STIP1", "ATF3", "JUN", "FOS", "EGR1", "IER2", "IER3", "KLF4", "KLF6", "DUSP1", "CDKN1A", "GADD45A",
                        "GADD45B", "DDIT3", "PPP1R15A", "ATF4", "ZFAND5", "ZFAND2A", "UBC", "UBB", "HERPUD1", "SERPINH1"],
    "Lysosome": ["CTSL", "CTSB", "CTSD", "CTSK", "CTSZ", "LGMN", "LIPA", "GRN", "CREG1", "LAMP1", "LAMP2", "NPC2",
                 "PSAP", "HEXB", "HEXA", "GBA", "ASAH1", "SCARB2", "ATP6V0D1", "ATP6V1A", "GAA", "GLB1", "GLA",
                 "MAN2B1", "NEU1", "TPP1", "CD63", "LAPTM5", "LAPTM4B", "GUSB", "AP1S2", "NAGLU", "FUCA1", "GNS", "SGSH",
                 "CLN3", "MCOLN1", "TFEB"],
    "Cholesterol synthesis": ["HMGCR", "HMGCS1", "SQLE", "DHCR24", "DHCR7", "LDLR", "INSIG1", "FDFT1", "SC5D", "MVD", "IDI1",
                              "CYP51A1", "LSS", "FASN", "SCD", "SREBF2", "ACAT2", "MSMO1", "FDPS", "STARD4", "MVK", "PMVK",
                              "ACLY", "ACACA", "TM7SF2", "NSDHL", "HSD17B7", "EBP", "ACSS2"],
    # macrophage-specific genes: a cluster enriched for them carries transcripts of neighbouring myeloid cells
    # (segmentation spill-over into the cytoplasmic compartment), the myeloid counterpart of the stroma border
    "Myeloid border": ["ITGB2", "SIRPA", "MSR1", "CD68", "CYBB", "SLC40A1", "CD163", "C1QA", "C1QB", "C1QC", "CD14", "LYZ",
                       "AIF1", "TYROBP", "FCGR3A", "CSF1R", "MS4A7", "MS4A6A", "PTPRC", "LAPTM5", "CTSS", "FCER1G", "CD53",
                       "MRC1", "MARCO", "VSIG4", "FOLR2", "LILRB4", "SPI1", "CD86", "ITGAM", "ITGAX", "HLA-DRA", "CD74",
                       "LST1", "CORO1A", "MPEG1", "SLCO2B1", "APOE", "APOC1", "FCGR2A", "CD4", "LYVE1", "STAB1", "TREM2",
                       "GPNMB", "MS4A4A"],
}

# Program name -> (source, set name) for the GMT-derived programs ("hallmark" / "gobp")
GMT_PROGRAMS = {
    "S-phase": ("hallmark", "HALLMARK_E2F_TARGETS"),
    "G2/M": ("hallmark", "HALLMARK_G2M_CHECKPOINT"),
    "Hypoxia": ("hallmark", "HALLMARK_HYPOXIA"),
    "Glycolysis": ("hallmark", "HALLMARK_GLYCOLYSIS"),
    "EMT/stroma": ("hallmark", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
    "ER stress": ("gobp", "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"),
    "Oxidative stress": ("gobp", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS"),
    "Autophagy": ("gobp", "GOBP_MACROAUTOPHAGY"),
    "TNFA/NFKB": ("hallmark", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    "Apoptosis": ("hallmark", "HALLMARK_APOPTOSIS"),
    "p53": ("hallmark", "HALLMARK_P53_PATHWAY"),
    "mTORC1": ("hallmark", "HALLMARK_MTORC1_SIGNALING"),
    "MYC targets": ("hallmark", "HALLMARK_MYC_TARGETS_V1"),
    # scores only, never used as a label (tissue identity or sets that overlap the curated programs)
    "Estrogen early": ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_EARLY"),
    "Estrogen late": ("hallmark", "HALLMARK_ESTROGEN_RESPONSE_LATE"),
    "Interferon gamma (hallmark)": ("hallmark", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    "Complement (hallmark)": ("hallmark", "HALLMARK_COMPLEMENT"),
}

# Order of the full dictionary (defines LABEL_ORDER and the colours; must not change between runs)
PROGRAM_NAMES = ["S-phase", "G2/M", "Interferon", "Complement / acute phase", "Stress response", "Lysosome",
                 "Cholesterol synthesis", "Hypoxia", "Glycolysis", "EMT/stroma", "ER stress", "Oxidative stress",
                 "Autophagy", "TNFA/NFKB", "Apoptosis", "p53", "mTORC1", "MYC targets", "Myeloid border",
                 "Estrogen early", "Estrogen late", "Interferon gamma (hallmark)", "Complement (hallmark)"]
score_only = ["Estrogen early", "Estrogen late", "Interferon gamma (hallmark)", "Complement (hallmark)"]

# Fixed label vocabulary and colours, shared by all samples: the label programs plus three non-program classes
LABEL_ORDER = [name for name in PROGRAM_NAMES if name not in score_only + ["Myeloid border"]] + ["Baseline", "Stroma border", "Myeloid border", "Low depth"]
palette_programs = [color_dct[i] for i in range(20) if i not in (14, 15)]  # tab20 without the two greys
LABEL_COLORS = {name: palette_programs[i % len(palette_programs)] for i, name in enumerate(LABEL_ORDER[:-4])}
LABEL_COLORS.update({"Baseline": "#999999", "Stroma border": "#8c564b", "Myeloid border": "#c49c94", "Low depth": "#d9d9d9"})


def load_programs_raw(utils_path):
    """The full program dictionary {name: [genes]} in PROGRAM_NAMES order (GMT sets from data/_utils + curated lists)."""
    gmt = {"hallmark": read_gmt(utils_path + "hallmark_pathways_filtered.gmt"),
           "gobp": read_gmt(utils_path + "all_pathways_filtered.gmt")}
    programs_raw = {}
    for name in PROGRAM_NAMES:
        if name in CURATED_PROGRAMS:
            programs_raw[name] = list(CURATED_PROGRAMS[name])
        else:
            source, set_name = GMT_PROGRAMS[name]
            programs_raw[name] = gmt[source][set_name]
    return programs_raw


def slug(label):
    """Label -> file-name-safe string (e.g. 'Complement / acute phase' -> 'complement_acute_phase')."""
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def build_programs(gene_set, programs_raw):
    """Restrict every program to the embedding gene set; drop the ones that become too small.
    Returns (programs, label_programs): the restricted dictionary and the names that may be used as a label."""
    programs = {}
    for name, genes in programs_raw.items():
        genes = [g for g in genes if g in gene_set]
        if len(genes) < MIN_PROGRAM_GENES:
            print(f"{name}: only {len(genes)} genes in the gene set, dropped")
            continue
        programs[name] = genes
        print(f"{name:28s} {len(genes):3d} genes: {', '.join(genes[:10])}{', ...' if len(genes) > 10 else ''}")
    label_programs = [name for name in programs if name not in score_only]
    return programs, label_programs


def score_programs(adata_tumor_cyto, adata_tumor, programs):
    """Per-cell program scores on the log-normalized cytoplasmic expression -> z-scored DataFrame (cells x programs).
    score_genes has no layer argument, so it runs on a light copy whose X is the lognorm layer (adata_tumor_cyto.X holds
    the regressed and scaled values used for the embedding, which are not suitable). Raw scores are also stored in
    adata_tumor.obs as score_<program>."""
    adata_scores = sc.AnnData(X=adata_tumor_cyto.layers["lognorm"], obs=adata_tumor_cyto.obs[[]].copy(), var=adata_tumor_cyto.var[[]].copy())
    for name, genes in programs.items():
        sc.tl.score_genes(adata_scores, genes, score_name=name, ctrl_size=50, n_bins=25, use_raw=False, random_state=0)
    scores = adata_scores.obs[list(programs)].astype(float)
    for name in programs:
        adata_tumor.obs[f"score_{name}"] = scores[name].values
    return (scores - scores.mean()) / scores.std()


def top_program(scores_z, label_programs, cluster_labels):
    """Mean program z per cluster (label programs only) -> DataFrame(top_program, top_z) indexed by cluster."""
    cl_z = scores_z[label_programs].groupby(np.asarray(cluster_labels)).mean()
    return pd.DataFrame({"top_program": cl_z.idxmax(axis=1), "top_z": cl_z.max(axis=1)})


# ----------------------------------------------------------------------------------------------------------------
# Resolution sweep and clustering
# ----------------------------------------------------------------------------------------------------------------
def resolution_sweep(adata_tumor_cyto, scores_z, label_programs, cfg, paths):
    """Louvain at every resolution of RES_GRID on the stored kNN graph. For each value: number of clusters, smallest
    cluster, and the number of distinct programs that the score criterion alone (mean z > Z_SCORE) resolves. The
    suggested resolution is the smallest one at which that number stops increasing, i.e. the coarsest clustering
    that still separates all detectable programs; cfg["res_cyto"] is the reviewed value actually used."""
    sweep = []
    for res in RES_GRID:
        key = f"louvain_{res}"
        sc.tl.louvain(adata_tumor_cyto, resolution=res, key_added=key)
        tp = top_program(scores_z, label_programs, adata_tumor_cyto.obs[key].values)
        sizes = adata_tumor_cyto.obs[key].value_counts()
        resolved = tp.loc[tp["top_z"] > Z_SCORE, "top_program"]
        sweep.append({"resolution": res, "n_clusters": len(sizes), "min_cluster_size": int(sizes.min()),
                      "n_program_labels": resolved.nunique(), "programs": ", ".join(sorted(resolved.unique()))})
    sweep = pd.DataFrame(sweep)
    sweep.to_csv(paths["clustering"] + "resolution_sweep_cytoplasmic.csv", index=False)
    res_suggested = sweep.loc[sweep["n_program_labels"] == sweep["n_program_labels"].max(), "resolution"].min()
    print(sweep.to_string(index=False))
    print(f"Suggested resolution: {res_suggested} (config uses {cfg['res_cyto']})")
    return sweep


def cluster(adata_tumor_cyto, adata_tumor, res_cyto):
    """Louvain at the reviewed resolution -> subtype_louvain (embedded object) and subtype_louvain_cyto (metadata object)."""
    sc.tl.louvain(adata_tumor_cyto, resolution=res_cyto, key_added="subtype_louvain")
    adata_tumor.obs["subtype_louvain_cyto"] = adata_tumor_cyto.obs["subtype_louvain"].values
    print(f"Number of clusters (cytoplasmic, resolution {res_cyto}): {adata_tumor.obs['subtype_louvain_cyto'].nunique()}")


# ----------------------------------------------------------------------------------------------------------------
# Cluster diagnostics: spatial context and depth
# ----------------------------------------------------------------------------------------------------------------
def spatial_context(obs_all, obs_tumor, xy=("global_x", "global_y")):
    """Spatial context of every tumor cell within ONE sample (coordinates must share one frame).
    obs_all: all cells of the sample (columns xy and cell_type_merged); obs_tumor: the tumor cells (columns xy).
    Returns a DataFrame indexed like obs_tumor with, per reference group of spatial_refs, the distance (um) to the
    nearest reference cell and the proximity = -log10(distance + 1) (higher = closer; the log compresses the long
    tail of far cells), and the neighbourhood within RADIUS um: number of cells (crowding) and the immune / T-cell
    fraction. A group that is absent from the sample gets NaN distances / proximities."""
    xy_all = obs_all[list(xy)].values
    xy_tumor = obs_tumor[list(xy)].values
    out = pd.DataFrame(index=obs_tumor.index)
    for name, cell_types in spatial_refs.items():
        is_ref = obs_all["cell_type_merged"].isin(cell_types).values
        if is_ref.sum() == 0:
            print(f"{name}: no cells of types {cell_types} in this sample, distance set to NaN")
            out[f"dist_{name}"] = np.nan
        else:
            out[f"dist_{name}"] = NearestNeighbors(n_neighbors=1).fit(xy_all[is_ref]).kneighbors(xy_tumor)[0].ravel()
        out[f"proximity_{name}"] = -np.log10(out[f"dist_{name}"] + 1)

    # neighbourhood composition within RADIUS um (all cells, the tumor cell itself included)
    neighbors = NearestNeighbors(radius=RADIUS).fit(xy_all).radius_neighbors(xy_tumor, return_distance=False)
    out[f"n_cells_{RADIUS}um"] = [len(idx) for idx in neighbors]
    for name in ["immune", "tcell"]:
        is_ref = obs_all["cell_type_merged"].isin(spatial_refs[name]).values
        out[f"frac_{name}_{RADIUS}um"] = [is_ref[idx].mean() if len(idx) else np.nan for idx in neighbors]
    return out


def cluster_statistics(adata_tumor, adata_tumor_cyto, paths):
    """Depth columns, the depth check on the embedding and the per-cluster statistics table.
    Cytoplasmic depth (full panel) comes from obs["compartment_total_counts"] of the embedded object; the raw
    target-gene counts are recovered from the log-normalized layer. Spearman(PC1-5, log depth) should stay below 0.3
    in absolute value (the embedding regressed depth out). The table (clustering/cluster_stats_cytoplasmic.csv)
    feeds the annotation overrides: median target-gene counts (Low depth) and median fibroblast distance (Stroma
    border); it requires the spatial-context columns dist_immune and dist_fibroblast in adata_tumor.obs."""
    depth = adata_tumor_cyto.obs["compartment_total_counts"].values
    lognorm = adata_tumor_cyto.layers["lognorm"]
    lognorm = lognorm.toarray() if sparse.issparse(lognorm) else lognorm
    adata_tumor.obs["cyto_depth"] = depth
    adata_tumor.obs["cyto_target_counts"] = np.expm1(lognorm).sum(axis=1) * depth / 1e4

    for k in range(5):
        rho = spearmanr(adata_tumor_cyto.obsm["X_pca"][:, k], np.log10(depth + 1)).correlation
        print(f"Spearman(PC{k + 1}, log10 cytoplasmic depth) = {rho:+.2f}")

    groups = adata_tumor.obs.groupby("subtype_louvain_cyto", observed=True)
    cluster_stats = pd.DataFrame({
        "n_cells": groups.size(),
        "median_cyto_depth": groups["cyto_depth"].median(),
        "median_cyto_target_counts": groups["cyto_target_counts"].median(),
        "median_cell_area": groups["cell_area"].median(),
        "median_dist_immune": groups["dist_immune"].median(),
        "median_dist_fibroblast": groups["dist_fibroblast"].median(),
    })
    cluster_stats["low_depth"] = cluster_stats["median_cyto_target_counts"] < MIN_TARGET_COUNTS
    cluster_stats = cluster_stats.round(1)
    cluster_stats.to_csv(paths["clustering"] + "cluster_stats_cytoplasmic.csv")
    print(cluster_stats.to_string())
    return cluster_stats


def diagnostics(adata, adata_tumor, adata_tumor_cyto, paths):
    """Per-sample diagnostics: spatial context from the full sample object, then the cluster statistics."""
    context = spatial_context(adata.obs, adata_tumor.obs)
    for col in context.columns:
        adata_tumor.obs[col] = context[col].values
    return cluster_statistics(adata_tumor, adata_tumor_cyto, paths)


# ----------------------------------------------------------------------------------------------------------------
# Differential expression
# ----------------------------------------------------------------------------------------------------------------
def differential_expression(adata_tumor_cyto, paths):
    """Wilcoxon test of each cluster against all other tumor cells on layers["lognorm"]. With ~1e5 cells almost every
    gene is significant, so a marker must also pass log2 fold change > LOGFC_THR and be expressed in more than
    PCT_THR of the cluster's cells. One table per cluster in DE/, sorted by the Wilcoxon score; returns {cluster: table}."""
    sc.tl.rank_genes_groups(adata_tumor_cyto, groupby="subtype_louvain", method="wilcoxon", layer="lognorm", use_raw=False, pts=True)
    de_dict = {}
    for cl in adata_tumor_cyto.obs["subtype_louvain"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata_tumor_cyto, group=cl)
        df = df[(df["logfoldchanges"] > LOGFC_THR) & (df["pvals_adj"] < 0.05) & (df["pct_nz_group"] > PCT_THR)]
        df = df.sort_values("scores", ascending=False)
        df.to_csv(paths["DE"] + f"DE_cytoplasmic_cluster_{cl}.csv", index=False)
        de_dict[cl] = df
        print(f"Cluster {cl}: {df.shape[0]} marker genes; top: {', '.join(df['names'].head(8))}")
    return de_dict


# ----------------------------------------------------------------------------------------------------------------
# Annotation
# ----------------------------------------------------------------------------------------------------------------
def annotate(adata_tumor_cyto, adata_tumor, de_dict, cluster_stats, scores_z, programs, label_programs, cfg, paths):
    """Label every cluster from its enriched genes with two independent criteria:
    1. marker overlap: hypergeometric test of the cluster's marker list against each label program (universe = gene
       set); the program with the smallest p labels the cluster if p < P_OVERLAP;
    2. program score: if no program passes the overlap test, the top-scoring program (mean z over the cluster's
       cells) labels the cluster if its mean z > Z_SCORE.
    Neither -> "Baseline" (tissue identity without a process program). Overrides: "Low depth" if the median
    target-gene count is below MIN_TARGET_COUNTS; "Stroma border" if the label is EMT/stroma and the median distance
    to fibroblasts is below BORDER_DIST um (spill-over of adjacent fibroblast transcripts). The curated "Myeloid
    border" program is a label of its own for macrophage spill-over.
    confidence: "high" when both criteria name the same program (or a Baseline / Low depth cluster has fewer than
    10 markers), "medium" when only one criterion fires, "review" when the criteria disagree or a Baseline / Low
    depth cluster has 10 or more markers. Reviewed decisions in cfg["manual_labels"] (cluster id -> label) are
    applied last and marked "reviewed".
    Writes annotation/subtype_annotation_cytoplasmic.csv, sets subtype_label on both objects (categories in
    LABEL_ORDER order), prints the table and the top markers per label; returns (annotation, categories)."""
    n_universe = adata_tumor_cyto.shape[1]
    tp = top_program(scores_z, label_programs, adata_tumor_cyto.obs["subtype_louvain"].values)

    rows = []
    for cl in adata_tumor_cyto.obs["subtype_louvain"].cat.categories:
        markers = set(de_dict[cl]["names"])

        # criterion 1: overlap of the marker list with each label program (hypergeometric, universe = gene set)
        overlap = {}
        for name in label_programs:
            k = len(markers & set(programs[name]))
            overlap[name] = (k, hypergeom.sf(k - 1, n_universe, len(programs[name]), len(markers)) if markers else 1.0)
        best = min(overlap, key=lambda name: overlap[name][1])
        by_overlap = best if overlap[best][1] < P_OVERLAP else None

        # criterion 2: top program score
        by_score = tp.loc[cl, "top_program"] if tp.loc[cl, "top_z"] > Z_SCORE else None

        # rule and confidence
        if by_overlap and by_score:
            label, confidence = by_overlap, ("high" if by_overlap == by_score else "review")
        elif by_overlap or by_score:
            label, confidence = (by_overlap or by_score), "medium"
        else:
            label, confidence = "Baseline", ("review" if len(markers) >= 10 else "high")

        # overrides
        stats = cluster_stats.loc[cl]
        if stats["median_cyto_target_counts"] < MIN_TARGET_COUNTS:
            label, confidence = "Low depth", ("review" if len(markers) >= 10 else "high")
        elif label == "EMT/stroma" and stats["median_dist_fibroblast"] < BORDER_DIST:
            label = "Stroma border"

        # reviewed decisions
        if cl in cfg["manual_labels"]:
            label, confidence = cfg["manual_labels"][cl], "reviewed"

        rows.append({"cluster": cl, "n_cells": int(stats["n_cells"]), "label": label, "confidence": confidence,
                     "n_markers": len(markers),
                     "overlap_program": best, "overlap_genes": overlap[best][0], "overlap_p": overlap[best][1],
                     "score_program": tp.loc[cl, "top_program"], "score_z": tp.loc[cl, "top_z"],
                     "top_markers": ", ".join(de_dict[cl]["names"].head(10))})
    annotation = pd.DataFrame(rows).set_index("cluster")
    annotation.to_csv(paths["annotation"] + "subtype_annotation_cytoplasmic.csv")

    # cluster -> label; the category order and colours are fixed across samples
    label_map = annotation["label"].to_dict()
    categories = [label for label in LABEL_ORDER if label in set(label_map.values())]
    assert len(categories) == len(set(label_map.values())), "a label is missing from LABEL_ORDER / LABEL_COLORS"
    adata_tumor.obs["subtype_label"] = pd.Categorical(adata_tumor.obs["subtype_louvain_cyto"].astype(str).map(label_map), categories=categories)
    adata_tumor_cyto.obs["subtype_label"] = adata_tumor.obs["subtype_label"].values

    to_review = annotation.index[annotation["confidence"] == "review"].tolist()
    print(f"{len(annotation)} clusters -> {len(categories)} labels; clusters to review: {to_review if to_review else 'none'}")
    pd.set_option("display.width", 200)
    print(annotation.drop(columns="top_markers").round({"overlap_p": 12, "score_z": 2}).to_string())

    # top markers per label (union of its clusters, ordered by Wilcoxon score)
    for label in categories:
        clusters = annotation.index[annotation["label"] == label].tolist()
        df = pd.concat([de_dict[c] for c in clusters]).sort_values("scores", ascending=False)
        print(f"{label} (clusters {', '.join(clusters)}; {adata_tumor.obs['subtype_label'].value_counts()[label]} cells): "
              f"{', '.join(list(dict.fromkeys(df['names']))[:15])}")
    return annotation, categories


def save_cell_subtypes(adata_tumor, annotation, sample, res_cyto, path):
    """Compact per-cell record of the subtyping for later stages, keyed by sample + cell_id:
    sample, cell_id, subtype_louvain_cyto, subtype_label, confidence (of the cluster's label), resolution."""
    table = pd.DataFrame({"sample": sample,
                          "cell_id": adata_tumor.obs["cell_id"].values,
                          "subtype_louvain_cyto": adata_tumor.obs["subtype_louvain_cyto"].astype(str).values,
                          "subtype_label": adata_tumor.obs["subtype_label"].astype(str).values})
    table["confidence"] = table["subtype_louvain_cyto"].map(annotation["confidence"].to_dict())
    table["resolution"] = res_cyto
    table.to_csv(path, index=False)
    return table


# ----------------------------------------------------------------------------------------------------------------
# Subtype characteristics
# ----------------------------------------------------------------------------------------------------------------
def cell_characteristics(adata_tumor, paths):
    """Per-cell characteristics that need the sample's processed data or whole-cell counts: the ssGSEA hypoxia and
    interferon-gamma scores of 2_outcome (cytoplasmic), the 8-TF hypoxia stress score of 2_outcome/5_score.ipynb
    (2nd-highest TF on whole-cell normalize_total + log1p; ~0 for most cells at Xenium depth, kept for
    comparability) and the nucleus-to-cell area ratio."""
    ssgsea = pd.read_parquet(paths["processed"] + "ssgsea_hallmark_cytoplasmic.parquet")
    ssgsea = ssgsea.loc[adata_tumor.obs["cell_id"].values]
    adata_tumor.obs["hypoxia_ssgsea_cyto"] = ssgsea["HALLMARK_HYPOXIA"].values
    adata_tumor.obs["interferon_ssgsea_cyto"] = ssgsea["HALLMARK_INTERFERON_GAMMA_RESPONSE"].values

    hypoxia_tf = ["HIF1A", "EPAS1", "NFE2L2", "CREB1", "RELA", "RELB", "NFKB1", "NFKB2"]
    adata_wc = adata_tumor.copy()
    sc.pp.normalize_total(adata_wc, target_sum=1e4)
    sc.pp.log1p(adata_wc)
    X_tf = adata_wc[:, [g for g in hypoxia_tf if g in adata_wc.var_names]].X
    X_tf = X_tf.toarray() if sparse.issparse(X_tf) else np.asarray(X_tf)
    adata_tumor.obs["hypoxia_score_tf"] = np.sort(X_tf, axis=1)[:, -2]
    del adata_wc

    adata_tumor.obs["nucleus_cell_area_ratio"] = adata_tumor.obs["nucleus_area"] / adata_tumor.obs["cell_area"]

    for a, b in [("score_Hypoxia", "hypoxia_ssgsea_cyto"), ("score_Hypoxia", "hypoxia_score_tf"), ("hypoxia_ssgsea_cyto", "hypoxia_score_tf")]:
        print(f"Spearman({a}, {b}) = {spearmanr(adata_tumor.obs[a], adata_tumor.obs[b]).correlation:+.2f}")


CHARACTERISTIC_COLS = ["score_Hypoxia", "hypoxia_ssgsea_cyto", "hypoxia_score_tf", "score_Interferon", "interferon_ssgsea_cyto",
                       "score_Complement / acute phase", "score_S-phase", "score_G2/M", "score_Stress response", "score_Lysosome"] + \
                      [f"proximity_{name}" for name in spatial_refs] + \
                      [f"frac_immune_{RADIUS}um", f"frac_tcell_{RADIUS}um", f"n_cells_{RADIUS}um",
                       "cyto_depth", "cyto_target_counts", "cell_area", "nucleus_cell_area_ratio"]


def summarize_characteristics(adata_tumor, paths, extra_cell_cols=()):
    """Per-label summary (median and IQR of CHARACTERISTIC_COLS) -> characteristics/subtype_characteristics_cytoplasmic.csv,
    and the per-cell table (scores, distances and proximities, neighbourhood, depth, morphology) ->
    characteristics/cell_annotations_cytoplasmic.csv. extra_cell_cols: additional obs columns for the per-cell table."""
    char_cols = [c for c in CHARACTERISTIC_COLS if c in adata_tumor.obs]
    groups = adata_tumor.obs.groupby("subtype_label", observed=True)
    summary = groups[char_cols].median()
    summary.insert(0, "n_cells", groups.size())
    iqr = (groups[char_cols].quantile(0.75) - groups[char_cols].quantile(0.25)).add_suffix("_iqr")
    pd.concat([summary, iqr], axis=1).round(4).to_csv(paths["characteristics"] + "subtype_characteristics_cytoplasmic.csv")

    cell_cols = list(extra_cell_cols) + ["cell_id", "subtype_louvain_cyto", "subtype_label"] + \
                [c for c in adata_tumor.obs.columns if c.startswith("score_")] + \
                [f"dist_{name}" for name in spatial_refs] + [c for c in char_cols if not c.startswith("score_")]
    adata_tumor.obs[cell_cols].to_csv(paths["characteristics"] + "cell_annotations_cytoplasmic.csv", index=False)

    pd.set_option("display.max_columns", 40)
    print(summary.round(3).to_string())
    return summary


def characteristics(adata_tumor, paths):
    """Per-sample characteristics: per-cell quantities, then the per-label summary and the per-cell table."""
    cell_characteristics(adata_tumor, paths)
    return summarize_characteristics(adata_tumor, paths)


# ----------------------------------------------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------------------------------------------
def make_figures(obs_all, adata_tumor, adata_tumor_cyto, de_dict, annotation, categories, scores_z, title, paths, fig_size,
                 xy=("global_x", "global_y")):
    """All figures of one analysis (see the notebook's figure table). obs_all: all cells (tumor and other) with the
    coordinate columns xy and cell_type_merged, used as the immune / T-cell backgrounds of the per-label tissue maps;
    xy: coordinate columns of the tissue maps (the sample's global coordinates, or canvas coordinates in the pooled
    notebook); title: sample name for the figure titles; fig_size: (width, height) of the tissue maps."""
    label_palette = [LABEL_COLORS[label] for label in categories]
    x, y = xy

    # UMAPs and PCA loadings
    set_figsize((8, 8))
    sc.pl.umap(adata_tumor_cyto, color="subtype_louvain", size=10, palette=color_dct, show=False, title="Cytoplasmic expression")
    plt.savefig(paths["clustering"] + "umap_clustering.png", dpi=300, bbox_inches="tight")
    plt.close()
    sc.pl.umap(adata_tumor_cyto, color="subtype_label", size=10, palette=label_palette, show=False, title="Cytoplasmic expression")
    plt.savefig(paths["clustering"] + "umap_subtype_label.png", dpi=300, bbox_inches="tight")
    plt.close()
    sc.pl.pca_loadings(adata_tumor_cyto, components=[1, 2, 3, 4, 5], show=False)
    plt.savefig(paths["clustering"] + "pca_loadings.png", dpi=200, bbox_inches="tight")
    plt.close()

    # tissue map by label
    set_figsize(fig_size)
    ax = sc.pl.scatter(adata_tumor, x=x, y=y, color="subtype_label", palette=label_palette, size=2.5, alpha=1, show=False)
    ax.grid(False)
    plt.savefig(paths["clustering"] + "scatter_subtype_label.png", dpi=300, bbox_inches="tight")
    plt.close()

    # label composition
    counts = adata_tumor.obs["subtype_label"].value_counts().loc[categories]
    fig, ax = plt.subplots(figsize=(0.6 * len(categories) + 2, 4))
    ax.bar(categories, counts.values / counts.sum() * 100, color=label_palette)
    ax.grid(False)
    ax.set_ylabel("% of tumor cells")
    ax.set_title(title)
    ax.set_xticklabels(categories, rotation=60, ha="right")
    plt.savefig(paths["clustering"] + "subtype_composition.png", dpi=300, bbox_inches="tight")
    plt.close()

    # program-score heatmap (labels x programs)
    label_z = scores_z.groupby(adata_tumor.obs["subtype_label"].values).mean().loc[categories]
    fig, ax = plt.subplots(figsize=(0.55 * label_z.shape[1] + 3, 0.5 * label_z.shape[0] + 2))
    sns.heatmap(label_z, cmap="RdBu_r", center=0, vmin=-1.5, vmax=1.5, annot=True, fmt=".1f", annot_kws={"size": 7},
                linewidths=0, cbar_kws={"label": "mean z of program score"}, ax=ax)
    ax.grid(False)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title)
    plt.savefig(paths["annotation"] + "heatmap_program_scores.png", dpi=300, bbox_inches="tight")
    plt.close()

    # marker dotplot: top 5 markers of each label (union of its clusters)
    marker_genes = {}
    for label in categories:
        clusters = annotation.index[annotation["label"] == label].tolist()
        df = pd.concat([de_dict[c] for c in clusters]).sort_values("scores", ascending=False)
        genes = list(dict.fromkeys(df["names"]))[:5]
        if genes:
            marker_genes[label] = genes
    sc.pl.dotplot(adata_tumor_cyto, var_names=marker_genes, groupby="subtype_label", layer="lognorm", use_raw=False,
                  standard_scale="var", show=False)
    plt.savefig(paths["DE"] + "dotplot_subtype_markers.png", dpi=300, bbox_inches="tight")
    plt.close()

    # one tissue map per label: alone on the other tumor cells, on the immune cells, on the T cells
    backgrounds = {"": adata_tumor.obs, "immune_cell": obs_all[obs_all["cell_type_merged"].isin(spatial_refs["immune"])],
                   "tcell": obs_all[obs_all["cell_type_merged"].isin(spatial_refs["tcell"])]}
    for label in categories:
        is_target = (adata_tumor.obs["subtype_label"] == label).values
        for background_name, background in backgrounds.items():
            background = background if background_name else background.loc[~is_target]
            fig, ax = plt.subplots(figsize=fig_size)
            ax.scatter(background[x], background[y], c="lightgray", s=0.5, alpha=0.8, linewidths=0)
            ax.scatter(adata_tumor.obs.loc[is_target, x], adata_tumor.obs.loc[is_target, y],
                       c=LABEL_COLORS[label], s=1, alpha=0.8, linewidths=0)
            ax.grid(False)
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            ax.set_title(f"{title}: {label} ({is_target.sum()} cells)" + (f" on {background_name.replace('_', ' ')}s" if background_name else ""))
            suffix = f"_{background_name}" if background_name else ""
            plt.savefig(paths["spatial"] + f"scatter_{slug(label)}{suffix}.png", dpi=300, bbox_inches="tight")
            plt.close()

    # box plots of the characteristics per label (outliers hidden)
    box_cols = {"score_Hypoxia": ("Hypoxia score (hallmark, cytoplasmic)", ""),
                "hypoxia_ssgsea_cyto": ("Hypoxia score (ssGSEA, cytoplasmic)", ""),
                "score_Interferon": ("Interferon score", ""),
                "proximity_immune": ("Proximity to nearest immune cell", "-log10(um + 1)"),
                "proximity_tcell": ("Proximity to nearest T cell", "-log10(um + 1)"),
                "proximity_myeloid": ("Proximity to nearest myeloid cell", "-log10(um + 1)"),
                f"frac_immune_{RADIUS}um": (f"Immune fraction within {RADIUS} um", ""),
                f"n_cells_{RADIUS}um": (f"Cells within {RADIUS} um", ""),
                "cell_area": ("Cell area", "um2"),
                "nucleus_cell_area_ratio": ("Nucleus / cell area", "")}
    fig, axes = plt.subplots(2, 5, figsize=(26, 11))
    for ax, (col, (box_title, unit)) in zip(axes.ravel(), box_cols.items()):
        if adata_tumor.obs[col].notna().any():
            sns.boxplot(data=adata_tumor.obs, x="subtype_label", y=col, order=categories, palette=LABEL_COLORS, showfliers=False, ax=ax)
        ax.grid(False)
        ax.set_title(box_title)
        ax.set_xlabel("")
        ax.set_ylabel(unit)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=60, ha="right")
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(paths["characteristics"] + "boxplot_characteristics.png", dpi=200, bbox_inches="tight")
    plt.close()

    # tissue maps of per-cell values (clipped at the 1st / 99th percentile for the colour scale; high = bright)
    set_figsize(fig_size)
    for col, fname in {"score_Hypoxia": "scatter_score_hypoxia", "score_Interferon": "scatter_score_interferon",
                       "proximity_immune": "scatter_proximity_immune", "proximity_tcell": "scatter_proximity_tcell"}.items():
        if adata_tumor.obs[col].isna().all():
            continue
        lo, hi = adata_tumor.obs[col].quantile([0.01, 0.99])
        adata_tumor.obs["_plot"] = adata_tumor.obs[col].clip(lo, hi)
        ax = sc.pl.scatter(adata_tumor, x=x, y=y, color="_plot", color_map=color_cts, size=2.5, title=f"{title}: {col}", show=False)
        ax.grid(False)
        plt.savefig(paths["characteristics"] + f"{fname}.png", dpi=300, bbox_inches="tight")
        plt.close()
    adata_tumor.obs.drop(columns="_plot", inplace=True, errors="ignore")
