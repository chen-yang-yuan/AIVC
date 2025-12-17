import matplotlib.colors as clr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from mcDETECT.utils import *
from mcDETECT.model import *

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# Colors
color_cts = clr.LinearSegmentedColormap.from_list("bwr", ["#3B4CC0", "#4F69C6", "#FFFFFF", "#D24E4E", "#B40426"], N=256)
color_reds = plt.get_cmap("Reds")

# ==================== Main operations ==================== #

settings = {"Xenium_5K_BC": {"cell_type_label": True},
            "Xenium_5K_OC": {"cell_type_label": True},
            "Xenium_5K_CC": {"cell_type_label": True},
            "Xenium_5K_LC": {"cell_type_label": False},
            "Xenium_5K_Prostate": {"cell_type_label": False},
            "Xenium_5K_Skin": {"cell_type_label": False}}

plot_coords = ["global_x", "global_y"]

# for data in settings.keys():
for data in ["Xenium_5K_BC"]:
    
    print(f"========== Processing {data}... ==========")
    
    # paths
    data_dir = f"../../data/{data}/"
    utils_dir = "../../data/_utils/"
    output_dir = f"../../output/{data}/"
    
    # Read data
    genes = pd.read_csv(data_dir + "processed_data/genes.csv")
    genes = list(genes.iloc[:, 0])

    adata = sc.read_h5ad(data_dir + "intermediate_data/adata.h5ad")
    adata_tumor = adata[adata.obs["cell_type_merged"] == "Malignant cell"].copy()

    transcripts = pd.read_parquet(data_dir + "processed_data/transcripts.parquet")
    transcripts["global_z"] = 0
    print(f"Number of transcripts: {transcripts.shape[0]}")
    print("-" * 30)
    
    # Determine plot size
    x_range = adata.obs["global_x"].max() - adata.obs["global_x"].min()
    y_range = adata.obs["global_y"].max() - adata.obs["global_y"].min()
    short_edge = min(x_range, y_range)

    scale = 5 / short_edge
    plot_figsize = (int(x_range * scale), int(y_range * scale))
    print(f"Plot size: {plot_figsize}")
    print("-" * 30)
    
    # Read SG marker genes
    sg_markers_df = pd.read_excel(utils_dir + "SG_markers.xlsx")
    sg_markers_df = sg_markers_df.sort_values(by = "Fraction of RNA molecules in SGs", ascending = False)

    thr = 0.25
    sg_marker_genes = sg_markers_df[sg_markers_df["Fraction of RNA molecules in SGs"] > thr]["gene"].to_list()
    overlap_genes = [i for i in sg_marker_genes if i in genes]

    print(f"Number of SG marker genes (fraction > {thr}): {len(sg_marker_genes)}")
    print(f"Number of overlapping genes in the dataset: {len(overlap_genes)}")
    print("-" * 30)
    
    # Select transcripts only within tumor cells
    transcripts = transcripts[transcripts["cell_id"].isin(adata_tumor.obs["cell_id"])]
    print(f"Number of transcripts in tumor cells: {transcripts.shape[0]}")

    # SG transcript counts at each level
    transcripts_sg = transcripts[transcripts["target"].isin(overlap_genes)].copy()
    print(f"Number of SG transcripts in tumor cells: {transcripts_sg.shape[0]}")

    transcripts_sg_in_nucleus = transcripts_sg[transcripts_sg["in_nucleus"].astype(int) == 1].copy()
    print(f"Number of SG transcripts in tumor nucleus: {transcripts_sg_in_nucleus.shape[0]}")

    transcripts_sg_in_cytoplasm = transcripts_sg[transcripts_sg["overlaps_nucleus"].astype(int) == 0].copy()
    print(f"Number of SG transcripts in tumor cytoplasm: {transcripts_sg_in_cytoplasm.shape[0]}")

    # Non-SG transcript counts at each level
    transcripts_non_sg = transcripts[transcripts["target"].isin(overlap_genes) == False].copy()
    print(f"Number of non-SG transcripts in tumor cells: {transcripts_non_sg.shape[0]}")

    transcripts_non_sg_in_nucleus = transcripts_non_sg[transcripts_non_sg["in_nucleus"].astype(int) == 1].copy()
    print(f"Number of non-SG transcripts in tumor nucleus: {transcripts_non_sg_in_nucleus.shape[0]}")

    transcripts_non_sg_in_cytoplasm = transcripts_non_sg[transcripts_non_sg["overlaps_nucleus"].astype(int) == 0].copy()
    print(f"Number of non-SG transcripts in tumor cytoplasm: {transcripts_non_sg_in_cytoplasm.shape[0]}")

    # (Optional) Merge SG genes
    transcripts_merged = transcripts.copy()
    transcripts_merged.loc[transcripts_merged["target"].isin(overlap_genes), "target"] = "Merged"
    print("-" * 30)
    
    # SG detection
    mc = mcDETECT(type = "Xenium", transcripts = transcripts_merged, gnl_genes = ["Merged"], nc_genes = None, eps = 1,
                minspl = 3, grid_len = 1, cutoff_prob = 0.95, alpha = 10, low_bound = 3, size_thr = 4,
                in_nucleus_thr = (0.1, 0.9), l = 1, rho = 0.1, s = 1, nc_top = 15, nc_thr = 0.1)
    granules = mc.detect(record_cell_id = True)
    print(f"Granules detected: {granules.shape[0]}")
    print("-" * 30)
    
    # Assign each granule to the nearest cell
    granules["nearest_cell_type"] = adata_tumor.obs.set_index("cell_id").loc[granules["cell_id"], "cell_type_merged"].values
    granules["nearest_cell_type"] = pd.Categorical(granules["nearest_cell_type"], categories = ["Malignant cell"], ordered = True)
    
    # SG profiling
    mc = mcDETECT(type = "Xenium", transcripts = transcripts, gnl_genes = overlap_genes, nc_genes = None, eps = 1,
                minspl = 3, grid_len = 1, cutoff_prob = 0.95, alpha = 10, low_bound = 3, size_thr = 4,
                in_nucleus_thr = (0.1, 0.9), l = 1, rho = 0.1, s = 1, nc_top = 15, nc_thr = 0.1)
    granule_adata = mc.profile(granules, genes = genes, buffer = 0.05)
    
    # Average transcripts per granule
    transcripts_per_granule = np.asarray(granule_adata.X.sum(axis=1)).ravel()
    print(f"Average number of transcripts per granule: {transcripts_per_granule.mean()}")
    print("-" * 30)
    
    # SG and non-SG transcript counts in granules
    sg_gene_counts = np.asarray(granule_adata.X.sum(axis=0)).ravel()
    gene_names = granule_adata.var_names.to_numpy()

    # SG genes
    sg_mask = np.isin(gene_names, overlap_genes)
    sg_gene_counts_total = sg_gene_counts[sg_mask].sum().astype(int)
    print(f"Number of SG transcripts in granules: {sg_gene_counts_total}")

    # Non-SG genes
    non_sg_mask = sg_mask == False
    non_sg_gene_counts_total = sg_gene_counts[non_sg_mask].sum().astype(int)
    print(f"Number of non-SG transcripts in granules: {non_sg_gene_counts_total}")
    print("-" * 30)
    
    # Gene-level summary
    gene_total = transcripts.groupby("target").size().rename("total_transcripts")
    gene_sg = pd.Series(sg_gene_counts, index=gene_names, name="sg_transcripts")
    df = pd.concat([gene_total, gene_sg], axis=1).fillna(0).reset_index().rename(columns={"index": "gene"})
    df["in_SG_ratio"] = df["sg_transcripts"] / df["total_transcripts"]
    df["is_SG_marker"] = df["gene"].isin(overlap_genes).astype(int)
    df = df.sort_values("in_SG_ratio", ascending=False)
    df.to_csv(output_dir + "in_SG_ratio.csv", index = False)
    
    # Plot all granules
    sc.set_figure_params(scanpy = True, figsize = plot_figsize)
    ax = sc.pl.scatter(granule_adata, alpha = 1, x = plot_coords[0], y = plot_coords[1], color = "nearest_cell_type", palette = ["#9864bc"], size = 0.5, title = " ", show = False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_linewidth(False)
    plt.savefig(output_dir + "granules_by_cell_type.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    
    # Number of granules held by each cell
    granule_counts = granules.groupby("cell_id").size()
    adata_tumor.obs["granule_count"] = adata_tumor.obs["cell_id"].map(granule_counts).fillna(0).astype(int)
    adata_tumor.obs["log_granule_count"] = np.log1p(adata_tumor.obs["granule_count"])

    sc.set_figure_params(scanpy = True, figsize = plot_figsize)
    ax = sc.pl.scatter(adata_tumor, alpha = 1, x = plot_coords[0], y = plot_coords[1], color = "log_granule_count", color_map = color_cts, size = 1, title = " ", show = False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_linewidth(False)
    plt.savefig(output_dir + "granule_count_per_cell.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    
    # Save granules and granule_adata
    granules.to_csv(output_dir + "granules.csv", index = False)
    granule_adata.write_h5ad(output_dir + "granule_adata.h5ad")