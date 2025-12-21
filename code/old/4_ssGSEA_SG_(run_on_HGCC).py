import matplotlib.colors as clr
import matplotlib.pyplot as plt
import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from matplotlib.patches import Patch
from scipy import sparse
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# Color
color_cts = clr.LinearSegmentedColormap.from_list("magma", ["#000003", "#3B0F6F", "#8C2980", "#F66E5B", "#FD9F6C", "#FBFCBF"], N=256)

# ==================== ssGSEA functions ==================== #

# Read GMT file into dict: {pathway: [genes]}
def read_gmt(gmt_path: str) -> dict:
    gene_sets = {}
    with open(gmt_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gs_name = parts[0]
            genes = [g for g in parts[2:] if g]
            gene_sets[gs_name] = genes
    return gene_sets

# Convert gseapy ssGSEA res (res.res2d, long format) to scores matrix (sample by pathway)
def res2d_to_scores(res, score_col = "NES"):
    
    df = res.res2d.copy()

    col_map = {c.lower(): c for c in df.columns}
    name_col = col_map.get("name", "Name")
    term_col = col_map.get("term", "Term")

    score_col_actual = None
    for c in df.columns:
        if c.upper() == score_col.upper():
            score_col_actual = c
            break
    if score_col_actual is None:
        raise ValueError(f"Score column {score_col} not found. Available: {list(df.columns)}")

    scores = df.pivot(index=name_col, columns=term_col, values=score_col_actual)
    scores.index.name = "cell_id"
    return scores

# ssGSEA from cell by gene matrix (npz format)
def ssGSEA_from_cellxgene_npz_filtered(npz_path: str, cell_ids: list, gene_ids: list, gmt_path: str, out_path: str, chunk_size: int = 2000, min_geneset_size: int = 5, max_geneset_size: int = 5000, do_log1p: bool = True, do_cpm: bool = True, min_total_counts: int = 5, min_nnz_genes: int = 20):
    
    # load cell by gene matrix
    X = sparse.load_npz(npz_path).tocsr()
    if X.shape != (len(cell_ids), len(gene_ids)):
        raise ValueError(f"Shape mismatch: X {X.shape} vs {(len(cell_ids), len(gene_ids))}")

    # parse GMT into dict: {pathway: [genes]}
    gene_sets = read_gmt(gmt_path)
    pathway_names = list(gene_sets.keys())

    # ---------
    # 1) filter cells (SG-positive + enough genes)
    # ---------
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    nnz_genes = np.diff(X.indptr)  # number of nonzero genes per row in CSR

    keep = (total_counts >= min_total_counts) & (nnz_genes >= min_nnz_genes)
    keep_idx = np.where(keep)[0]
    
    print(f"ssGSEA filtering: keeping {keep_idx.size} / {X.shape[0]} cells ({keep_idx.size / X.shape[0] * 100:.2f}%)")

    # pre-allocate full scores (all cells) as zeros
    scores_full = pd.DataFrame(
        0.0,
        index=pd.Index(cell_ids, name="cell_id"),
        columns=pathway_names,
        dtype=np.float32,
    )

    # if nothing passes filtering, just write zeros and return
    if keep_idx.size == 0:
        scores_full.to_parquet(out_path)
        return scores_full

    # ---------
    # 2) run ssGSEA on kept cells only (chunk over kept_idx)
    # ---------
    for start in range(0, keep_idx.size, chunk_size):
        end = min(start + chunk_size, keep_idx.size)
        idx = keep_idx[start:end]

        Xb = X[idx, :].astype(np.float32)

        # optional: CPM + log1p to reduce ties (many zeros) and depth effects
        if do_cpm:
            libsize = np.asarray(Xb.sum(axis=1)).ravel()
            libsize[libsize == 0] = 1.0
            Xb = Xb.multiply(1e6 / libsize[:, None])
        if do_log1p:
            Xb = Xb.copy()
            Xb.data = np.log1p(Xb.data)

        # gseapy wants genes by samples (DataFrame)
        expr = pd.DataFrame(
            Xb.toarray().T,
            index=gene_ids,
            columns=[cell_ids[i] for i in idx],
        )

        res = gp.ssgsea(
            data=expr,
            gene_sets=gene_sets,
            sample_norm_method="rank",
            min_size=min_geneset_size,
            max_size=max_geneset_size,
            outdir=None,
            verbose=False,
            processes=1,
        )

        # sample by pathway
        scores = res2d_to_scores(res, score_col="NES")
        scores = scores.reindex(index=[cell_ids[i] for i in idx], columns=pathway_names)

        # write back into the full matrix (others stay 0)
        scores_full.loc[scores.index, scores.columns] = scores.astype(np.float32)

    scores_full.to_parquet(out_path)
    return scores_full

# ==================== Main operations ==================== #

settings = {"Xenium_5K_BC": {"cell_type_label": True},
            "Xenium_5K_OC": {"cell_type_label": True},
            "Xenium_5K_CC": {"cell_type_label": True},
            "Xenium_5K_LC": {"cell_type_label": False},
            "Xenium_5K_Prostate": {"cell_type_label": False},
            "Xenium_5K_Skin": {"cell_type_label": False}}

for data in settings.keys():
    
    print(f"========== Processing {data}... ==========")
    
    # paths
    data_dir = f"../../data/{data}/"
    utils_dir = "../../data/_utils/"
    output_dir = f"../../output/{data}/"
    
    # read data
    adata = sc.read_h5ad(data_dir + "intermediate_data/adata.h5ad")
    adata_tumor = adata[adata.obs["cell_type_merged"] == "Malignant cell"].copy()
    granule_adata = sc.read_h5ad(data_dir + "processed_data/granule_adata.h5ad")
    
    # determine plot size
    x_range = adata.obs["global_x"].max() - adata.obs["global_x"].min()
    y_range = adata.obs["global_y"].max() - adata.obs["global_y"].min()
    short_edge = min(x_range, y_range)

    scale = 5 / short_edge
    plot_figsize = (int(x_range * scale), int(y_range * scale))
    print(f"Plot size: {plot_figsize}")
    
    # check cell and gene IDs
    cell_ids = list(adata_tumor.obs["cell_id"])
    gene_ids = list(adata_tumor.var.index)
    
    cell_ids_npz = np.load(data_dir + "processed_data/cell_ids.npy", allow_pickle = True).tolist()
    gene_ids_npz = np.load(data_dir + "processed_data/gene_ids.npy", allow_pickle = True).tolist()
    
    if cell_ids_npz != cell_ids:
        raise ValueError("Cell ID order mismatch between NPZ and current adata_tumor!")

    if gene_ids_npz != gene_ids:
        raise ValueError("Gene order mismatch between NPZ and current adata_tumor!")
    
    # prepare data
    cell_ids = adata_tumor.obs["cell_id"].astype(str).to_numpy()
    granule_cell_ids = granule_adata.obs["cell_id"].astype(str).to_numpy()
    
    cell2row = {cid: i for i, cid in enumerate(cell_ids)}
    rows = np.fromiter((cell2row.get(cid, -1) for cid in granule_cell_ids),
                       dtype=np.int64, count=len(granule_cell_ids))
    keep = rows >= 0
    rows = rows[keep]
    
    Xg = granule_adata.X
    if not sp.isspmatrix(Xg):
        Xg = sp.csr_matrix(Xg)
    else:
        Xg = Xg.tocsr()
    Xg = Xg[keep, :]

    n_cells = adata_tumor.n_obs
    n_granules = Xg.shape[0]

    A = sp.csr_matrix(
        (np.ones(n_granules, dtype=np.float32),
        (rows, np.arange(n_granules, dtype=np.int64))),
        shape=(n_cells, n_granules),
    )
    
    X_sg = (A @ Xg).tocsr()
    sparse.save_npz(data_dir + "processed_data/SG_expression_matrix.npz", X_sg)
    print(f"Shape of the SG expression matrix: {X_sg.shape}")
    
    # run ssGSEA on SG expression
    gmt_path = utils_dir + "hallmark_pathways_filtered.gmt"
    scores = ssGSEA_from_cellxgene_npz_filtered(
            npz_path = data_dir + "processed_data/SG_expression_matrix.npz",
            cell_ids = cell_ids,
            gene_ids = gene_ids,
            gmt_path = gmt_path,
            out_path = data_dir + "processed_data/ssgsea_hallmark_sg.parquet",
        )
    
    # long format
    scores_long = scores.reset_index().melt(id_vars="cell_id", var_name="Pathway", value_name="NES")
    order = scores_long.groupby("Pathway")["NES"].median().sort_values(ascending=False).index
    order_labels = [" ".join(s.capitalize() for s in i.split("_")[1:]) for i in order]
    
    # statistical tests
    stats = []
    for pathway, df in scores_long.groupby("Pathway"):
        nes = df["NES"].dropna().to_numpy(dtype=float)
        if np.allclose(nes, 0):
            pval = 1.0
            stat = 0.0
        else:
            stat, pval = wilcoxon(nes, alternative="two-sided")
        stats.append({"Pathway": pathway, "median": np.median(nes), "pval": pval})
    stats_df = pd.DataFrame(stats)
    stats_df["qval"] = multipletests(stats_df["pval"], method="fdr_bh")[1]
    
    # determine significance
    alpha = 0.05
    stats_df["significance"] = "nonsignificant"
    stats_df.loc[(stats_df["qval"] < alpha) & (stats_df["median"] > 0), "significance"] = "positive"
    stats_df.loc[(stats_df["qval"] < alpha) & (stats_df["median"] < 0), "significance"] = "negative"
    
    scores_long = scores_long.merge(stats_df[["Pathway", "significance"]], on="Pathway", how="left")
    palette = {row.Pathway: "#d73027" if row["significance"] == "positive" else "#4575b4" if row["significance"] == "negative" else "lightgray" for _, row in stats_df.iterrows()}
    legend_handles = [Patch(facecolor="#d73027", edgecolor="black", label="Positive median NES"),
                        Patch(facecolor="#4575b4", edgecolor="black", label="Negative median NES"),
                        Patch(facecolor="lightgray", edgecolor="black", label="Not significant")]

    # boxplot of all pathways
    plt.figure(figsize=(15, 6))
    ax = sns.boxplot(data=scores_long, x="Pathway", y="NES", order=order, showfliers=False, palette=palette)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xticklabels(order_labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(handles=legend_handles, loc="upper right", frameon=True, fontsize=12)
    plt.savefig(output_dir + "ssgsea_hallmark_sg.jpeg", dpi = 300, bbox_inches = "tight")
    plt.close()
    
    # plot top pathways
    n_top = 5

    scores_mean = scores.mean(axis = 0).sort_values(ascending = False)
    top_scores = scores_mean.head(n_top)
    
    for pathway in top_scores.index:
        
        # add pathway to adata_tumor.obs
        pathway_label = f"{pathway}"
        adata_tumor.obs[pathway_label] = scores[pathway].values
        
        # plot pathway score
        sc.set_figure_params(figsize = plot_figsize)
        ax = sc.pl.scatter(adata_tumor, x="global_x", y="global_y", color=pathway_label, color_map=color_cts, size=1, show=False)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("")
        for spine in ax.spines.values():
            spine.set_visible(False)
        plt.savefig(output_dir + f"sg_{pathway_label}.jpeg", dpi = 300, bbox_inches = "tight")
        plt.close()