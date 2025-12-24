import anndata
import numpy as np
import pandas as pd
from collections import Counter
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def make_tree(d1 = None, d2 = None, d3 = None):
    active_dimensions = [dimension for dimension in [d1, d2, d3] if dimension is not None]
    if len(active_dimensions) == 1:
        points = np.c_[active_dimensions[0].ravel()]
    elif len(active_dimensions) == 2:
        points = np.c_[active_dimensions[0].ravel(), active_dimensions[1].ravel()]
    elif len(active_dimensions) == 3:
        points = np.c_[active_dimensions[0].ravel(), active_dimensions[1].ravel(), active_dimensions[2].ravel()]
    return cKDTree(points)


def profile(granules, transcripts, genes = None, radius_alpha = 2, radius_beta = 10, radius_shift = 0.6, print_itr = False):
    
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    gene_array = transcripts["target"].to_numpy()
    tree = make_tree(d1 = np.array(transcripts["global_x"]), d2 = np.array(transcripts["global_y"]), d3 = np.array(transcripts["global_z"]))
    
    n_gnl = granules.shape[0]
    n_gene = len(genes)
    data, row_idx, col_idx = [], [], []
    
    # randomly sample radius from Beta distribution
    np.random.seed(42)
    radius = np.random.beta(radius_alpha, radius_beta, size = n_gnl) 
    radius += radius_shift
    
    # iterate over all granules to count nearby transcripts
    for i in range(n_gnl):
        temp = granules.iloc[i]
        target_idx = tree.query_ball_point([temp["global_x"], temp["global_y"], temp["global_z"]], radius[i])
        if not target_idx:
            continue
        local_genes = gene_array[target_idx]    # extract genes for those nearby transcripts
        counts = Counter(local_genes)           # count how many times each gene occurs
        for g, cnt in counts.items():           # append nonzero entries to sparse matrix lists
            j = gene_to_idx[g]                  # get gene column index
            data.append(cnt)                    # nonzero count
            row_idx.append(i)                   # row index = granule index
            col_idx.append(j)                   # column index = gene index
        if print_itr and (i % 5000 == 0):
            print(f"{i} out of {n_gnl} granules profiled!")
    
    # construct sparse spatial transcriptome profile, (n_granules × n_genes)
    X = csr_matrix((data, (row_idx, col_idx)), shape = (n_gnl, n_gene), dtype = np.float32)
    adata = anndata.AnnData(X = X, obs = granules.copy())
    adata.obs["granule_id"] = [f"gnl_{i}" for i in range(n_gnl)]
    adata.obs = adata.obs.astype({"granule_id": str})
    adata.var["genes"] = genes
    adata.var_names = genes
    adata.var_keys = genes
    return adata