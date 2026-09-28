import numpy as np
import os
import scanpy as sc
from scipy import sparse

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

sample_list = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate", "Xenium_5K_Skin"]

for sample in sample_list:
    
    # paths
    intermediate_path = f"../../data/{sample}/intermediate_data/"
    processed_path = f"../../data/{sample}/processed_data/"
    output_path = f"../../output/1_exploration/{sample}/"
    os.makedirs(output_path, exist_ok=True)
    
    # read data
    adata = sc.read_h5ad(intermediate_path + "adata.h5ad")
    cell_ids = np.load(processed_path + "cell_ids.npy", allow_pickle=True)
    assert adata.shape[0] > cell_ids.shape[0], "cell_ids should be a subset of adata.obs_names"
    adata_tumor = adata[adata.obs["cell_id"].isin(cell_ids)].copy()

    # construct data
    for label in ["nuclear", "cytoplasmic"]:
        
        compartment_expression = sparse.load_npz(processed_path + f"{label}_expression_matrix.npz")
        assert compartment_expression.shape[0] == adata_tumor.shape[0], "compartment_expression should have the same number of cells as adata_tumor"
        adata_tmp = sc.AnnData(X=compartment_expression, obs=adata_tumor.obs, var=adata_tumor.var)
    
        # embedding
        sc.pp.normalize_total(adata_tmp, target_sum = 1e4)
        sc.pp.log1p(adata_tmp)
        sc.tl.pca(adata_tmp, n_comps = 100, svd_solver = "auto")
        sc.tl.tsne(adata_tmp, n_pcs = 50)
        sc.pp.neighbors(adata_tmp, n_neighbors = 50, n_pcs = 50)
        sc.tl.umap(adata_tmp)
        adata_tmp.write(output_path + f"adata_{label}_embedded.h5ad")