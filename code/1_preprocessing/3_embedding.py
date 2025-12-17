import scanpy as sc

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

data_dir = "../../data/merged_data/"

# Read raw data
adata = sc.read_h5ad(data_dir + "adata_all_raw.h5ad")

# t-SNE and UMAP embedding for granule_adata
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, n_comps = 100, svd_solver = "auto")
sc.tl.tsne(adata, n_pcs = 50)
sc.pp.neighbors(adata, n_neighbors = 50, n_pcs = 50)
sc.tl.umap(adata)
adata.write_h5ad(data_dir + "adata_embedded.h5ad")

print("Embedding done.")