import matplotlib.pyplot as plt
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# Read data
adata = sc.read_h5ad("adata_embedded.h5ad")

# Plot UMAP and t-SNE
for key in ["batch", "cell_type_merged"]:
    
    sc.set_figure_params(figsize = (12, 12))
    ax = sc.pl.umap(adata, alpha=1, color=key, size=0.5, title=" ", show=False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.savefig(f"{key}_umap.jpeg", dpi = 300, bbox_inches = "tight")
    plt.close()
    
    sc.set_figure_params(figsize = (12, 12))
    ax = sc.pl.tsne(adata, alpha=1, color=key, size=0.5, title=" ", show=False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.savefig(f"{key}_tsne.jpeg", dpi = 300, bbox_inches = "tight")
    plt.close()

print("Plotting done.")