import matplotlib.colors as clr
import matplotlib.pyplot as plt
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

color_cts = clr.LinearSegmentedColormap.from_list("magma", ["#000003", "#3B0F6F", "#8C2980", "#F66E5B", "#FD9F6C", "#FBFCBF"], N=256)

data_dir = "../../data/merged_data/"
output_dir = "../../output/merged_data/"

# Read data
adata = sc.read_h5ad(data_dir + "adata_embedded.h5ad")

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
    plt.savefig(output_dir + f"{key}_umap.jpeg", dpi = 300, bbox_inches = "tight")
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
    plt.savefig(output_dir + f"{key}_tsne.jpeg", dpi = 300, bbox_inches = "tight")
    plt.close()

for key in ["EPCAM", "KRT20"]:
    
    sc.set_figure_params(figsize = (12, 12))
    ax = sc.pl.umap(adata, alpha=1, color=key, cmap=color_cts, size=0.5, title=" ", show=False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.savefig(output_dir + f"{key}_umap.jpeg", dpi = 300, bbox_inches = "tight")
    plt.close()
    
    sc.set_figure_params(figsize = (12, 12))
    ax = sc.pl.tsne(adata, alpha=1, color=key, cmap=color_cts, size=0.5, title=" ", show=False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.savefig(output_dir + f"{key}_tsne.jpeg", dpi = 300, bbox_inches = "tight")
    plt.close()

print("Plotting done.")