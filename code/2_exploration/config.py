"""Constants for the Fig. 1 exploratory analysis (nuclear vs cytoplasmic expression in tumor cells).

Single source of truth for dataset names, paths, thresholds and gene lists. Paths are relative to this folder
(code/2_exploration), matching the repo rule that scripts run with their own directory as the working directory.
"""

# ---------------------------------------------------------------- datasets
DATASETS = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC",
            "Xenium_5K_LC", "Xenium_5K_Prostate", "Xenium_5K_Skin"]

# The shared `settings` dict convention used across code/; kept for continuity with other notebooks.
settings = {ds: {} for ds in DATASETS}

CANCER_TYPE = {"Xenium_5K_BC": "Breast", "Xenium_5K_OC": "Ovarian", "Xenium_5K_CC": "Cervical",
               "Xenium_5K_LC": "Lung", "Xenium_5K_Prostate": "Prostate", "Xenium_5K_Skin": "Melanoma"}


def short(ds):
    return ds.replace("Xenium_5K_", "")


DATA_DIR = "../../data/"
UTILS_DIR = "../../data/_utils/"
OUT_DIR = "../../output/2_exploration/"

# ---------------------------------------------------------------- cell types
# Hard-coded union: each dataset carries only 9-16 of these categories, so never derive the vocabulary from one
# dataset's categories.
CELL_TYPES_18 = ["Adipocyte", "B cell", "CD4+ T cell", "CD8+ T cell",
                 "Dendritic cell", "Endothelial cell",
                 "Epithelial cell (non-malignant)", "Fibroblast (CAF)",
                 "Lymphatic endothelial cell", "Malignant cell", "Mast cell",
                 "Mesothelial cell", "Myeloid cell", "Pericyte",
                 "Smooth muscle cell", "T cell", "Mixed", "Unknown"]
TUMOR_TYPE = "Malignant cell"
IMMUNE_TYPES = ["B cell", "CD4+ T cell", "CD8+ T cell", "T cell",
                "Dendritic cell", "Myeloid cell", "Mast cell"]
STROMAL_TYPES = ["Fibroblast (CAF)", "Endothelial cell", "Lymphatic endothelial cell",
                 "Pericyte", "Smooth muscle cell", "Adipocyte", "Mesothelial cell"]

# ---------------------------------------------------------------- positive controls
# MALAT1 / NEAT1 / XIST are not on the Xenium 5K panel. These panel lncRNAs are nuclear-retained in the data.
POS_CTRL_NUCLEAR = ["MEG3", "MIAT", "PVT1", "CRNDE", "HOTAIR"]
POS_CTRL_CYTO = ["NORAD"]
# Candidate cytoplasm-dominated mRNAs; intersected with the panel at runtime (most are absent from the 5K panel).
POS_CTRL_CYTO_CANDIDATES = ["KRT8", "KRT18", "KRT19", "EPCAM", "TPT1", "EEF1A1", "FTH1", "FTL", "B2M",
                            "TMSB4X", "ACTB", "GAPDH", "HSPA8", "EEF1G", "PKM", "LDHA"]
N_TOP_ABUNDANT_CTRL = 10   # data-driven cytoplasmic controls: the most abundant genes

# ---------------------------------------------------------------- thresholds
D_MIN = 30              # minimum matched depth (min of nuclear, cytoplasmic) for rarefied analyses
REQUIRE_ONE_NUCLEUS = True
N_MIN_LOC = 10          # per-unit trials n_ig required for the localisation reliability
MIN_UNITS = 200         # units needed before a per-gene reliability is reported
REL_MIN = 0.1           # within-compartment reliability needed for the disattenuated concordance
MEAN_MIN_CELL = 0.5     # cell-level per-gene statistics: mean in-cell count per cell
MEAN_MIN_TILE = 0.05    # tile-level per-gene statistics
MEAN_MIN_ASSOC = 0.2    # genes entering the association GLMs
MIN_TOTAL_COUNTS_GENE = 200   # pooled transcripts needed for a per-gene log-OR class
LOG2 = 0.6931471805599453     # |beta| threshold for nuclear-retained / cytoplasm-enriched classes
TILE_TARGET_CELLS = 15
TILE_MIN_CELLS = 8
TILE_CANDIDATES_UM = (30, 40, 50, 60, 75, 90, 110, 130, 160, 200, 250)
OFFSET_CLIP = 6.0
N_EXPR_BINS = 20
CTRL_PERCENTILE = 90    # candidates must reach this percentile within their abundance bin
N_MATCHED_SETS = 50
N_BOOT = 200
B_NULL = 5
N_PERM = 2
SEED = 0

# clustering
DETECT_FRAC_MIN = 0.005   # shared gene set: detected in >= this fraction of cells in either compartment
PCA_COMPS = 30
KNN = 15
RESOLUTIONS = (0.5, 1.0)
PRIMARY_RES = 0.5
JACCARD_CYTO_ONLY = 0.3
JACCARD_REPRO = 0.5
PURITY_K = 10
CLUSTER_MATRICES = ["nuc_m", "cyto_m", "total_m", "total_full", "nuc_h1", "nuc_h2", "cyto_h1", "cyto_h2"]

# DE
LOGFC_DE = 1.0
LOGFC_NULL = 0.5
TOPK_DE = 50

# neighbourhood
SIGMA = 30.0
KERNEL_CUT = 90.0
RADII_COUNT = (50.0, 90.0)
N_NICHES = 8
N_BINS = 5

# pathways
GMT_FILE = "all_pathways_filtered.gmt"
SG_THR = 0.4
PATHWAY_MIN_GENES = 10

# quick (development) mode overrides
QUICK_N_CELLS = 20_000
QUICK_N_BOOT = 50
QUICK_B_NULL = 3
QUICK_N_MATCHED_SETS = 20

# candidate calls (step 9)
ASSOC_MIN_EXCESS = 0.005     # deviance explained beyond the permutation baseline
R_TRUE_DIVERGENT = 0.5       # disattenuated nuclear-cytoplasmic correlation below this = divergent
REL_CANDIDATE = 0.1          # localisation reliability required for a per-gene target
CANDIDATE_AXES = ["subtype_nuc", "niche", "tumor_frac_bin", "immune_frac_bin", "morph_ratio_bin"]
MIN_SAMPLES_RECURRENT = 3

# step 3b: cytoplasm beyond nucleus (residual analysis); pre-registered thresholds
RESID_N_PCS_NUC = 30         # nuclear PCs offered to the design
RESID_PC_REL_MIN = 0.2       # nuclear PCs kept in the design need this split-half reliability
RESID_SPLINE_KNOTS = 3       # cubic spline of log depth: n_knots + degree - 1 = 5 columns
RESID_N_PCS = 30             # residual PCs
RESID_STRATA_MIN = 20        # cells per stratum (nuclear cluster x depth quintile x segmentation method)
RESID_MIN_CLUSTER = 200
RESID_MIN_CLUSTER_FRAC = 0.01
RESID_REPRO_F1 = 0.6
RESID_PURITY_EXCESS = 0.10
RESID_EFFECT_SD = 0.3
RESID_MIN_CYTO_ONLY_DE = 5
RESID_DEPTH_RATIO_MIN = 0.5
RESID_SEG_SHARE_MAX = 0.9
RESID_MORPH_SD_MAX = 0.5
RESID_PRIMARY_FEATURE = "tumor_frac"          # kept for the score dose-response tables
RESID_INTERFACE_FEATURES = ["tumor_frac", "immune_kernel"]   # either may satisfy the interface criterion (boundary OR immune-rich)
RESID_FEATURES = ["tumor_frac", "immune_kernel", "log1p_dist_nontumor", "log1p_crowding", "log_cell_area", "nuc_cell_ratio"]
RESID_SCORE_REL_MIN = 0.3
RESID_MORAN_MIN = 0.05
RESID_N_SCORE_PCS = 3
N_PERM_STRAT = 200
QUICK_N_PERM_STRAT = 50
RESID_ZOOM_UM = 600

# step 0b: annotation QC (report only; the tumor definition stays cell_type_merged == "Malignant cell")
MARKER_SETS = {
    "epithelial": ["EPCAM", "CDH1", "MUC1"],
    "basal": ["TP63"],
    "proliferation": ["MKI67", "TOP2A", "PCNA"],
    "immune": ["CD3E", "CD68"],
    "endothelial": ["PECAM1"],
    "oncogenic": ["MYC", "CCND1", "EGFR", "MET", "CDKN2A", "MDM2", "TERT"],
}
LINEAGE_MARKERS = {
    "Xenium_5K_BC": ["ESR1", "PGR", "ERBB2", "GATA3", "FOXA1"],
    "Xenium_5K_OC": ["PAX8", "WT1", "MUC16", "MSLN"],
    "Xenium_5K_CC": ["TP63", "SOX2", "CDKN2A"],
    "Xenium_5K_LC": ["NKX2-1"],
    "Xenium_5K_Prostate": ["AR", "NKX3-1", "AMACR", "TMPRSS2", "ACPP"],
    "Xenium_5K_Skin": ["SOX10", "MLANA", "MITF", "DCT"],
}
EPITHELIAL_TUMOR = {"Xenium_5K_BC": True, "Xenium_5K_OC": True, "Xenium_5K_CC": True, "Xenium_5K_LC": True,
                    "Xenium_5K_Prostate": True, "Xenium_5K_Skin": False}   # melanoma is EPCAM-negative
MANUAL_ANNOTATION = ["Xenium_5K_LC", "Xenium_5K_Prostate", "Xenium_5K_Skin"]
TRANSFER_TRAIN = ["Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC"]
TRANSFER_EXCLUDE = ["Mixed", "Unknown"]
TRANSFER_EPITHELIAL_CLASS = "Epithelial (any)"   # Malignant + non-malignant epithelium merged for training
TRANSFER_MAX_PER_CLASS_DS = 5000
TRANSFER_CONF = 0.8
