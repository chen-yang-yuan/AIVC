#!/bin/bash
#SBATCH --job-name=ncd_ratio_0_counts
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

# NCD ratio analysis, step 0: counts check and choice of N (0_counts.py; plans/ncd_ratio_exploration.md).
# The script is light (about 40 s locally in AIVC-env); this wrapper exists for HGCC.
# Needs: the compartment matrices, cell_ids.npy and gene_ids.npy of the five samples (data/<sample>/processed_data/),
#   data/_utils/hallmark_pathways_filtered.gmt, output/1_exploration/enrichment/overlap_genes_3plus_samples.npy,
#   data/_utils/gene_biotype.csv (protein-coding genes; written once by gene_biotype.py, needs the network),
#   code: 0_counts.py, ncd_ratio_utils.py, ../subtype_utils.py
# Submit from this folder:  sbatch 0_counts.sh      -> output/1_exploration/ncd_ratio/0_counts/ (emptied first)

set -euo pipefail

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate preprocessing-env

cd ~/hulab/projects/AIVC/code/1_exploration/2_ncd_ratio
mkdir -p logs

echo "Host: $(hostname)"
echo "Job:  $SLURM_JOB_ID"
echo "PWD:  $(pwd)"

python3 0_counts.py "$@"

echo "Job finished at $(date)"
