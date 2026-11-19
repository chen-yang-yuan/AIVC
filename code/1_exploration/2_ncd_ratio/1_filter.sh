#!/bin/bash
#SBATCH --job-name=ncd_ratio_1_filter
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

# NCD ratio analysis, step 1: gene filter (1_filter.py; plans/ncd_ratio_exploration.md).
# The script is light (a few seconds locally in AIVC-env); this wrapper exists for HGCC.
# Needs: output/1_exploration/ncd_ratio/0_counts/gene_table.csv and by_n.csv (step 0, full run),
#   data/_utils/hallmark_pathways_filtered.gmt, output/1_exploration/enrichment/overlap_genes_3plus_samples.npy,
#   code: 1_filter.py, ncd_ratio_utils.py, ../subtype_utils.py
# Submit from this folder:  sbatch 1_filter.sh      -> output/1_exploration/ncd_ratio/1_filter/ (emptied first)

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

python3 1_filter.py "$@"

echo "Job finished at $(date)"
