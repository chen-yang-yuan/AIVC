#!/bin/bash
#SBATCH --job-name=ncd_ratio_2_selection
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

# NCD ratio analysis, step 2: selection statistics and candidate gene lists (2_selection.py;
# plans/ncd_ratio_exploration.md section 4). About 4 min locally in AIVC-env; this wrapper exists for HGCC.
# Needs: the compartment matrices, cell_ids.npy and gene_ids.npy of the five samples (data/<sample>/processed_data/),
#   output/1_exploration/ncd_ratio/1_filter/filtered_genes.npy (step 1), 0_counts/gene_table.csv (step 0, for a check),
#   data/_utils/hallmark_pathways_filtered.gmt, output/1_exploration/enrichment/overlap_genes_3plus_samples.npy,
#   data/_utils/gene_biotype.csv (protein-coding genes; written once by gene_biotype.py, needs the network),
#   code: 2_selection.py, ncd_ratio_utils.py, ../subtype_utils.py
# Submit from this folder:  sbatch 2_selection.sh      -> output/1_exploration/ncd_ratio/2_selection/ (emptied first)

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

python3 2_selection.py "$@"

echo "Job finished at $(date)"
