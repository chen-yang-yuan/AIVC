#!/bin/bash
#SBATCH --job-name=ncd_ratio_3_audit
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

# NCD ratio analysis, step 3: audit of the candidate gene lists V-null and V-abs (3_audit.py;
# plans/ncd_ratio_exploration.md section 4.5). Runs locally in AIVC-env in about 5 min; this wrapper exists for HGCC.
# Needs: the compartment matrices, cell_ids.npy and gene_ids.npy of the five samples (data/<sample>/processed_data/),
#   output/1_exploration/ncd_ratio/1_filter/filtered_genes.npy (step 1),
#   output/1_exploration/ncd_ratio/2_selection/gene_stats.csv and gene_summary.csv (step 2),
#   data/_utils/hallmark_pathways_filtered.gmt,
#   data/_utils/gene_biotype.csv (protein-coding genes; written once by gene_biotype.py, needs the network),
#   code: 3_audit.py, ncd_ratio_utils.py, ../subtype_utils.py
# Submit from this folder:  sbatch 3_audit.sh      -> output/1_exploration/ncd_ratio/3_audit/ (emptied first)

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

python3 3_audit.py "$@"

echo "Job finished at $(date)"
