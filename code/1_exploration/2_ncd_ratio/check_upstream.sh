#!/bin/bash
#SBATCH --job-name=ncd_ratio_check_upstream
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

# NCD ratio analysis: audit of the compartment count matrices against the raw 10x transcripts (check_upstream.py;
# plans/ncd_ratio_exploration.md section 10.8). Not a step of the analysis. This wrapper exists for HGCC.
# Needs: raw_data/transcripts.parquet, intermediate_data/adata.h5ad and the compartment matrices, cell_ids.npy,
#   gene_ids.npy of the five samples (data/<sample>/), output/1_exploration/ncd_ratio/1_filter/filtered_genes.npy,
#   output/1_exploration/ncd_ratio/3_audit/target_genes.npy,
#   data/_utils/gene_biotype.csv (protein-coding genes; written once by gene_biotype.py, needs the network),
#   code: check_upstream.py, ncd_ratio_utils.py
# Submit from this folder:  sbatch check_upstream.sh   -> output/1_exploration/ncd_ratio/upstream_check/ (emptied first)

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

python3 check_upstream.py "$@"

echo "Job finished at $(date)"
