#!/bin/bash
#SBATCH --job-name=embedding
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --array=0-23
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# One array task per (sample, compartment, n_pcs); the n_neighbors loop runs inside the task.
# Task table:  python3 embedding.py --list   (ids 0-3 BC, 4-7 OC, 8-11 CC, 12-15 LC, 16-19 Prostate, 20-23 Skin)
# Submit from this folder:
#   sbatch embedding.sh                 # all 24 tasks
#   sbatch --array=0-3 embedding.sh     # one sample (BC)
#   sbatch --array=0-23%8 embedding.sh  # at most 8 tasks at a time
# Existing outputs are skipped; add --force to the python line to recompute.

set -euo pipefail

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate preprocessing-env

cd ~/hulab/projects/AIVC/code/1_exploration
mkdir -p logs

echo "Host: $(hostname)"
echo "Job:  $SLURM_JOB_ID  task $SLURM_ARRAY_TASK_ID"
echo "PWD:  $(pwd)"
which python
python --version

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

python3 embedding.py --task "$SLURM_ARRAY_TASK_ID"

echo "Job finished at $(date)"
