#!/bin/bash
#SBATCH --job-name=fig1_heavy
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --array=0-5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# One array task per dataset (index into config.DATASETS: 0 BC, 1 OC, 2 CC, 3 LC, 4 Prostate, 5 Skin).
# Submit from this folder after `mkdir -p logs`:  sbatch 1_heavy.sh
# A single dataset:  sbatch --array=3 1_heavy.sh

set -euo pipefail

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate preprocessing-env

cd ~/hulab/projects/AIVC/code/2_exploration
mkdir -p logs

echo "Host: $(hostname)"
echo "Job:  $SLURM_JOB_ID  task $SLURM_ARRAY_TASK_ID"
echo "PWD:  $(pwd)"
which python
python --version

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

python3 1_heavy.py --dataset "$SLURM_ARRAY_TASK_ID"

echo "Job finished at $(date)"
