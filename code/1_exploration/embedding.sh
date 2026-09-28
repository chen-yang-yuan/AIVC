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

# Stagger the starts so 24 tasks do not import the same conda env in the same second (BeeGFS metadata flakes).
sleep $(( (SLURM_ARRAY_TASK_ID % 8) * 15 ))

# Import check with retries: transient "No such file or directory: ...egg-info/PKG-INFO" errors from
# pkg_resources while importing scanpy happen under load; a genuine env problem still fails after 5 tries.
for attempt in 1 2 3 4 5; do
    if python3 -c "import scanpy" 2>/dev/null; then break; fi
    echo "import scanpy failed (attempt $attempt), retrying in $((attempt * 30))s"
    sleep $((attempt * 30))
    if [ "$attempt" -eq 5 ]; then echo "import scanpy failed 5 times, giving up"; exit 1; fi
done

python3 embedding.py --task "$SLURM_ARRAY_TASK_ID"

echo "Job finished at $(date)"
