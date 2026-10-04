#!/bin/bash
#SBATCH --job-name=embedding
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --array=0-9
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# One array task per (sample, compartment): nuclear/cytoplasmic expression of the five tumor samples (no Skin)
# restricted to the cytoplasm-enriched HVG set from 1_enrichment.ipynb (overlap_genes_3plus_samples.npy).
# Task table:  python3 2_embedding.py --list   (ids 0-1 BC, 2-3 OC, 4-5 CC, 6-7 LC, 8-9 Prostate)
# Before submitting copy the gene set to HGCC:  make enrich-push  (from the local repo root)
# Submit from this folder:
#   sbatch 2_embedding.sh                 # all 10 tasks
#   sbatch --array=0-1 2_embedding.sh     # one sample (BC)
#   sbatch --array=0-9%5 2_embedding.sh   # at most 5 tasks at a time
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

# The gene set is produced locally (output/ is gitignored); it must be rsynced here first (make enrich-push).
GENE_SET=~/hulab/projects/AIVC/output/1_exploration/enrichment/overlap_genes_3plus_samples.npy
if [ ! -f "$GENE_SET" ]; then
    echo "Missing $GENE_SET: run 'make enrich-push' locally before submitting"; exit 1
fi

# Stagger the starts so 10 tasks do not import the same conda env in the same second (BeeGFS metadata flakes).
sleep $(( (SLURM_ARRAY_TASK_ID % 8) * 15 ))

# Import check with retries: transient "No such file or directory: ...egg-info/PKG-INFO" errors from
# pkg_resources while importing scanpy happen under load; a genuine env problem still fails after 5 tries.
for attempt in 1 2 3 4 5; do
    if python3 -c "import scanpy" 2>/dev/null; then break; fi
    echo "import scanpy failed (attempt $attempt), retrying in $((attempt * 30))s"
    sleep $((attempt * 30))
    if [ "$attempt" -eq 5 ]; then echo "import scanpy failed 5 times, giving up"; exit 1; fi
done

python3 2_embedding.py --task "$SLURM_ARRAY_TASK_ID"

echo "Job finished at $(date)"
