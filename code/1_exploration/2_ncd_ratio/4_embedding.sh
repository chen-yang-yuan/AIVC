#!/bin/bash
#SBATCH --job-name=ncd_ratio_4_embedding
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-23
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# NCD ratio analysis, step 4: pooled embedding and clustering (4_embedding.py; plans/ncd_ratio_exploration.md
# section 5). SLURM array, one task per (gene list, correction, run): 3 x 2 x 4 = 24 tasks, about 596,000 cells each.
# Task ids: `python3 4_embedding.py --status` prints the table (0-7 primary, 8-15 larger, 16-23 union; within a list
# zscore then depth; within a correction data, null, half_a, half_b).
#
# Needs on HGCC:
#   data/<sample>/processed_data/{cytoplasmic,nuclear}_expression_matrix.npz, cell_ids.npy, gene_ids.npy and
#   data/<sample>/intermediate_data/adata.h5ad of the five samples (already there from the earlier jobs),
#   data/_utils/gene_biotype.csv (comes with git),
#   output/1_exploration/ncd_ratio/3_audit/target_genes.npy and target_genes_union.npy (`make ncd-push`, local),
#   code: 4_embedding.py, ncd_ratio_utils.py (git).
#
# Steps
#   1  local, repo root:  make push
#   2  HGCC:              cd ~/hulab/projects/AIVC && git pull
#   3  local, repo root:  make ncd-push-dry ; make ncd-push
#   4  HGCC:              cd ~/hulab/projects/AIVC/code/1_exploration/2_ncd_ratio && sbatch 4_embedding.sh
#   5  HGCC:              squeue -u $USER ; tail logs/ncd_ratio_4_embedding_<job>_<task>.out
#   6  HGCC:              python3 4_embedding.py --status        (after: conda activate preprocessing-env)
#   7  HGCC, if a task failed:  sbatch --array=<ids> 4_embedding.sh      (complete tasks are skipped)
#   8  local, repo root:  make ncd-pull-dry ; make ncd-pull      -> output/1_exploration/ncd_ratio/4_embedding/
# A complete task is skipped; add --force to the python call below to recompute.

set -euo pipefail

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate preprocessing-env

cd ~/hulab/projects/AIVC/code/1_exploration/2_ncd_ratio
mkdir -p logs

echo "Host: $(hostname)"
echo "Job:  ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "PWD:  $(pwd)"
which python
python --version

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# inputs
for f in ../../../output/1_exploration/ncd_ratio/3_audit/target_genes.npy \
         ../../../output/1_exploration/ncd_ratio/3_audit/target_genes_union.npy \
         ../../../data/_utils/gene_biotype.csv; do
    if [ ! -f "$f" ]; then
        echo "Missing $f: run 'git pull' here and 'make ncd-push' locally before submitting"; exit 1
    fi
done
for s in Xenium_5K_BC Xenium_5K_OC Xenium_5K_CC Xenium_5K_LC Xenium_5K_Prostate; do
    for f in processed_data/cytoplasmic_expression_matrix.npz processed_data/nuclear_expression_matrix.npz \
             processed_data/cell_ids.npy processed_data/gene_ids.npy intermediate_data/adata.h5ad; do
        if [ ! -f "../../../data/$s/$f" ]; then echo "Missing data/$s/$f"; exit 1; fi
    done
done

# Import check with retries (transient pkg_resources errors while importing scanpy under load, see
# ../1_cytoplasmic/4_pooled_embedding.sh); tasks start a few seconds apart for the same reason
sleep $(( (SLURM_ARRAY_TASK_ID % 8) * 5 ))
for attempt in 1 2 3 4 5; do
    if python3 -c "import scanpy" 2>/dev/null; then break; fi
    echo "import scanpy failed (attempt $attempt), retrying in $((attempt * 30))s"
    sleep $((attempt * 30))
    if [ "$attempt" -eq 5 ]; then echo "import scanpy failed 5 times, giving up"; exit 1; fi
done
python3 -c "
try:
    import louvain; print('clustering: louvain')
except Exception:
    import leidenalg; print('clustering: leiden (louvain is not installed)')
"

python3 4_embedding.py --task "$SLURM_ARRAY_TASK_ID"

echo "Job finished at $(date)"
