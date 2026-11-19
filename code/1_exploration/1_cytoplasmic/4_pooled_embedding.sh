#!/bin/bash
#SBATCH --job-name=pooled_embedding
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# Pooled embedding of the tumor cells of the five Xenium 5K samples (4_pooled_embedding.py), the pooled counterpart
# of 2_embedding.sh: one job, ~624k cells x 634 genes. Needs the same inputs as 2_embedding.sh (data of the five
# samples on HGCC, gene set from 'make enrich-push').
# Submit from this folder: first argument = correction mode (default none), second = compartment (default cytoplasmic):
#   sbatch 4_pooled_embedding.sh                   # uncorrected pool                 -> output/1_exploration/pooled/
#   sbatch 4_pooled_embedding.sh sample            # per-sample transform, pooled PCA -> output/1_exploration/pooled_sample_corrected/
#   sbatch 4_pooled_embedding.sh sample nuclear    # the same on the nuclear compartment -> output/1_exploration/pooled_sample_corrected_nuclear/
# Afterwards, locally:  make out-pull   (pulls <folder>/adata_pooled_<compartment>_embedded.h5ad)
# The output is skipped if it exists; add --force to the python call below to recompute, --no-tsne to skip t-SNE.

set -euo pipefail
CORRECTION=${1:-none}
COMPARTMENT=${2:-cytoplasmic}

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate preprocessing-env

cd ~/hulab/projects/AIVC/code/1_exploration/1_cytoplasmic
mkdir -p logs

echo "Host: $(hostname)"
echo "Job:  $SLURM_JOB_ID"
echo "PWD:  $(pwd)"
echo "Correction: $CORRECTION"
echo "Compartment: $COMPARTMENT"
which python
python --version

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

GENE_SET=~/hulab/projects/AIVC/output/1_exploration/enrichment/overlap_genes_3plus_samples.npy
if [ ! -f "$GENE_SET" ]; then
    echo "Missing $GENE_SET: run 'make enrich-push' locally before submitting"; exit 1
fi

# Import check with retries (transient pkg_resources errors while importing scanpy under load, see 2_embedding.sh)
for attempt in 1 2 3 4 5; do
    if python3 -c "import scanpy" 2>/dev/null; then break; fi
    echo "import scanpy failed (attempt $attempt), retrying in $((attempt * 30))s"
    sleep $((attempt * 30))
    if [ "$attempt" -eq 5 ]; then echo "import scanpy failed 5 times, giving up"; exit 1; fi
done

python3 4_pooled_embedding.py --correction "$CORRECTION" --compartment "$COMPARTMENT"

echo "Job finished at $(date)"
