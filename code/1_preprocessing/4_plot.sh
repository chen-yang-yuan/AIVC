#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --output=plot_%j.out
#SBATCH --error=plot_%j.err
#SBATCH --time=240:00:00
#SBATCH --mem=1000G
#SBATCH --cpus-per-task=16
#SBATCH --partition=nodes
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cyuan36@emory.edu

# Load Python 3.11 module (same one you used interactively)
module load python/3.11.7-gcc-13.2.0-ljcvqdx

# Print some debugging info
echo "Running on $(hostname)"
python3 --version
python3 -m site --user-site
echo "Job started at $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"

# Move to your project directory
cd ~/hulab/projects/AIVC

# Run your script
python3 4_plot.py

echo "Job finished at $(date)"