#!/bin/bash
#SBATCH --job-name=RCTD_testing
#SBATCH --array=0-36
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/RCTD_testing_%A_%a.out
#SBATCH --error=logs/RCTD_testing_%A_%a.err

# Optional: load R module or environment
module load R/RStudio_R-4.4.1

/home/janzules/Spatial/dietary_project/data/cell_typing_reference/rctd_parameter_grid.csv

# Optional: bind to a specific R script or command
Rscript your_rscript.R $SLURM_ARRAY_TASK_ID
