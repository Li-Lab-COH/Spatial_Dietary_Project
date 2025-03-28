#!/bin/bash
#SBATCH --job-name=RCTD_optimizing
#SBATCH --array=0-35
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=/home/janzules/Spatial/dietary_project/slurmOutput/RCTD_parameter_tests/RCTD_testing_%A_%a.out
#SBATCH --error=/home/janzules/Spatial/dietary_project/slurmOutput/RCTD_parameter_tests/RCTD_testing_%A_%a.err

module load R/RStudio_R-4.4.1
Rscript /home/janzules/Spatial/dietary_project/code/CellTyping/cell_typing_RCTD_loop.R $SLURM_ARRAY_TASK_ID
