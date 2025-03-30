#!/bin/bash
#SBATCH --job-name=RCTD
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --output=/home/janzules/Spatial/dietary_project/slurmOutput/RCTD_full_annotation_tests/test1/RCTD_%a.out
#SBATCH --error=/home/janzules/Spatial/dietary_project/slurmOutput/RCTD_full_annotation_tests/test1/RCTD_%a.err

module load R/RStudio_R-4.4.1
Rscript /home/janzules/Spatial/dietary_project/code/CellTyping/cell_typing_RCTD_loop.R $SLURM_ARRAY_TASK_ID