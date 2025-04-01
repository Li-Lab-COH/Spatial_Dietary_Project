#!/bin/bash
#SBATCH --job-name=RCTD_full
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=70G
#SBATCH --time=24:00:00
#SBATCH --output=/home/janzules/Spatial/dietary_project/slurmOutput/RCTD_full_annotation_tests/test3/RCTD_%a.out
#SBATCH --error=/home/janzules/Spatial/dietary_project/slurmOutput/RCTD_full_annotation_tests/test3/RCTD_%a.err

module load R/RStudio_R-4.4.1
Rscript /home/janzules/Spatial/dietary_project/code/CellTyping/RCTD_full_reference_full_dataset_perSample.R $SLURM_ARRAY_TASK_ID

