#!/bin/bash
#SBATCH --job-name=copy_data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/home/janzules/Spatial/dietary_project/slurmOutput/copy_data.out
#SBATCH --error=/home/janzules/Spatial/dietary_project/slurmOutput/copy_data.err

# Do the copy
cp -r /home/rogaguilar/Spatial/data/Rose_Li_VisiumHD/Analysis/BANKSY_Normalized_QC_Filtered_minUMI_50_minGene50_MT5 \
      /home/janzules/Spatial/dietary_project/data/
