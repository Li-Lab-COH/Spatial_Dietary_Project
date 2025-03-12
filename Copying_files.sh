#!/bin/bash
#SBATCH --job-name=data_copy
#SBATCH --output=./slurmOutput/data_transfer.log
#SBATCH --error=./slurmOutput/data_transfer.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Define source and destination paths
SOURCE="/scratch/emee/Rose_Li_VisiumHD"
DEST="/coh_labs/yunroseli/spatial_data/12_2024_samples/data/"

rsync -avzh --progress $SOURCE $DEST
