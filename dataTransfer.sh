#!/bin/bash
#SBATCH --job-name=data_transfer
#SBATCH --output=./slurmOutput/data_transfer.log
#SBATCH --error=./slurmOutput/data_transfer.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Define source and destination paths
SOURCE="janzules@gemini-data1.coh.org:/scratch/emee/Rose_Li_VisiumHD"
DEST="/home/janzules/Spatial/Li_lab_samples_12-24/data/"

rsync -avzh --progress $SOURCE $DEST
