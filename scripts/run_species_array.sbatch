#!/bin/bash
#SBATCH --job-name=fish
#SBATCH --account=co_genomicdata
#SBATCH --partition=savio3_htc
#SBATCH --array=1-200
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/per_sp_%A_%a.out
#SBATCH --error=logs/per_sp_%A_%a.err

module load r/4.4.0-gcc-11.4.0

species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" data/fishes_list.txt)

scripts/run_per_species_data.sh "$species"