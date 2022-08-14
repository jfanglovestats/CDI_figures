#!/bin/bash
#
#SBATCH --job-name=lung_all
#SBATCH --mail-user=jf243[AT]duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH -c 2
#SBATCH --output=slurm-R-job-%J.stdout
#SBATCH --error=slurm-R-job-%J.stderr

module load R
R CMD BATCH time_lung_cluster_all_labels.R
