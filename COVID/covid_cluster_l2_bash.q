#!/bin/bash
#
#SBATCH --job-name=covid_l2
#SBATCH --mail-user=jf243[AT]duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=300G
#SBATCH --time=50:00:00
#SBATCH -c 11
#SBATCH --output=slurm-R-job-%J.stdout
#SBATCH --error=slurm-R-job-%J.stderr

module load R
R CMD BATCH covid_cluster_l2.R
