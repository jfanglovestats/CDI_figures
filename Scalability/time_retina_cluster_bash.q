#!/bin/bash
#
#SBATCH --job-name=time_retina
#SBATCH --mail-user=jf243[AT]duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=300G
#SBATCH --time=100:00:00
#SBATCH -c 1
#SBATCH --output=slurm-R-job-%J.stdout
#SBATCH --error=slurm-R-job-%J.stderr

module load R
R CMD BATCH time_retina_cluster.R