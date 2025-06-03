#!/bin/sh

#SBATCH --job-name=codalasso_24
#SBATCH --time=00:30:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=8g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load gsl
module load Rtidyverse/4.4.0

cd $SLURM_SUBMIT_DIR

Rscript --vanilla codalasso_24.R ${SLURM_ARRAY_TASK_ID}



