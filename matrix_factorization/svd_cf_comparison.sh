#!/bin/sh

#SBATCH --job-name=svd_cf_comparison
#SBATCH --time=01:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=10g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load gsl
module load Rtidyverse/4.4.0

cd $SLURM_SUBMIT_DIR

Rscript --vanilla logisticSVD_logisticCF_comparison.R ${SLURM_ARRAY_TASK_ID}



