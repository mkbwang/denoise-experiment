#!/bin/sh

#SBATCH --job-name=HMP_ZINB_simulate
#SBATCH --time=00:20:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=10g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load Rtidyverse/4.4.0

cd $SLURM_SUBMIT_DIR


Rscript --vanilla HMP_simulate_ZINB.R -n 50 -d 10 -s ${SLURM_ARRAY_TASK_ID}
Rscript --vanilla HMP_simulate_ZINB.R -n 50 -d 5 -s ${SLURM_ARRAY_TASK_ID}
Rscript --vanilla HMP_simulate_ZINB.R -n 50 -d 2 -s ${SLURM_ARRAY_TASK_ID}
Rscript --vanilla HMP_simulate_ZINB.R -n 50 -d 1 -s ${SLURM_ARRAY_TASK_ID}

