#!/bin/sh

#SBATCH --job-name=AutoClass_ZINB
#SBATCH --time=00:30:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=10g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=3
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out

source /home/wangmk/.bashrc
micromamba activate Autoclass

cd $SLURM_SUBMIT_DIR


python run_AutoClass.py -n 50 -d 10 -s ${SLURM_ARRAY_TASK_ID} -m zinb
python run_AutoClass.py -n 50 -d 5 -s ${SLURM_ARRAY_TASK_ID} -m zinb
python run_AutoClass.py -n 50 -d 2 -s ${SLURM_ARRAY_TASK_ID} -m zinb
python run_AutoClass.py -n 50 -d 1 -s ${SLURM_ARRAY_TASK_ID} -m zinb

