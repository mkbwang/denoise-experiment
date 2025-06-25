#!/bin/sh

#SBATCH --job-name=DCA_ZINB
#SBATCH --time=00:30:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=10g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


source /home/wangmk/.bashrc
micromamba activate DCA

cd $SLURM_SUBMIT_DIR


input_folder="/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/data/zinb"
output_folder="/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/DCA"
input_filename_1="sim_count_zinb_n50_d10_${SLURM_ARRAY_TASK_ID}.csv"
output_filename_1="denoised_count_zinb_n50_d10_${SLURM_ARRAY_TASK_ID}.csv"
input_filename_2="sim_count_zinb_n50_d5_${SLURM_ARRAY_TASK_ID}.csv"
output_filename_2="denoised_count_zinb_n50_d5_${SLURM_ARRAY_TASK_ID}.csv"
input_filename_3="sim_count_zinb_n50_d2_${SLURM_ARRAY_TASK_ID}.csv"
output_filename_3="denoised_count_zinb_n50_d2_${SLURM_ARRAY_TASK_ID}.csv"
input_filename_4="sim_count_zinb_n50_d1_${SLURM_ARRAY_TASK_ID}.csv"
output_filename_4="denoised_count_zinb_n50_d1_${SLURM_ARRAY_TASK_ID}.csv"

temp_folder="temp_output_zinb_n50_${SLURM_ARRAY_TASK_ID}"
mkdir $temp_folder

# Run DCA
dca ${input_folder}/${input_filename_1} $temp_folder
mv ${temp_folder}/mean.tsv zinb/${output_filename_1}

dca ${input_folder}/${input_filename_2} $temp_folder
mv ${temp_folder}/mean.tsv zinb/${output_filename_2}

dca ${input_folder}/${input_filename_3} $temp_folder
mv ${temp_folder}/mean.tsv zinb/${output_filename_3}

dca ${input_folder}/${input_filename_4} $temp_folder
mv ${temp_folder}/mean.tsv zinb/${output_filename_4}



