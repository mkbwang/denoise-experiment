from AutoClass import AutoClassImpute
import argparse
import os
import pandas as pd
import numpy as np
import time

if __name__ == "__main__":

    input_folder = "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/data"
    # parse arguments
    parser = argparse.ArgumentParser(description="Input arguments for MAGIC")

    # add arguments
    parser.add_argument('-n', '--number', type=int, default=50, help="Number of samples for each phenotype [default=50]")
    parser.add_argument('-d', '--dispersion', type=int, default=5, help="Dispersion size parameter for negative binomial distribution [default=5]")
    parser.add_argument('-s', '--seed', type=int, default=1, help='Simulation seed for replication [default=1]')
    parser.add_argument('-m', '--model', type=str, default="zinb", help="Simulation distribution [default=zinb]")

    # Parse the arguments
    args = parser.parse_args()

    input_file = os.path.join(input_folder, args.model, 
        f"sim_count_{args.model}_n{args.number}_d{args.dispersion}_{args.seed}.csv")
    
    count_df = pd.read_csv(input_file, index_col=0)
    count_mat = count_df.to_numpy().T

    smoothed_result = AutoClassImpute(data=count_mat, num_cluster = [4,5,6])
    imputed_matrix = smoothed_result["imp"]

    imputed_matrix = 2**imputed_matrix - 1
    imputed_matrix = np.round(imputed_matrix, decimals=5)
    imputed_matrix = pd.DataFrame(imputed_matrix.T, index=count_df.index, columns=count_df.columns)

    output_folder = "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/AutoClass"
    output_file = os.path.join(output_folder, args.model,
                        f"denoised_count_{args.model}_n{args.number}_d{args.dispersion}_{args.seed}.csv")

    imputed_matrix.to_csv(output_file)
