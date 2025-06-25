import pandas as pd
import argparse
import os
import pickle
import magic
import scprep

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

    # filename = "s3://sign-series/Simulation/HMP16S/NB_sim.txt"
    counts_df = pd.read_csv(input_file, index_col=0)
    counts_df = counts_df.transpose()

    # normalize by median library size
    # then take square root
    normalized_counts_df = scprep.normalize.library_size_normalize(counts_df)
    normalized_counts_df = scprep.transform.sqrt(normalized_counts_df)

    # counts_mat = counts_df.to_numpy()
    magic_operator = magic.MAGIC()

    denoised_counts_df = magic_operator.fit_transform(normalized_counts_df)
    denoised_counts_df = denoised_counts_df ** 2
    denoised_counts_df = denoised_counts_df.round(6)

    output_folder = "/nfs/turbo/sph-ligen/wangmk/denoise_experiment/parametric_simulation/competitors/MAGIC"
    output_file = os.path.join(output_folder, args.model,
                        f"denoised_count_{args.model}_n{args.number}_d{args.dispersion}_{args.seed}.csv")

    denoised_counts_df.to_csv(output_file)

