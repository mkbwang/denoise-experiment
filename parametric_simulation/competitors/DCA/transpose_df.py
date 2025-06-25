import pandas as pd

if __name__ == "__main__":

    counts_df = pd.read_csv("input.tsv", index_col=0, sep='\t')
    counts_df.T.to_csv("input_transpose.tsv", sep='\t')
    