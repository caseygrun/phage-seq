import pandas as pd
input = pd.read_csv(snakemake.input[0], sep="\t")
input.to_csv(snakemake.output[0], sep=",", index=False)
