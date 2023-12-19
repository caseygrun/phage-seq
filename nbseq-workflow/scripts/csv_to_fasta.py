import pandas as pd
from nbseq.utils import dataframe_to_fasta, index_to_column

feature_data = pd.read_csv(snakemake.input[0])
columns = list(feature_data.columns)

dataframe_to_fasta(feature_data,
	snakemake.output[0], name_col=columns[0], seq_col=columns[1])
