import pandas as pd
from nbseq.utils import dataframe_to_fasta, index_to_column

uniqs = pd.read_csv(snakemake.input['csv'], index_col=0)
uniqs = index_to_column(uniqs, 'pair_id')

for direction in ['fwd','rev']:
	dataframe_to_fasta(uniqs,
		snakemake.output[direction], name_col='pair_id', seq_col=direction)
