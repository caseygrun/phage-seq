from common import *
import pandas as pd
import nbseq
import nbseq.utils

with snakemake_log(snakemake) as logfile:

	min_cdr_lens = snakemake.params.library['min_CDR_length']
	min_aa_length = snakemake.params.library['min_aa_length']
	print("Minimum segment lengths:")
	print(min_cdr_lens)
	print(f"Minimum overall length: {min_aa_length}")

	cdrs = pd.read_csv(snakemake.input['cdrs'], index_col=0)
	print("Available segments:")
	print(cdrs.columns)
	print(cdrs.head())

	cdr_lens = nbseq.segment_lens(cdrs, min_cdr_lens.keys())

	asvs = pd.read_csv(snakemake.input['aa'], index_col=0)
	asvs['length'] = asvs['aligned'].apply(nbseq.utils.len_ungapped)
	print(asvs.head())


	df = asvs.join(cdr_lens) # pd.merge(left=cdr_lens, right=asvs, on='aaSVID')

	# keep = df.loc[(
	# 	(df['length'] >= snakemake.params.library['min_aa_length']) &
	# 	(df['CDR1'] >= snakemake.params.library['min_CDR_length']['CDR1']) &
	# 	(df['CDR2'] >= snakemake.params.library['min_CDR_length']['CDR2']) &
	# 	(df['CDR3'] >= snakemake.params.library['min_CDR_length']['CDR3'])),
	# 	'aaSVID'
	# ]


	query = (f"(length >= {min_aa_length}) & " + 
		" & ".join(f"({k} >= {v})" for k,v in min_cdr_lens.items()))

	print(f"Filtering {len(cdrs)} features")

	print("Filtering on query:")
	print(query)

	keep = df.query(query).index
	print(f"Kept {len(keep)} features.")

	# keep = df.query(
	# 	f"(length >= {snakemake.params.library['min_aa_length']}) & "
	# 	f"(CDR1   >= {snakemake.params.library['min_CDR_length']['CDR1']}) & "
	# 	f"(CDR2   >= {snakemake.params.library['min_CDR_length']['CDR2']}) & "
	# 	f"(CDR3   >= {snakemake.params.library['min_CDR_length']['CDR3']}) "
	# ).loc[:,'aaSVID']

	# df_filtered = df.query(
	# 	'(library == "alpaca" & length > 75 & CDR1 > 4 & CDR2 > 6 & CDR3 > 15) |' +
	# 	'(library == "synthetic" & length > 85 & CDR1 > 4 & CDR2 > 4 & CDR3 > 8) ')

	asvs.loc[asvs.index.isin(keep),:].to_csv(snakemake.output['aa'],index=True)
	cdrs.loc[cdrs.index.isin(keep),:].to_csv(snakemake.output['cdrs'],index=True)

	# asvs.set_index('aaSVID').loc[keep,:].to_csv(snakemake.output['aa'],index=True)
	# cdrs.set_index('aaSVID').loc[keep,:].reset_index().to_csv(snakemake.output['aa'],index=True)
