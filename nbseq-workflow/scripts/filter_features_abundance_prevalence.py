from common import *
import pandas as pd
import nbseq.ft
import nbseq.utils

ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])

with snakemake_log(snakemake):
	print(f"Input feature table: {repr(ft)}")

	ids_to_keep = nbseq.ft.filter_ids_min_abundance_prevalence(ft,
		min_abundance=snakemake.params['min_abundance'],
		min_prevalence=snakemake.params['min_prevalence']
	)
	print(f"IDs to keep: {len(ids_to_keep)}")
	print(ids_to_keep)
	print()

	ftf = ft.filter(ids_to_keep, axis='observation')
	nbseq.ft.write_feature_table_biom(ftf, snakemake.output['feature_table'])

	print(f"Filtered feature table: {repr(ftf)}")
	del ft
	del ftf

	feature_data = read_delim_auto(snakemake.input['feature_data'])
	
	print("Feature data:")
	print(feature_data.head())
	idx_col = feature_data.columns[0]
	feature_data = feature_data.set_index(idx_col)
	feature_data.loc[ids_to_keep,:].to_csv(snakemake.output['feature_data'], index=True)
	del feature_data

	print("Done.")