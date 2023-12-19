	# input:
	# 	feature_table = 'intermediate/{run}/nochimeras/feature_table.biom'
	# 	aligned_asvs = expand('intermediate/{run}/merge/{library}_na.csv', library=get_libraries()) #'results/tables/vhh_asvs.csv'
	# output:
	# 	feature_table = 'intermediate/{run}/align/feature_table.biom',
	# 	summary = 'intermediate/{run}/align/summary.txt'

import pandas as pd
import nbseq.ft

from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	print("Reading feature table...")
	ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])

	print("Reading aligned ASVs...")
	aligned_asvs = pd.concat([pd.read_csv(f) for f in snakemake.input['aligned_asvs']])
	
	print(f"Feature table before filtering: {repr(ft)} \n = {ft.sum(axis='whole')} reads")
	print(ft.head())

	# remove any rows corresponding to reads that didn't map to reference
	ft = ft.filter(aligned_asvs.pair_id, axis='observation')
	print(f"Feature table after filtering: {repr(ft)} \n = {ft.sum(axis='whole')} reads")
	print(ft.head())

	# re-label the features according to the MD5 hash of the sequence (ASVID)
	# instead of MD5 hash of the concatenated read pair (paid_id)
	mapping = aligned_asvs.set_index('pair_id')['ASVID']
	print("Mapping: ")
	print(mapping)
	print(ft.ids('observation'))

	ft = nbseq.ft.group_feature_table(ft, mapping, axis='observation')
	print(f"Grouped feature table: {repr(ft)}")

	# this can be handled by summarize_feature_table
	if 'summary' in snakemake.output:
		# summarize how many reads we lost in the process
		summary = nbseq.ft.summarize(ft, summary_name='aligned')
		summary.to_csv(snakemake.output['summary'], index=False, sep="\t")

	# write updated feature table
	nbseq.ft.write_feature_table_biom(ft, snakemake.output['feature_table'])
