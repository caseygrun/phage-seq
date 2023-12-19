import pandas as pd
import nbseq.ft


from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	ft = nbseq.ft.read_feature_table_biom(snakemake.input[0])
	print(f"Summarizing feature table {snakemake.input[0]}: {repr(ft)}...")
	# print(f"Writing to output files:")
	# print("\n".join("\t".join(fs) for fs in snakemake.output.items()))

	# sample features, e.g. # of features per sample
	if 'features' in snakemake.output.keys():
		print(f"Writing feature summary to {snakemake.output['features']}")
		summary_features = nbseq.ft.summarize(ft, summary_name=snakemake.params['header'], method='count')
		summary_features.to_csv(snakemake.output['features'], index=False, sep="\t")

	# total # of features
	if 'total_features' in snakemake.output.keys():
		total_features = ft.shape[0]
		(pd.DataFrame([total_features],
				columns= [snakemake.params['header']],
				index  = [snakemake.params['total_row']])
			.to_csv(snakemake.output['total_features'],
				index=True, sep="\t"))

	if 'reads' in snakemake.output.keys():
		summary = nbseq.ft.summarize(ft, summary_name=snakemake.params['header'])
		print(f"Writing read summary to {snakemake.output['reads']}")
		summary.to_csv(snakemake.output['reads'], index=False, sep="\t")
