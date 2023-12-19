import pandas as pd

dfs = []
for f in snakemake.input:
	df = pd.read_csv(f, sep="\t")
	if 'row' in snakemake.params.keys():
		df.rename(columns={df.columns[0]: snakemake.params['row']}, inplace=True)
	dfs.append(df)

df = pd.concat(dfs)
df.to_csv(snakemake.output[0], sep="\t", index=False)
