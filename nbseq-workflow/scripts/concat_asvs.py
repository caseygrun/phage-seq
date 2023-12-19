import pandas as pd

dfs = []
for library in snakemake.input.keys():
	runs = snakemake.input[library]
	for run_path in runs:
		df = pd.read_csv(run_path)
		df['library'] = library
		dfs.append(df)

df = pd.concat(dfs)
# duplicates may appear from analysis of two runs in parallel
df.drop_duplicates().to_csv(snakemake.output[0], index=False)
