import pandas as pd

# columns = ['primers_trimmed', 'filtered', 'denoised', 'nochimeras', 'aligned', 'aas']

dfs = []
for col in snakemake.input.keys():
	df = pd.read_csv(snakemake.input[col], sep="\t")
	df.iloc[:,0] = df.iloc[:,0].str.replace('.fastq.gz','', regex=False)
	df = df.set_index(df.columns[0])
	for col in df.columns:
		if not ('percent' in col):
			df[col] = df[col].round(0).astype(int)
	dfs.append(df)

summary = pd.concat(dfs, axis=1).reset_index()
# summary.columns = columns

summary.to_csv(snakemake.output['summary'], sep="\t", index=False)
