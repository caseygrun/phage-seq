import pandas as pd
input = pd.read_csv(snakemake.input[0])

# TODO: implement clustering
pd.DataFrame({
	'aaSVID':input['aaSVID'],
	'cluster':input['aaSVID']}).to_csv(snakemake.output[0])
