import pandas as pd

cdrs = pd.read_csv(snakemake.input[0])


cols = []
if 'library' in list(cdrs.columns):
	cols = ['library']

if snakemake.wildcards['space'].lower() == 'cdr3':
	id_cols = ['CDR3ID']
	cols = ['CDR3ID','CDR3'] + cols
	cdrs = cdrs[cols]
elif snakemake.wildcards['space'].lower() == 'cdrs':
	cols = ['clonotypeID','clonotype','CDR1','CDR2','CDR3'] + cols
	id_cols = ['clonotypeID']
	cdrs = cdrs[cols]

cdrs = cdrs.drop_duplicates(subset=id_cols)
cdrs.to_csv(snakemake.output[0],index=False)
