import nbseq
import pandas as pd

from common import snakemake_log

with snakemake_log(snakemake) as logfile:

	# read dataframe, index should be ASV name/sequence, column should be
	# translated, aligned AA sequence (we ignore the name here and set it to the index)
	df = pd.read_csv(snakemake.input[0], index_col=0)

	cdr_dfs = []

	# for rule extract_cdrs:
	if 'library_CDRs' in snakemake.params.keys():
		library_CDRs = snakemake.params['library_CDRs']
	# for rule extract_library_cdrs:
	elif 'library' in snakemake.wildcards.keys():
		library = snakemake.wildcards.library
		library_CDRs = { library : snakemake.params['CDRs'] }

	for library, cdr_positions in library_CDRs.items():
		if 'library' in df.columns:
			seqs = df.iloc[(df['library'] == library).values,0]#, :].iloc[:,0]
		else:
			seqs = df.iloc[:,0]

		print(f"Extracting CDRs for {len(seqs)} ASVs from {library}")
		print("CDR positions:")
		print(cdr_positions)
		cdr_names = sorted(list(cdr_positions.keys()))

		cdrs = nbseq.extract_CDRs(seqs, library=library, CDRs=cdr_positions)
		cdrs['library'] = library

		cdrs = nbseq.identify_clonotype(cdrs,id_col='clonotypeID',seq_col='clonotype')
		cdrs = nbseq.identify_clonotype(cdrs,cdrs=['CDR3'],id_col='CDR3ID',seq_col=None)

		column_order = ['clonotypeID', 'CDR3ID', 'clonotype'] + cdr_names + ['library']
		columns = [c for c in column_order if c in cdrs.columns]
		print("Output columns:")
		print(columns)
		cdrs = cdrs[columns]
		cdr_dfs.append(cdrs)

	# output with index which should be the aaSVID
	pd.concat(cdr_dfs).to_csv(snakemake.output[0], index=True)
