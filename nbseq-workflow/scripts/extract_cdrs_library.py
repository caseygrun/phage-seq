import nbseq
import pandas as pd

from common import snakemake_log

with snakemake_log(snakemake) as logfile:

	# read dataframe, index should be ASV name/sequence, column should be
	# translated, aligned AA sequence (we ignore the name here and set it to the index)
	df = pd.read_csv(snakemake.input[0], index_col=0)

	library = snakemake.wildcards.library
	seqs = df.iloc[:,0]
	print(f"Extracting CDRs for {len(seqs)} ASVs from {library}")

	cdr_positions = snakemake.params.CDRs
	print("CDR positions:")
	print(cdr_positions)

	cdr_names = sorted(list(cdr_positions.keys()))

	cdrs = nbseq.extract_CDRs(seqs, library=library, CDRs=cdr_positions)
	cdrs['library'] = library

	cdrs = nbseq.identify_clonotype(cdrs,id_col='clonotypeID',seq_col='clonotype')
	cdrs = nbseq.identify_clonotype(cdrs,cdrs=['CDR3'],id_col='CDR3ID',seq_col=None)
	cdrs = cdrs[['clonotypeID', 'CDR3ID', 'clonotype'] + cdr_names + ['library']]

	# output with index which should be the aaSVID
	# pd.concat(cdrs).to_csv(snakemake.output[0], index=True)
	cdrs.to_csv(snakemake.output[0], index=True)
