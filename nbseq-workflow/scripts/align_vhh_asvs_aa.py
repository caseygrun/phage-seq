import os
import pandas as pd
import nbseq
import nbseq.utils
import nbseq.asvs.align

from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	# read dataframe, index should be ASV name/sequence, column should be
	# translated AA sequence
	input = pd.read_csv(snakemake.input['translation'])
	input = input.set_index('aaSVID')
	if 'translated' in input.columns:
		seqs = input.loc[:, 'translated'].values
	else:
		seqs = input.iloc[:, 0].values

	print(f"Sequences to align ({len(seqs)}):")
	print(seqs)

	# aligned is a dict of AA_sequence:aligned_AA_sequence
	# aligned = nbseq.align_translated_peptides(seqs,
	#     processes = snakemake.threads,
	#     library = snakemake.wildcards['library'],
	#     reference_path = os.path.abspath(snakemake.input['reference']),
	#     reference_length = snakemake.params.library['reference_length_aa'],
	#     reference_frame_start=snakemake.params.library['reference_frame_start_nt'],
	#     min_length = snakemake.params.library['min_aa_length'])

	aligned = nbseq.asvs.align.align_to_reference(seqs,
		reference_path = os.path.abspath(snakemake.input['reference']),
		library = snakemake.wildcards['library'],
		aa=True, mode='global',
		reference_length = snakemake.params.library['reference_length_aa'],
		reference_frame_start = snakemake.params.library['reference_frame_start_nt'],
		min_length = snakemake.params.library['min_aa_length'],
		threads = snakemake.threads
		# , verbose=True
		)

	# construct new DataFrame, index is ASV name/sequence, column 'aligned' is the
	# aligned AA sequence. drop any rows which didn't align or were too short (as
	# indicated by a blank `aligned`)
	# df_out = pd.DataFrame({ 'aligned': (aligned[x] for x in seqs) },
	#     index = input.index)

	df_out = pd.DataFrame({
		'aligned':aligned
	}, index=input.index)

	print("Done. Preview of output:")
	print(df_out)

	nbseq.utils.index_to_column(
		df_out.query('aligned != ""'),
		name='aaSVID').to_csv(snakemake.output[0], index=False)
