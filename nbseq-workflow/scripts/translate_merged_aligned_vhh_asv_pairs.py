import pandas as pd
import nbseq.asvs.align
import nbseq.utils

from common import *
with snakemake_log(snakemake) as logfile:
# with open('.test','w') as logfile:

	# input = pd.read_csv(snakemake.input['na']).set_index('ASVID')
	#
	# translated = (nbseq.asvs.align.translate_aligned(input['seq'],
	# 		reference_path=snakemake.input['reference'],
	# 		reference_frame_start=snakemake.params.library['reference_frame_start_nt'],
	# 		reference_length_aa=snakemake.params.library['reference_length_aa'],
	# 		suppress_amber=True
	# 		,verbose=True
	# 		)
	# 	.to_frame(name='translated'))
	#
	# translated['aaSVID'] = translated['translated'].apply(nbseq.utils.md5_ungapped)
	# translated = nbseq.utils.index_to_column(translated,'ASVID')
	#
	# # maps nucleic acid ASV IDs to amino acid ASV IDs
	# translated[['ASVID','aaSVID']].to_csv(snakemake.output['na_to_aa'], index=False)
	#
	# # amino acid ASV IDs to translations
	# translated[['aaSVID','translated']].to_csv(snakemake.output['aa'], index=False)

	# snakemake_log_end()


	chunksize = 10 ** 6
	print(f"Using chunksize = {chunksize} rows")
	print(f"Using parallel processing with {snakemake.threads} threads")
	with pd.read_csv(snakemake.input['na'], chunksize=chunksize) as reader:
		first = True
		nrows = 0
		for i, input in enumerate(reader):
			print(f"Reading chunk {i} = {len(input)} rows")

			input = input.set_index('ASVID')

			translated = (nbseq.asvs.align.translate_aligned(input['seq'],
					reference_path=snakemake.input['reference'],
					reference_frame_start=snakemake.params.library['reference_frame_start_nt'],
					reference_length_aa=snakemake.params.library['reference_length_aa'],
					suppress_amber=True,
					threads=snakemake.threads
					# ,verbose=True
					)
				.to_frame(name='translated'))

			translated['aaSVID'] = translated['translated'].apply(nbseq.utils.md5_ungapped)
			translated = nbseq.utils.index_to_column(translated,'ASVID', inplace=True)

			# maps nucleic acid ASV IDs to amino acid ASV IDs
			translated[['ASVID','aaSVID']].to_csv(snakemake.output['na_to_aa'],
				mode='w' if first else 'a', header=first,
				index=False)

			# amino acid ASV IDs to translations
			translated[['aaSVID','translated']].to_csv(snakemake.output['aa'],
				mode='w' if first else 'a', header=first,
				index=False)


			nrows += len(input)
			print(f"Wrote chunk {i}; {nrows} total written")
			first = False

	print("re-reading aa sequences to remove duplicates")
	translated = pd.read_csv(snakemake.output['aa'])

	# drop rows that are entirely duplicated (e.g. same aaSVID, same aa seq);
	# can be due to multiple degenerate NA sequences producing the same AA
	# sequence
	translated = translated.drop_duplicates()


	# find rows where the same aaSVID mapped to multiple translations. Since
	# aaSVID is the hash of the ungapped sequence, these rows differ only in
	# how the nucleic acid sequence was aligned. We will be re-aligning the
	# AA sequences anyway, so arbitrarily choose the first of such sequences
	# here
	translated.groupby('aaSVID')['translated'].first().reset_index()

	translated.to_csv(snakemake.output['aa'], index=False)
