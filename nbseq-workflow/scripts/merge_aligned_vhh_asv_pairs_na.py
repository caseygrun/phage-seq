import pandas as pd
import nbseq.asvs.align
import nbseq.utils

from common import *
with snakemake_log(snakemake) as logfile:

	alignment_table = nbseq.asvs.align.export_aligned_reads_paired(
		samfile_path = snakemake.input['alignment'],
		reference_path = snakemake.input['reference'],
		reference_frame_start = snakemake.params.library['reference_frame_start_nt'],
		min_fwd_end           = snakemake.params.library['min_fwd_end'],
		max_rev_start         = snakemake.params.library['max_rev_start'],
		verbose = True)

	alignment_table['ASVID'] = alignment_table['seq'].apply(nbseq.utils.md5)
	alignment_table = nbseq.utils.index_to_column(alignment_table,'pair_id')

	alignment_table[['ASVID','seq','fwd','int','rev','pair_id']].to_csv(snakemake.output[0],index=False)

	# snakemake_log_end()
