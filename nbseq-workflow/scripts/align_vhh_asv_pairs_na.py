from nbseq.asvs.align import align_local_bowtie
from common import snakemake_log

with snakemake_log(snakemake) as logfile:
	align_local_bowtie(
		input_fasta_fwd = snakemake.input['fwd'],
		input_fasta_rev = snakemake.input['rev'],
		input_index = snakemake.input['index'],
		library = snakemake.wildcards['library'],
		output_bam_coord_sorted = snakemake.output['coord_sorted'],
		output_bam_name_sorted = snakemake.output['name_sorted'],
		log = snakemake.log[0],
		n_penalty = snakemake.params['n_penalty'],
		n_ceil = snakemake.params['n_ceil'],
		threads=snakemake.threads)
