import shutil
from nbseq.asvs.align import index_reference_bowtie
from common import snakemake_log

# shutil.copy(snakemake.input[0], snakemake.output['fasta'])

with snakemake_log(snakemake) as logfile:
	index_reference_bowtie(
		reference_path = snakemake.input[0],
		output_path = snakemake.output['dir'],
		index_name=snakemake.wildcards['library'])
