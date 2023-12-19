# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


# report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
# singularity: "docker://continuumio/miniconda3"

include: 'rules/common.smk'


rule all:
	input:
		# summaries
		'results/tables/sample_summary_reads.txt',
		'results/tables/sample_summary_features.txt',
		'results/tables/summary_total_features.txt',

		# feature tables
		# expand('results/tables/{space}/feature_table.biom',space=['aa','cdrs','cdr3']),
		# expand('results/tables/{space}/asvs.csv',space=['aa','cdrs','cdr3']),
		# expand('results/tables/aa/cdrs.csv'),
		# expand('results/tables/{space}/cdrs.csv',space=['aa','cdrs','cdr3']),


rule guids:
	input:
		samples='config/samples.tsv',
		runs='config/runs.tsv'
	output: 'config/guids.tsv'
	run:
		samples = read_delim_auto(input.samples)
		runs = read_delim_auto(input.runs)
		guids_to_samples(samples, runs).to_csv(output[0],index=False, sep="\t")

rule augment_metadata:
	input:
		samples='config/samples.tsv',
		metadata='config/metadata.csv'
	output: 'config/metadata_full.tsv'
	conda: 'envs/dada2.yaml'
	script: 'scripts/augment_metadata.R'

rule augment_metadata_csv:
	input:
		samples='config/samples.tsv',
		metadata='config/metadata.csv'
	output: 'config/metadata_full.csv'
	conda: 'envs/dada2.yaml'
	script: 'scripts/augment_metadata.R'



rule concat_sequences:
	input: lambda wc: str(Path(config['raw_sequence_dir']) / wc.run / f"Sample_{guid_to_sample_run(wc.sample)['sample']}/")
	output:
		fwd='intermediate/{run}/concat/fwd/{sample}.fastq.gz',
		rev='intermediate/{run}/concat/rev/{sample}.fastq.gz'
	shell:"""
	mkdir -p `dirname {output.fwd}`
	# cat {input}/*_R1_*.fastq.gz > {output.fwd}
	# cp {input}/*_R1_*.fastq.gz {output.fwd}
	ln -sr {input}/*_R1_*.fastq.gz {output.fwd}

	mkdir -p `dirname {output.rev}`
	# cat {input}/*_R2_*.fastq.gz > {output.rev}
	# cp {input}/*_R2_*.fastq.gz {output.rev}
	ln -sr {input}/*_R2_*.fastq.gz {output.rev}
	"""

rule trim_primers:
	group: 'trim_primers'
	# input: #unpack(trim_primers_input)
	# 	fwd='intermediate/partitioned/{sample}-R1-{partition}.fastq.gz',
	# 	rev='intermediate/partitioned/{sample}-R2-{partition}.fastq.gz'
	input: #unpack(trim_primers_input)
		fwd=ancient('intermediate/{run}/concat/fwd/{sample}.fastq.gz'),
		rev=ancient('intermediate/{run}/concat/rev/{sample}.fastq.gz')
	output:
		fwd='intermediate/{run}/primers_trimmed/fwd/{sample}.fastq.gz',
		rev='intermediate/{run}/primers_trimmed/rev/{sample}.fastq.gz',
		info='intermediate/{run}/primers_trimmed/info/{sample}.txt'
		# directory('intermediate/primers_trimmed/')
	params:
		primer_fwd=lambda wc: get_sample_library_param(wc.run, wc.sample, 'primer_fwd'),
		primer_rev=lambda wc: get_sample_library_param(wc.run, wc.sample, 'primer_rev')
	threads: 2
	conda: 'envs/cutadapt.yaml'
	shell: """
	PATH_FWD={input.fwd:q}
	PATH_REV={input.rev:q}
	echo $PATH_FWD
	echo $PATH_REV

	# reads may be in one of two orientations due to library prep. try both
	# and concat results
	cutadapt -g {params.primer_fwd} -G {params.primer_rev} \
		--cores={threads} \
		--discard-untrimmed \
		-o {output.fwd} -p {output.rev} \
		$PATH_FWD $PATH_REV > "{output.info}"
	"""

rule summarize_trim_primers:
	input: lambda wc: expand('intermediate/{run}/primers_trimmed/info/{sample}.txt', run=wc.run, sample=get_samples_per_run(wc.run))
	output:
		'intermediate/{run}/primers_trimmed/summary_reads.txt'
	conda:
		'envs/biopython-pysam.yaml'
	script: 'scripts/summarize_trim_primers.py'

# DENOISING
# -----------------------------------------------------------------------------

# filter and trim sequences from a sample using dada2; also plot the read quality
rule filter_trim_sequences:
	group: 'filter_trim_sequences'
	input: #'intermediate/primers_trimmed/'
		# unpack(input_filter_trim_sequences)
		fwd = 'intermediate/{run}/primers_trimmed/fwd/{sample}.fastq.gz',
		rev = 'intermediate/{run}/primers_trimmed/rev/{sample}.fastq.gz'
	output: #directory('intermediate/filtered_trimmed/')
		fwd     ='intermediate/{run}/filtered_trimmed/fwd/{sample}.fastq.gz',
		rev     ='intermediate/{run}/filtered_trimmed/rev/{sample}.fastq.gz',
		summary ='intermediate/{run}/filtered_trimmed/summary/{sample}.rds',
		plot_fwd='results/plots/dada2/quality/{run}/fwd/{sample}.png',
		plot_rev='results/plots/dada2/quality/{run}/rev/{sample}.png'
	log: 'results/logs/{run}/filter_trim_sequences/{sample}.log'
	threads: 2
	conda:
		'envs/dada2.yaml'
	params:
		maxEE     = 3,
		truncQ    = 2,
		minLen    = 100
	script:
		'scripts/filter_trim_sequences.R'

# learn the forward and reverse error model using dada2
def input_learn_error_model(wc):
	return dict(
		fwd=ancient(expand('intermediate/{run}/filtered_trimmed/fwd/{sample}.fastq.gz',run=wc.run, sample=get_samples_per_run(wc.run))),
		rev=ancient(expand('intermediate/{run}/filtered_trimmed/rev/{sample}.fastq.gz',run=wc.run, sample=get_samples_per_run(wc.run)))
	)

rule learn_error_model:
	input: unpack(input_learn_error_model)
	output:
		#'intermediate/error_model/error_model.rds' #directory('intermediate/error_model/')
		error_fwd = 'results/runs/{run}/error_model/error_model_fwd.rds',
		error_rev = 'results/runs/{run}/error_model/error_model_rev.rds',
		plot_fwd  = 'results/plots/dada2/{run}_error-model-fwd.png',
		plot_rev  = 'results/plots/dada2/{run}_error-model-rev.png'
	conda:
		'envs/dada2.yaml'
	threads: workflow.cores
	params:
		nbases=100000000
	log: "results/logs/{run}/learn_error_model.log"
	script: 'scripts/learn_error_model.R'
	# shell: "touch {output}"


rule dereplicate_sequences:
	input:
		fwd = 'intermediate/{run}/filtered_trimmed/fwd/{sample}.fastq.gz',
		rev = 'intermediate/{run}/filtered_trimmed/rev/{sample}.fastq.gz'
	output:
		fwd = 'intermediate/{run}/filtered_trimmed_dereplicated/fwd/{sample}.rds',
		rev = 'intermediate/{run}/filtered_trimmed_dereplicated/rev/{sample}.rds'
	conda:
		'envs/dada2.yaml'
	params:
		max_partition_size=75000
	log: "results/logs/{run}/dereplicate_sequences/{sample}.log"
	# script: 'scripts/dereplicate_sequences.R'
	script: 'scripts/dereplicate_sequences_partition.R'


# def input_denoise_sequences(wc):
# 	return dict(
# 		error_fwd = expand('intermediate/{run}/error_model/error_model_fwd.rds',run=wc.run),
# 		error_rev = expand('intermediate/{run}/error_model/error_model_rev.rds',run=wc.run),
# 		dereps_fwd = expand('intermediate/{run}/filtered_trimmed_dereplicated/fwd/{sample}.rds',run=wc.run, sample=get_samples_per_run(wc.run)),
# 		dereps_rev = expand('intermediate/{run}/filtered_trimmed_dereplicated/rev/{sample}.rds',run=wc.run, sample=get_samples_per_run(wc.run))
# 	)
#
# rule denoise_sequences:
# 	input: unpack(input_denoise_sequences)
# 	output:
# 		fwd = 'intermediate/{run}/denoise/fwd.rds',
# 		rev = 'intermediate/{run}/denoise/rev.rds'
# 	conda: 'envs/dada2.yaml'
# 	threads: workflow.cores,
# 	log: 'results/logs/{run}/denoise_sequences.log'
# 	script: 'scripts/denoise_sequences.R'

#
# rule partition_sequences:
# 	input:
# 	output: 'intermediate/{run}/partitioned/{direction}'

def input_denoise_sequences(wc):
	return dict(
		error = expand('results/runs/{run}/error_model/error_model_{direction}.rds',
			run=wc.run, direction=wc.direction),
		dereps = expand('intermediate/{run}/filtered_trimmed_dereplicated/{direction}/{sample}.rds',
			run=wc.run, sample=get_samples_per_run(wc.run), direction=wc.direction)
	)

rule denoise_sequences:
	input: unpack(input_denoise_sequences)
	output:
		'intermediate/{run}/denoise/{direction}.rds'
	conda: 'envs/dada2.yaml'
	# if more than 16 cores, use half; if less than 16, use all
	threads: (workflow.cores//2 if (workflow.cores >= 16) else workflow.cores)
	log: 'results/logs/{run}/denoise_sequences/{direction}.log'
	script: 'scripts/denoise_sequences_partition.R'

# def input_export_denoised_pairs(wc):
# 	return dict(
# 		dadas_fwd  = expand('intermediate/{run}/denoise/fwd.rds', run=wc.run),
# 		dadas_rev  = expand('intermediate/{run}/denoise/rev.rds', run=wc.run),
# 		dereps_fwd = expand('intermediate/{run}/filtered_trimmed_dereplicated/fwd/{sample}.rds',run=wc.run, sample=get_samples_per_run(wc.run)),
# 		dereps_rev = expand('intermediate/{run}/filtered_trimmed_dereplicated/rev/{sample}.rds',run=wc.run, sample=get_samples_per_run(wc.run))
# 	)
def input_export_denoised_pairs(wc):
	return dict(
		# temporary hack
		fwd = expand('intermediate/{run}/denoise/fwd.rds', run=wc.run),
		rev = expand('intermediate/{run}/denoise/rev.rds', run=wc.run),

		# dadas_fwd  = ancient(expand('intermediate/{run}/denoise/fwd/initial/{sample}.rds',              run=wc.run, sample=get_samples_per_run(wc.run))),
		# dadas_rev  = ancient(expand('intermediate/{run}/denoise/rev/initial/{sample}.rds',              run=wc.run, sample=get_samples_per_run(wc.run))),
		dereps_fwd = expand('intermediate/{run}/filtered_trimmed_dereplicated/fwd/{sample}.rds',run=wc.run, sample=get_samples_per_run(wc.run)),
		dereps_rev = expand('intermediate/{run}/filtered_trimmed_dereplicated/rev/{sample}.rds',run=wc.run, sample=get_samples_per_run(wc.run))
	)

rule export_denoised_pairs:
	input:
		unpack(input_export_denoised_pairs)
	output:
		mergers = 'intermediate/{run}/denoise/pairs.rds'
	log: 'results/logs/{run}/export_dadas_paired.log',
	threads: workflow.cores,
	conda: 'envs/dada2.yaml'
	script: 'scripts/export_denoised_pairs.R'

rule export_unique_pairs:
	input:
		mergers = 'intermediate/{run}/denoise/pairs.rds'
	output:
		unique_pairs = 'results/runs/{run}/denoise/unique_pairs.csv'
	log: 'results/logs/{run}/export_unique_pairs.log'
	threads: workflow.cores
	conda: 'envs/dada2.yaml'
	script: 'scripts/export_unique_pairs.R'

rule export_feature_table_folder:
	input:
		mergers = 'intermediate/{run}/denoise/pairs.rds'
	output:
		feature_table = directory('intermediate/{run}/denoise/feature_table')
	log: 'results/logs/{run}/export_feature_table.log'
	threads: workflow.cores
	conda: 'envs/dada2.yaml'
	script: 'scripts/export_feature_table_folder.R'

# todo: combine this with the R script above
rule export_feature_table_folder_biom:
	input: 'intermediate/{run}/denoise/feature_table'
	output: 'results/runs/{run}/denoise/feature_table.biom'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/export_feature_table_folder_biom.py'



# SUMMARIZE
# ---------

rule summarize_dereps:
	input: lambda wc: expand('intermediate/{{run}}/filtered_trimmed_dereplicated/fwd/{sample}.rds',sample = get_samples_per_run(wc.run))
	output: 'intermediate/{run}/filtered_trimmed_dereplicated/summary_features.txt'
	conda: 'envs/dada2.yaml'
	script: 'scripts/summarize_dereps.R'

def input_summarize_filter(wc):
	return {
		'filter_summaries': expand(
			'intermediate/{{run}}/filtered_trimmed/summary/{sample}.rds',
			sample = get_samples_per_run(wc.run)
		)
	}
rule summarize_filter:
	input: unpack(input_summarize_filter)
	output:
		summary_filter = 'intermediate/{run}/filtered_trimmed/summary_reads.txt',
	conda: 'envs/dada2.yaml'
	script: 'scripts/summarize_filter.R'
