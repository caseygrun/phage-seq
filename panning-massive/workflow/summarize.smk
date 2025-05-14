# include: './rules/common.smk'

include: "preprocess.smk"
include: "feature_table.smk"
include: "downstream.smk"

rule summarize_all:
	input:
		# summaries
		expand('intermediate/{run}/sample_summary_reads.txt',    run=get_runs()),
		expand('intermediate/{run}/sample_summary_features.txt', run=get_runs()),
		expand('intermediate/{run}/summary_total_features.txt',  run=get_runs()),
		'intermediate/summary_total_features.txt',
		'results/tables/sample_summary_reads.txt',
		'results/tables/sample_summary_features.txt',
		'results/tables/summary_total_features.txt'

rule summarize_samples_reads:
	wildcard_constraints:
		run=list_to_regex(get_runs())
	input:
		primers_trimmed = 'intermediate/{run}/primers_trimmed/summary_reads.txt',
		filtered        = 'intermediate/{run}/filtered_trimmed/summary_reads.txt',
		denoised        = 'intermediate/{run}/denoise/summary_reads.txt',
		aligned         = 'intermediate/{run}/align/summary_reads.txt',
		nochimeras      = 'intermediate/{run}/nochimeras/summary_reads.txt',
		aa              = 'intermediate/{run}/aa/summary_reads.txt',
		aa_clustered    = 'intermediate/{run}/aa_cluster/summary_reads.txt',
		aa_filtered     = 'intermediate/{run}/aa_filter/summary_reads.txt'
	output:
		summary = 'intermediate/{run}/sample_summary_reads.txt'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/summarize_samples.py'

rule summarize_samples_features:
	wildcard_constraints:
		run=list_to_regex(get_runs())
	input:
		derep        = 'intermediate/{run}/filtered_trimmed_dereplicated/summary_features.txt',
		denoised     = 'intermediate/{run}/denoise/summary_features.txt',
		aligned      = 'intermediate/{run}/align/summary_features.txt',
		nochimeras   = 'intermediate/{run}/nochimeras/summary_features.txt',
		aa           = 'intermediate/{run}/aa/summary_features.txt',
		aa_clustered = 'intermediate/{run}/aa_cluster/summary_features.txt',
		aa_filtered  = 'intermediate/{run}/aa_filter/summary_features.txt',
		cdrs         = 'intermediate/{run}/cdrs/summary_features.txt',
		cdr3         = 'intermediate/{run}/cdr3/summary_features.txt'
	output:
		summary = 'intermediate/{run}/sample_summary_features.txt'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/summarize_samples.py'

rule summarize_total_features_run:
	wildcard_constraints:
		run=list_to_regex(get_runs())
	input:
		# derep = 'intermediate/{run}/filtered_trimmed_dereplicated/summary_total_features.txt',
		denoise     = 'intermediate/{run}/denoise/summary_total_features.txt',
		align       = 'intermediate/{run}/align/summary_total_features.txt',
		nochimeras  = 'intermediate/{run}/nochimeras/summary_total_features.txt',
		aa          = 'intermediate/{run}/aa/summary_total_features.txt',
		aa_cluster  = 'intermediate/{run}/aa_cluster/summary_total_features.txt',
		aa_filter   = 'intermediate/{run}/aa_filter/summary_total_features.txt',
		cdrs        = 'intermediate/{run}/cdrs/summary_total_features.txt',
		cdr3        = 'intermediate/{run}/cdr3/summary_total_features.txt'
	output:
		summary = 'intermediate/{run}/summary_total_features.txt'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/summarize_total_features_run.py'

rule summarize_total_features:
	input:
		# derep = 'intermediate/{run}/filtered_trimmed_dereplicated/summary_total_features.txt',
		denoise      = expand('results/runs/{run}/denoise/feature_table.biom'    , run = get_runs()) ,
		align        = expand('results/runs/{run}/align/feature_table.biom'      , run = get_runs()) ,
		nochimeras   = expand('intermediate/{run}/nochimeras/feature_table.biom' , run = get_runs()) ,
		aa           = expand('intermediate/{run}/aa/feature_table.biom'         , run = get_runs()) ,
		aa_cluster   = expand('intermediate/{run}/aa_cluster/feature_table.biom' , run = get_runs()) ,
		aa_filter    = expand('intermediate/{run}/aa_filter/feature_table.biom'  , run = get_runs()) ,
		cdrs         = expand('intermediate/{run}/cdrs/feature_table.biom'       , run = get_runs()) ,
		cdr3         = expand('intermediate/{run}/cdr3/feature_table.biom'       , run = get_runs()) ,
	output:
		summary = 'intermediate/summary_total_features.txt'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/summarize_total_features.py'

rule concat_sample_summaries:
	wildcard_constraints:
		kind=list_to_regex(['reads','features'])
	input: expand('intermediate/{run}/sample_summary_{{kind}}.txt', run=get_runs())
	output: 'results/tables/sample_summary_{kind}.txt'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/concat_sample_summaries.py'

rule concat_summaries_total_features:
	input: (expand('intermediate/{run}/summary_total_features.txt', run=get_runs()) + ['intermediate/summary_total_features.txt'])
	output: 'results/tables/summary_total_features.txt'
	params:
		row='run_id'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/concat_sample_summaries.py'
