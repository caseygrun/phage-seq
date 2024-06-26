
from pathlib import Path
import pandas as pd
from scripts.common import *

include: 'rules/common.smk'

rule downstream_all:
	input:
		# transformations
		expand('results/tables/{space}/transformed/{transform}/feature_table.biom', space=['aa','cdrs','cdr3'], transform=['log1p','scran']),

		# diversity
		expand('results/tables/{space}/alpha_diversity.txt',        space=['aa','cdrs','cdr3']),
		expand('intermediate/{space}/depth/sample_depth_stats.rds', space=['aa','cdrs','cdr3']),

		# phylogeny
		expand('intermediate/{space}/features/{partition}/{library}/asvs.nwk', space=['aa', 'cdr3'], partition=['top_asvs'], library=get_libraries()),

		# database
		expand('intermediate/{space}/features_db', space=['aa','cdr3']),

		# ordinations
		expand('intermediate/{space}/features/{partition}/{library}/transform/{transform}/ordination/{method}.ordination.gz',
			space=['cdr3', 'cdrs', 'aa'],partition=['all'],transform=['scran','log1p'],method=['TSVD','TSVD-TSNE'],library=get_libraries()),

		expand('results/tables/{space}/beta/all/pcoa/{metric}-{library}.ordination.gz',      space=['aa','cdrs','cdr3'], metric=['braycurtis','jaccard'], library=get_libraries()),
		
		# additional optional ordinations

		# expand('results/tables/{space}/beta/top_asvs/pcoa/{metric}-{library}.ordination.gz', space=['aa','cdrs','cdr3'], metric=['braycurtis','jaccard',
		# 	'weighted_unifrac',
		# 	'unweighted_unifrac'], library=get_libraries()),

		# expand('results/tables/{space}/beta/all/pcoa/{metric}-{library}.ordination.gz', space=['cdr3'], metric=[
		# 	'weighted_unifrac',
		# 	'unweighted_unifrac'], library=get_libraries()),

		# expand('results/tables/{space}/beta/top_asvs/pcoa/{metric}-{library}.ordination.gz', space=['cdr3'], metric=[
		# 	'deicode'], library=['alpaca']),

		# by experiment:
		# expand(
		# 	# perform two expansions: first get all valid combinations of `expt`
		# 	# and `phage_library`, so experiments with only only one library
		# 	# don't get represented twice;
		# 	#
		# 	# then get all combinations of {space}, {transform}, {method}
		# 	expand('intermediate/{{space}}/features/{expt}/{phage_library}/'
		# 			'transform/{{transform}}/ordination/{{method}}.ordination.gz',
		# 			zip,
		# 			# get all valid combinations of
		# 			**get_metadata_values(['expt','phage_library']).apply(lambda x: x.astype(str).str.lower())
		# 		),
		# 	space=['cdr3','aa'],
		# 	transform=['log1p'],
		# 	method=['TSVD','TSVD-TSNE'],
		# ),
		# enrichment model
		expand('results/tables/{space}/enrichment/null/ecdf.pickle', space=['cdr3'])

# TRANSFORMATIONS
# ---------------

rule transformations:
	input: 
		expand('results/tables/{space}/transformed/{transform}/feature_table.biom', space=['aa','cdrs','cdr3'], transform=['log1p','sqrt','scran']),

rule alpha_diversities:
	input:
		expand('results/tables/{space}/alpha_diversity.txt',        space=['aa','cdrs','cdr3']),
		expand('intermediate/{space}/depth/sample_depth_stats.rds', space=['aa','cdrs','cdr3'])

rule ordination_top:
	input:
		expand('results/tables/{space}/beta/top_asvs/pcoa/{metric}-{library}.ordination.gz', space=['aa','cdrs','cdr3'], metric=['braycurtis','jaccard',
			'weighted_unifrac',
			'unweighted_unifrac'], library=get_libraries()),

		expand('results/tables/{space}/beta/top_asvs/pcoa/{metric}-{library}.ordination.gz', space=['aa','cdrs','cdr3'], metric=[
			'deicode'], library=['alpaca']),

rule lme_antigens:
	input:
		expand('results/tables/cdr3/antigens/lme/top_enrichment_expt/{antigen}.csv', antigen=ANTIGENS.name)

rule norm_feature_table:
	wildcard_constraints:
		norm=list_to_regex(['scran'])
	input: 'results/tables/{space}/feature_table.biom'
	output: 'results/tables/{space}/transformed/{norm}/feature_table.biom'
	log: 'results/logs/{space}/transformed/{norm}.log'
	conda: 'envs/scran-deseq2.yaml',
	script: 'scripts/norm_ft.py'


rule transform_feature_table:
	wildcard_constraints:
		transform=list_to_regex(['log1p','sqrt'])
	input: 'results/tables/{space}/feature_table.biom'
	output: 'results/tables/{space}/transformed/{transform}/feature_table.biom'
	conda: 'envs/biopython-pysam.yaml',
	script: 'scripts/transform_ft.py'

rule norm_feature_table_library:
	wildcard_constraints:
		norm=list_to_regex(['scran'])
	input: 'intermediate/{space}/features/{partition}/{library}/feature_table.biom'
	output: 'intermediate/{space}/features/{partition}/{library}/transform/{norm}/feature_table.biom'
	log: 'results/logs/{space}/features/{partition}/{library}/transform/{norm}.log'
	conda: 'envs/scran-deseq2.yaml',
	script: 'scripts/norm_ft.py'


rule transform_feature_table_library:
	wildcard_constraints:
		transform=list_to_regex(['log1p','sqrt'])
	input: 'intermediate/{space}/features/{partition}/{library}/feature_table.biom'
	output: 'intermediate/{space}/features/{partition}/{library}/transform/{transform}/feature_table.biom'
	conda: 'envs/biopython-pysam.yaml',
	script: 'scripts/transform_ft.py'


# rule transform_feature_table_library:
# 	input: 'intermediate/{space}/features/{partition}/{library}/feature_table.biom'
# 	output: 'intermediate/{space}/features/{partition}/{library}/transformed/{transform}/feature_table.biom'
# 	conda: 'envs/biopython-pysam.yaml',
# 	script: 'scripts/transform_ft.py'


# EXTRACT CDRS AND CREATE DATABASES FOR INTERACTIVITY
# ---------------------------------------------------

rule extract_cdrs_:
	input: 'results/tables/aa/asvs.csv'
	output: 'results/tables/aa/cdrs.csv'
	conda: 'envs/biopython-pysam.yaml',
	params:
		library_CDRs = {library: get_library_param(library, 'CDRs') for library in get_libraries() }
	script: 'scripts/extract_cdrs.py'




rule cdrs_to_sqlite:
	input:
		cdrs='results/tables/aa/cdrs.csv',
		aa='results/tables/aa/asvs.csv'
	output: 'intermediate/aa/asvs.db'
	conda: 'envs/biopython-pysam.yaml'
	# quote the heredoc delimiter, e.g. <<"EOF" instead of <<EOF to prevent the
	# shell from trying to quote-expand the backticked SQL statements
	shell: """sqlite3 {output} <<"EOF"
.mode csv
.import {input.cdrs} cdrs
.import {input.aa} aa
CREATE INDEX cdrs_CDR3ID ON `cdrs` (`CDR3ID` COLLATE NOCASE);
CREATE UNIQUE INDEX cdrs_aaSVID ON `cdrs` (`aaSVID` COLLATE NOCASE);
CREATE UNIQUE INDEX aa_aaSVID ON `aa` (`aaSVID` COLLATE NOCASE);
CREATE INDEX aa_aaSVID_aligned ON `aa`(`aaSVID`, `aligned`);
CREATE INDEX aa_library_aaSVID_aligned ON `aa`(`library`, `aaSVID`, `aligned`);
EOF
	"""

rule link_asvs_for_mmseqs2:
	wildcard_constraints:
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input: 'results/tables/{space}/asvs.csv'
	output: 'intermediate/{space}/asvs.csv'
	shell: """
	cp "{input}" "{output}"
	"""

rule asvs_to_mmseqs2db:
	input: 'intermediate/{space}/asvs.fasta'
	output: directory('intermediate/{space}/features_db')
	conda: 'envs/mmseqs2-vsearch.yaml'
	shell: """
	rm -r -f {output}
	mkdir -p {output}
	mmseqs createdb {input} {output}/features
	mmseqs createindex {output}/features {config[scratch]}
	"""

# moved to feature_table.smk
# --------------------------
# rule collapse_feature_table_clonotype:
# 	input:
# 		feature_table = 'intermediate/{run}/aa_filter/feature_table.biom',
# 		mapping =       'intermediate/{run}/aa_filter/aa_cdrs.csv'
# 	output:
# 		feature_table = 'intermediate/{run}/cdrs/feature_table.biom'
# 	params:
# 		mapping_cols=['aaSVID','clonotypeID']
# 	log: 'results/logs/{run}/collapse_feature_table_clonotype.log'
# 	conda: 'envs/biopython-pysam.yaml'
# 	script: 'scripts/collapse_feature_table_aa.py'


# use rule collapse_feature_table_clonotype as collapse_feature_table_clonotype_all with:
# 	input:
# 		feature_table = 'results/tables/aa/feature_table.biom',
# 		mapping =       'results/tables/aa/cdrs.csv'
# 	output:
# 		feature_table = 'results/tables/cdrs/feature_table.biom'
# 	params:
# 		mapping_cols=['aaSVID','clonotypeID']
# 	log: 'results/logs/collapse_feature_table_clonotype.log'

# 	threads: workflow.cores//2

# use rule collapse_feature_table_clonotype_all as collapse_feature_table_cdr3_all with:
# 	output:
# 		feature_table = 'results/tables/cdr3/feature_table.biom'
# 	params:
# 		mapping_cols=['aaSVID','CDR3ID']
# 	log: 'results/logs/collapse_feature_table_cdr3.log'


# moved to common.smk
# --------------------------
# rule copy_cdrs_clonotype_:
# 	wildcard_constraints:
# 		space=list_to_regex(['cdrs','cdr3'])
# 	input: 'results/tables/aa/cdrs.csv'
# 	output: 'results/tables/{space}/asvs.csv'
# 	threads: workflow.cores//2
# 	conda: 'envs/biopython-pysam.yaml'
# 	script: 'scripts/copy_cdrs_clonotype.py'


# TODO: this should really operate separately for each phage_library; will require refactoring
rule simulate_enrichment:
	wildcard_constraints:
		space=list_to_regex(['aa', 'cdrs','cdr3'])
	input: 
		feature_table = 'results/tables/{space}/feature_table.biom',
		metadata = 'config/metadata_full.csv'
	output: "results/tables/{space}/enrichment/null/ecdf.pickle"
	params:
		query=config['input_query']
	conda: 'envs/biopython-pysam.yaml',
	script: 'scripts/simulate_enrichment.py'

# SUMMARIZE
# ---------

rule summarize_samples_reads:
	wildcard_constraints:
		run=list_to_regex(get_runs())
	input:
		primers_trimmed = 'intermediate/{run}/primers_trimmed/summary_reads.txt',
		filtered = 'intermediate/{run}/filtered_trimmed/summary_reads.txt',
		denoised = 'intermediate/{run}/denoise/summary_reads.txt',
		aligned = 'intermediate/{run}/align/summary_reads.txt',
		nochimeras = 'intermediate/{run}/nochimeras/summary_reads.txt',
		aa = 'intermediate/{run}/aa/summary_reads.txt',
		aa_clustered = 'intermediate/{run}/aa_cluster/summary_reads.txt',
		aa_filtered = 'intermediate/{run}/aa_filter/summary_reads.txt'
	output:
		summary = 'intermediate/{run}/sample_summary_reads.txt'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/summarize_samples.py'

rule summarize_samples_features:
	wildcard_constraints:
		run=list_to_regex(get_runs())
	input:
		derep = 'intermediate/{run}/filtered_trimmed_dereplicated/summary_features.txt',
		denoised = 'intermediate/{run}/denoise/summary_features.txt',
		aligned = 'intermediate/{run}/align/summary_features.txt',
		nochimeras = 'intermediate/{run}/nochimeras/summary_features.txt',
		aa = 'intermediate/{run}/aa/summary_features.txt',
		aa_clustered = 'intermediate/{run}/aa_cluster/summary_features.txt',
		aa_filtered = 'intermediate/{run}/aa_filter/summary_features.txt',
		cdrs = 'intermediate/{run}/cdrs/summary_features.txt',
		cdr3 = 'intermediate/{run}/cdr3/summary_features.txt'
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
		denoise      = expand('intermediate/{run}/denoise/feature_table.biom'    , run = get_runs()) ,
		align        = expand('intermediate/{run}/align/feature_table.biom'      , run = get_runs()) ,
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



# ALPHA DIVERSITY
# ---------------

rule alpha_diversity:
	input:
		feature_table='results/tables/{space}/feature_table.biom'
		# feature_table='intermediate/features/top_asvs/feature_table.biom
	output:
		alpha_diversity='results/tables/{space}/alpha_diversity.txt'
	conda: 'envs/breakaway.yaml',
	threads: workflow.cores
	log: 'results/logs/{space}/alpha_diversity.log'
	script: 'scripts/diversity.R'

# rule alpha_diversity_run:
# 	input:
# 		feature_table='intermediate/aa/feature_table.biom'
# 	output:
# 		alpha_diversity='intermediate/aa/alpha_diversity.txt'
# 	conda: 'envs/breakaway.yaml',
# 	threads: workflow.cores
# 	log: 'results/logs/aa/alpha_diversity_run.log'
# 	script: 'scripts/diversity.R'

rule depth_stats:
	input: 'results/tables/{space}/feature_table.biom'
	output: 'intermediate/{space}/depth/sample_depth_stats.rds'
	log:    'results/logs/{space}/depth_stats.log'
	conda: 'envs/breakaway.yaml'
	threads: workflow.cores
	script: 'scripts/depth_stats.R'


# BETA DIVERSITY
# --------------

ruleorder: filter_features_abundance_prevalence > features_all > filter_feature_table_by_expt > link_asvs_for_mmseqs2

rule filter_features_abundance_prevalence:
	input:
		feature_table='results/tables/{space}/feature_table.biom',
		feature_data='results/tables/{space}/asvs.csv'
	output:
		feature_table='intermediate/{space}/features/top_asvs/feature_table.biom',
		feature_data= 'intermediate/{space}/features/top_asvs/asvs.csv'
	params:
		min_abundance=10,
		min_prevalence=2
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/filter_features_abundance_prevalence.py'

rule features_all:
	input:
		feature_table='results/tables/{space}/feature_table.biom',
		feature_data= 'results/tables/{space}/asvs.csv'
	output:
		feature_table='intermediate/{space}/features/all/feature_table.biom',
		feature_data= 'intermediate/{space}/features/all/asvs.csv'
	shell:"""
	ln -sr {input.feature_table} {output.feature_table}
	ln -sr {input.feature_data} {output.feature_data}
	"""

rule filter_feature_table_by_expt:
	wildcard_constraints:
		expt=r'[^/]+'
	input:
		feature_table = 'intermediate/{space}/features/all/feature_table.biom',
		feature_data  = 'intermediate/{space}/features/all/asvs.csv',
		metadata      = 'config/samples.tsv'
	output:
		feature_table = 'intermediate/{space}/features/{expt}/feature_table.biom',
		feature_data =  'intermediate/{space}/features/{expt}/asvs.csv'
	params:
		param='expt'
	threads: workflow.cores
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/filter_feature_table_by_metadata.py'


rule filter_feature_table_by_library:
	wildcard_constraints:
		partition=r'[^/]+',
		library=list_to_regex(get_libraries())
	input:
		feature_table = 'intermediate/{space}/features/{partition}/feature_table.biom',
		feature_data =  'intermediate/{space}/features/{partition}/asvs.csv'
	output:
		feature_table = 'intermediate/{space}/features/{partition}/{library}/feature_table.biom',
		feature_data =  'intermediate/{space}/features/{partition}/{library}/asvs.csv'
	log: 'results/logs/{space}/features/{partition}/{library}/filter_feature_table_by_library.log'
	threads: workflow.cores
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/filter_feature_table_by_library.py'


# duplicate
# ruleorder: fasttree_asvs > fasttree_ref_align_asvs

# rule fasttree_asvs:
# 	wildcard_constraints:
# 		partition=list_to_regex(['top_asvs','all'])
# 	input:  'intermediate/{space}/features/{partition}/{library}/asvs.fasta'
# 	output: 'intermediate/{space}/features/{partition}/{library}/asvs.nwk'
# 	log:    'results/logs/{space}/features/{partition}/{library}/fasttree_asvs.log'
# 	conda: 'envs/msa.yaml'
# 	threads: 4
# 	shell: """
# 	{{
# 	OMP_NUM_THREADS={threads}
# 	FastTreeMP < {input} > {output}
# 	}} >{log} 2>&1
# 	"""

rule msa_asvs_mafft:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input:  'intermediate/{space}/features/{partition}/{library}/asvs.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/msa.fasta'
	log:    'results/logs/{space}/features/{partition}/{library}/msa_asvs.log'
	conda: 'envs/msa.yaml'
	threads: 4
	shell:"""
	mafft --ep 0.123 --lep -0.5 --leavegappyregion --auto --thread {threads} {input} > {output} 2> {log}
	"""

rule msa_asvs_clustalo:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input:  'intermediate/{space}/features/{partition}/{library}/asvs.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/msa-clustalo.fasta'
	log:    'results/logs/{space}/features/{partition}/{library}/msa-clustalo.log'
	conda: 'envs/msa.yaml'
	threads: 4
	shell:"""
	clustalo --threads {threads} --in {input} --out {output} --log {log}
	"""

rule msa_asvs_muscle:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input:  'intermediate/{space}/features/{partition}/{library}/asvs.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/msa-muscle.fasta'
	log:    'results/logs/{space}/features/{partition}/{library}/msa-muscle.log'
	conda: 'envs/msa.yaml'
	threads: 4
	shell:"""
	muscle -super5 {input} -output {output} -threads {threads} 2>&1 > {log}
	"""



rule trim_alignment:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input:  'intermediate/{space}/features/{partition}/{library}/msa.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/msa-trimmed.fasta'
	log:    'results/logs/{space}/features/{partition}/{library}/trim_alignment.log'
	conda: 'envs/msa.yaml'
	shell:"""
	trimal -in {input} -out {output} -automated1 2> {log}
	if ! [ -s "{output}" ] ; then
		cp {input} {output}
		echo "Alignment trimming on feature MSA {input} resulted in an empty alignment (all positions empty); will use un-trimmed alignment instead to construct ML phylogeny for this feature space" >> {log}
	fi
	"""

rule fasttree_msa_asvs:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input:  'intermediate/{space}/features/{partition}/{library}/msa-trimmed.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/msa.nwk'
	log:    'results/logs/{space}/features/{partition}/{library}/fasttree_msa_asvs.log'
	conda: 'envs/msa.yaml'
	threads: 4
	shell: """
	{{
	OMP_NUM_THREADS={threads}
	FastTreeMP < {input} > {output}
	}} >{log} 2>&1
	"""


rule fasttree_msa_asvs_method:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3']),
		method=list_to_regex(['clustalo','muscle'])
		# method=r'^(?!trimmed$)'
	input:  'intermediate/{space}/features/{partition}/{library}/msa-{method}.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/msa-{method}.nwk'
	log:    'results/logs/{space}/features/{partition}/{library}/fasttree_msa_asvs_{method}.log'
	conda: 'envs/msa.yaml'
	threads: 4
	shell: """
	{{
	OMP_NUM_THREADS={threads}
	FastTreeMP < {input} > {output}
	}} >{log} 2>&1
	"""

rule fasttree_ref_align_asvs:
	wildcard_constraints:
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['aa','cdrs','cdr3'])
	input:  'intermediate/{space}/features/{partition}/{library}/asvs.fasta'
	output: 'intermediate/{space}/features/{partition}/{library}/asvs.nwk'
	log:    'results/logs/{space}/features/{partition}/{library}/fasttree_ref_align_asvs.log'
	conda: 'envs/msa.yaml'
	threads: 4
	shell: """
	{{
	OMP_NUM_THREADS={threads}
	FastTreeMP < {input} > {output}
	}} >{log} 2>&1
	"""



rule beta_diversity_distance:
	wildcard_constraints:
		metric=list_to_regex(['braycurtis','jaccard','shared_reads'])
	input:
		feature_table='intermediate/{space}/features/{partition}/feature_table.biom',
		metadata='config/samples.tsv'
	output: 'results/tables/{space}/beta/{partition}/distance/{metric}.tsv'
	params:
		min_asv_abundance = 10
	log: 'results/logs/{space}/beta/{partition}/distance/{metric}.log'
	threads: workflow.cores,
	conda: 'envs/beta-diversity.yaml'
	script: 'scripts/beta_diversity.py'

rule beta_diversity_distance_library:
	wildcard_constraints:
		metric=list_to_regex(['braycurtis','jaccard','shared_reads'])
	input:
		feature_table = 'intermediate/{space}/features/{partition}/{library}/feature_table.biom',
		metadata='config/samples.tsv'
	output:  'results/tables/{space}/beta/{partition}/distance/{metric}-{library}.tsv'
	log: 'results/logs/{space}/beta/{partition}/distance/{metric}-{library}.log'
	conda: 'envs/beta-diversity.yaml'
	threads: workflow.cores
	script: 'scripts/beta_diversity.py'

rule beta_diversity_distance_phylogetic:
	wildcard_constraints:
		metric=list_to_regex(['weighted_unifrac','unweighted_unifrac']),
		partition=list_to_regex(['top_asvs','all']),
		space=list_to_regex(['na','aa','cdrs','cdr3'])
	input:
		feature_table='intermediate/{space}/features/{partition}/{library}/feature_table.biom',
		phylogeny=    'intermediate/{space}/features/{partition}/{library}/msa.nwk'
	output: 'results/tables/{space}/beta/{partition}/distance/{metric}-{library}.tsv'
	log:    'results/logs/{space}/beta/{partition}/distance/{metric}-{library}.log'
	conda: 'envs/beta-diversity.yaml'
	threads: workflow.cores
	script: 'scripts/beta_diversity_distance_phylogetic.py'

ruleorder: beta_diversity_rpcoa > beta_diversity_pcoa

rule beta_diversity_pcoa:
	input:
		distance_matrix='results/tables/{space}/beta/{partition}/distance/{metric}.tsv'
	output:
		pcoa='results/tables/{space}/beta/{partition}/pcoa/{metric}.ordination.gz'
		# ,biplot='results/tables/beta/pcoa/{metric}-biplot.ordination'
	conda: 'envs/beta-diversity.yaml'
	log:    'results/logs/{space}/beta/{partition}/pcoa/{metric}.log'
	script: 'scripts/beta_diversity_pcoa.py'

rule beta_diversity_rpcoa:
	wildcard_constraints:
		metric=list_to_regex(['deicode'])
	input:
		feature_table='intermediate/{space}/features/top_asvs/{library}/feature_table.biom'
	output:
		distance_matrix='results/tables/{space}/beta/top_asvs/distance/{metric}-{library}.tsv',
		pcoa=           'results/tables/{space}/beta/top_asvs/pcoa/{metric}-{library}.ordination.gz'
		# ,biplot='results/tables/beta/pcoa/{metric}-biplot.ordination'
	threads: workflow.cores
	log:    'results/logs/{space}/beta/top_asvs/pcoa/{metric}-{library}.log'
	conda: 'envs/beta-diversity.yaml'
	script: 'scripts/beta_diversity_rpcoa.py'


rule ordination:
	input: 'intermediate/{space}/features/{partition}/{library}/transform/{transform}/feature_table.biom'
	output:
		skbio='intermediate/{space}/features/{partition}/{library}/transform/{transform}/ordination/{method}.ordination.gz',
		sklearn='intermediate/{space}/features/{partition}/{library}/transform/{transform}/ordination/{method}.pickle'
	conda: 'envs/beta-diversity.yaml'
	# threads: workflow.cores//2
	log:    'results/logs/{space}/features/{partition}/{library}/transform/{transform}/ordination/{method}.log'
	threads: workflow.cores,
	script: 'scripts/ordinate.py'

# CLASSIFIER TUNING
# -----------------

rule lme_enrichment:
	input: 'intermediate/{space}/features/{partition}/feature_table.csv'
	output: 'results/tables/{space}/antigens/lme/{partition}/{antigen}.csv'
	conda: 'envs/lme4.yaml'
	threads: 2
	script: 'scripts/lme_enrichment.py'

rule design_ag_matrix:
	input:
		'results/tables/aa/feature_table.biom',
		'results/tables/cdr3/feature_table.biom',
		'config/metadata-phenotypes.csv',
		'config/phenotypes.csv',
	output:
		designs = 'intermediate/learning/designs.pickle',
		ag_matrix = 'intermediate/learning/ag_matrix.pickle',
	log: 'results/logs/learning/design_ag_matrix.log'
	threads: workflow.cores,
	conda: 'envs/xgb.yaml'
	script: 'scripts/design_ag_matrix.py'

rule xgb_tuning_manual:
	input:
		designs='intermediate/learning/designs.pickle',
		ag_matrix='intermediate/learning/ag_matrix.pickle',
	output: 'intermediate/learning/xgb/manual/{design}/{antigen}.pickle'
	log: 'results/logs/learning/xgb/manual/{design}/{antigen}.log'
	conda: 'envs/xgb.yaml'
	threads: workflow.cores // 2
	script: 'scripts/xgb_tuning_manual.py'

rule xgb_tuning_bayes:
	input:
		designs='intermediate/learning/designs.pickle',
		ag_matrix='intermediate/learning/ag_matrix.pickle',
	output: 'intermediate/learning/xgb/bayes/{design}/{antigen}.pickle'
	log: 'results/logs/learning/xgb/bayes/{design}/{antigen}.log'
	conda: 'envs/xgb.yaml'
	threads: workflow.cores // 2
	script: 'scripts/xgb_tuning_bayes.py'

rule xgb_tuning_random:
	input:
		designs='intermediate/learning/designs.pickle',
		ag_matrix='intermediate/learning/ag_matrix.pickle',
	output: 'intermediate/learning/xgb/random/{design}/{antigen}.pickle'
	log: 'results/logs/learning/xgb/random/{design}/{antigen}.log'
	conda: 'envs/xgb.yaml'
	threads: workflow.cores // 2
	script: 'scripts/xgb_tuning_random.py'


rule logistic_regression_tuning_bayes:
	input:
		designs='intermediate/learning/designs.pickle',
		ag_matrix='intermediate/learning/ag_matrix.pickle',
	output: 'intermediate/learning/logistic_regression/bayes/{design}/{antigen}.pickle'
	log: 'results/logs/learning/logistic_regression/bayes/{design}/{antigen}.log'
	conda: 'envs/xgb.yaml'
	threads: workflow.cores // 2
	script: 'scripts/logistic_regression_tuning_bayes.py'


test_antigens = config['test_antigens'] if 'test_antigens' in config else []

rule xgb_tunings_manual_cdrs:
	input: expand('intermediate/learning/xgb/manual/{design}/{antigen}.pickle', design=['CDR3:enr', 'CDR3:R5+enr', 'CDR3:R2345+enr'], antigen=test_antigens)

rule xgb_tunings_manual_aa:
	input: expand('intermediate/learning/xgb/manual/{design}/{antigen}.pickle', design=['aa:enr', 'aa:R5+enr'], antigen=test_antigens)

rule xgb_tunings_random:
	input: expand('intermediate/learning/xgb/random/{design}/{antigen}.pickle', design=['CDR3:enr', 'CDR3:R5+enr', 'CDR3:R2345+enr', 'aa:enr', 'aa:R5+enr'], antigen=test_antigens)

rule xgb_tunings_bayes:
	input: expand('intermediate/learning/xgb/bayes/{design}/{antigen}.pickle', design=['CDR3:enr', 'CDR3:R5+enr', 'CDR3:R2345+enr', 'aa:enr', 'aa:R5+enr'], antigen=test_antigens)

rule logistic_regression_tunings_bayes:
	input: expand('intermediate/learning/logistic_regression/bayes/{design}/{antigen}.pickle', design=['CDR3:enr', 'CDR3:R5+enr', 'CDR3:R2345+enr', 'aa:enr', 'aa:R5+enr'], antigen=test_antigens)
