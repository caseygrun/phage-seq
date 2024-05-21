include: 'rules/common.smk'

rule feature_table_all:
	input:
		# feature tables
		expand('results/tables/{space}/feature_table.biom',space=['aa','cdrs','cdr3']),
		expand('results/tables/{space}/asvs.csv',space=['aa','cdrs','cdr3']),
		expand('results/tables/aa/cdrs.csv'),

		expand('intermediate/{run}/{space}/summary_features.txt'      , run=get_runs(), space=['aa','cdrs','cdr3']),
		expand('intermediate/{run}/{space}/summary_total_features.txt', run=get_runs(), space=['aa','cdrs','cdr3']),
		expand('intermediate/{run}/{phase}/summary_reads.txt'         , run=get_runs(), phase=['align','aa', 'aa_cluster', 'aa_filter'])


# ASV FILTERING AND FEATURE TABLES
# -----------------------------------------------------------------------------

# TODO: implement
rule filter_feature_table_nochimeras:
	input:
		feature_table = 'results/runs/{run}/denoise/feature_table.biom',
		unique_pairs = 'results/runs/{run}/denoise/unique_pairs.csv'
	output:
		feature_table = 'intermediate/{run}/nochimeras/feature_table.biom',
		unique_pairs = 'intermediate/{run}/nochimeras/unique_pairs.csv',
		chimeras = 'intermediate/{run}/nochimeras/chimeras.csv'
	log: 'results/logs/{run}/filter_feature_table_nochimeras.log'
	threads: workflow.cores
	conda: 'envs/mmseqs2-vsearch.yaml'
	# script: 'scripts/filter_feature_table_nochimeras.py'
	shell: """
	cp {input.feature_table} {output.feature_table}
	cp {input.unique_pairs} {output.unique_pairs}
	echo "ID,seq" > {output.chimeras}
	"""


def build_alignment_index_input(wildcards):
	assert wildcards['library'] in config['libraries'], ("Reference sequence for "
		f"{wildcards['library']} not found in configuration file. "
		"Cannot build reference database to align nucleic acid ASVs.")
	reference = config['libraries'][wildcards['library']]['reference']
	return reference

rule build_alignment_index:
	wildcard_constraints:
		library=list_to_regex(get_libraries())#"(" + "|".join(get_libraries()) + ")"
	input: build_alignment_index_input
	output:
		dir = directory("resources/references/{library}/"),
		fasta = 'resources/references/{library}/{library}.fasta'
	log: 'results/logs/build_alignment_index/{library}.log'
	conda: 'envs/bowtie2.yaml'
	script: 'scripts/build_alignment_index.py'

rule unique_pairs_csv_to_fasta:
	input:
		csv = '{base}.csv'
	output:
		fwd = '{base}_fwd.fasta',
		rev = '{base}_rev.fasta'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/unique_pairs_csv_to_fasta.py'

rule align_vhh_asv_pairs_na:
	input: #unpack(align_vhh_asvs_nucleic_acid_input)
		fwd = 'intermediate/{run}/nochimeras/unique_pairs_fwd.fasta',
		rev = 'intermediate/{run}/nochimeras/unique_pairs_rev.fasta',
		index = "resources/references/{library}/"
	output:
		coord_sorted = "results/runs/{run}/align/ASVs/{library}.bam",
		name_sorted = "results/runs/{run}/align/ASVs/{library}-name-sorted.bam"
	log: 'results/logs/{run}/align_vhh_asv_pairs_na/{library}.log'
	params:
		n_penalty = 0,
		n_ceil = 64
	threads: workflow.cores,
	conda: 'envs/bowtie2.yaml'
	script: 'scripts/align_vhh_asv_pairs_na.py'

# use BAM alignment of ASV sequences to produce an implied full-length
# reference sequence. ASVs are renamed according to read pairs
rule merge_aligned_vhh_asv_pairs_na:
	input:
		alignment = 'results/runs/{run}/align/ASVs/{library}-name-sorted.bam',
		reference = 'resources/references/{library}/{library}.fasta'
	output: 'results/runs/{run}/merge/{library}_merged_reads.csv.gz',
	log: 'results/logs/{run}/merge_aligned_vhh_asv_pairs_na/{library}.log'
	params:
		library=lambda wc: get_library_params(wc.library)
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/merge_aligned_vhh_asv_pairs_na.py'

# perform zero-width greedy OTU clustering
rule cluster_aligned_vhh_asv_pairs_na:
	input: 'results/runs/{run}/merge/{library}_merged_reads.csv.gz'
	output:
		clusters = 'results/runs/{run}/merge/{library}_clusters_na.csv',
		na = 'results/runs/{run}/merge/{library}_na.csv',
	threads: workflow.cores,
	log: 'results/logs/{run}/cluster_aligned_vhh_asv_pairs_na/{library}.log'
	conda: 'envs/mmseqs2-vsearch.yaml'
	script: 'scripts/cluster_aligned_vhh_asv_pairs_na.py'

# collapse feature table columns according to complete nucleic acid ASVs. reads
# which did not map to the reference will be removed.
rule collapse_feature_table_aligned_na:
	input:
		feature_table = 'intermediate/{run}/nochimeras/feature_table.biom',
		aligned_asvs = expand('results/runs/{{run}}/merge/{library}_na.csv', library=get_libraries()) #'results/tables/vhh_asvs.csv'
	output:
		feature_table = 'results/runs/{run}/align/feature_table.biom'
		# , summary = 'intermediate/{run}/align/summary_reads.txt'
	log: 'results/logs/{run}/filter_feature_table_aligned.log'
	conda: 'envs/biopython-pysam.yaml'
	threads: workflow.cores
	script: 'scripts/collapse_feature_table_aligned_na.py'

# translate ASVs to amino acid ASVs
rule translate_merged_aligned_vhh_asv_pairs:
	input:
		na = 'results/runs/{run}/merge/{library}_na.csv',
		reference = 'resources/references/{library}/{library}.fasta'
	output:
		# ASVID (hash of nucleic acid sequence), aaSVID (hash of ungapped amino acid sequence)
		na_to_aa = 'results/runs/{run}/aa/{library}_na_to_aa.csv',

		# aaSVID, translated = amino acid sequence
		aa = 'intermediate/{run}/aa/{library}_pairs_aa.csv'
	params:
		library=lambda wc: get_library_params(wc.library)
	threads: 4
	log: 'results/logs/{run}/translate_merged_aligned_vhh_asv_pairs/{library}.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/translate_merged_aligned_vhh_asv_pairs.py'

# converts the nucleic acid ASV feature table to an amino acid ASV (aaSV)
# feature table
rule collapse_feature_table_aa_unfiltered:
	input:
		feature_table = 'results/runs/{run}/align/feature_table.biom',
		mapping = expand('results/runs/{{run}}/aa/{library}_na_to_aa.csv',library=get_libraries())
	output:
		feature_table = 'intermediate/{run}/aa/feature_table.biom'
	log: 'results/logs/{run}/collapse_feature_table_aa_unfiltered.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/collapse_feature_table_aa.py'

# perform greedy zero-width OTU clustering on AA ASVs
rule cluster_vhh_asv_pairs_aa:
	input: 'intermediate/{run}/aa/{library}_pairs_aa.csv'
	output:
		clusters = 'intermediate/{run}/aa_cluster/{library}_clusters_aa.csv',
		aa = 'intermediate/{run}/aa_cluster/{library}_aa_clustered.csv',
	threads: workflow.cores,
	log: 'results/logs/{run}/cluster_vhh_asv_pairs_aa/{library}.log'
	conda: 'envs/mmseqs2-vsearch.yaml'
	script: 'scripts/cluster_vhh_asv_pairs_aa.py'

# in case the same ASV is assigned to two different clusters in the alpaca and
# synthetic library, break the tie
rule resolve_aa_clusters:
	input: expand('intermediate/{{run}}/aa_cluster/{library}_clusters_aa.csv',library=get_libraries())
	output: 'intermediate/{run}/aa_cluster/clusters_aa.csv'
	log: 'results/logs/{run}/resolve_aa_clusters.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/resolve_clusters.py'

# group the feature table according to AA clustering
rule collapse_feature_table_aa_clusters:
	input:
		feature_table = 'intermediate/{run}/aa/feature_table.biom',
		mapping = 'intermediate/{run}/aa_cluster/clusters_aa.csv' #expand('intermediate/{{run}}/aa_cluster/{library}_clusters_aa.csv',library=get_libraries())
	output:
		feature_table = 'intermediate/{run}/aa_cluster/feature_table.biom'
	log: 'results/logs/{run}/collapse_feature_table_aa_clusters.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/collapse_feature_table_aa.py'

# align AA ASV sequences to the AA reference, so can extract CDRs and check
# lengths
rule align_vhh_asvs_aa:
	input:
		translation = 'intermediate/{run}/aa_cluster/{library}_aa_clustered.csv',
		reference = 'resources/references/{library}/{library}.fasta'
	output: 'intermediate/{run}/aa_align/{library}_aa_aligned.csv',
	params:
		library=lambda wc: get_library_params(wc.library)
	log: 'results/logs/{run}/align_vhh_asvs_aa/{library}.log'
	threads: workflow.cores
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/align_vhh_asvs_aa.py'

# extract CDRs to use in filtering
rule extract_cdrs_library:
	input: 'intermediate/{run}/aa_align/{library}_aa_aligned.csv'
	output: 'intermediate/{run}/aa_align/{library}_aa_cdrs.csv'
	conda: 'envs/biopython-pysam.yaml',
	log: 'results/logs/{run}/extract_cdrs_library/{library}.log'
	threads: workflow.cores
	params:
		CDRs=lambda wc: get_library_param(wc.library, 'CDRs')
	script: 'scripts/extract_cdrs.py'

# filter aaSVs according to critera, e.g. total lengh, length of CDRs, etc.
rule filter_vhh_asvs_aa:
	input:
		aa = 'intermediate/{run}/aa_align/{library}_aa_aligned.csv',
		cdrs = 'intermediate/{run}/aa_align/{library}_aa_cdrs.csv'
	output:
		aa = 'intermediate/{run}/aa_filter/{library}_aa.csv',
		cdrs = 'intermediate/{run}/aa_filter/{library}_aa_cdrs.csv',
	log: 'results/logs/{run}/filter_vhh_asvs_aa/{library}.log'
	params:
		library=lambda wc: get_library_params(wc.library)
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/filter_vhh_asvs_aa.py'

def input_concat_cdrs(wc):
	return {
		library: expand('intermediate/{run}/aa_filter/{library}_aa_cdrs.csv',
			library=library, run=wc.run)
		for library in get_libraries()
	}

rule concat_cdrs:
	#expand('intermediate/{run}/aa_filter/{library}_aa_cdrs.csv', library=get_libraries())
	input: unpack(input_concat_cdrs)
	output: 'intermediate/{run}/aa_filter/aa_cdrs.csv'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/concat_asvs.py'

# filters feature table to remove AAs filtered on basis of length, CDRs, etc.
rule collapse_feature_table_aa_filtered:
	input:
		feature_table = 'intermediate/{run}/aa_cluster/feature_table.biom',
		features = expand('intermediate/{{run}}/aa_filter/{library}_aa.csv',library=get_libraries()),
	output:
		feature_table = 'intermediate/{run}/aa_filter/feature_table.biom'
	log: 'results/logs/{run}/collapse_feature_table_aa_filtered.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/collapse_feature_table_aa.py'

rule collapse_feature_table_clonotype:
	input:
		feature_table = 'intermediate/{run}/aa_filter/feature_table.biom',
		mapping =       'intermediate/{run}/aa_filter/aa_cdrs.csv'
	output:
		feature_table = 'intermediate/{run}/cdrs/feature_table.biom'
	params:
		mapping_cols=['aaSVID','clonotypeID']
	log: 'results/logs/{run}/collapse_feature_table_clonotype.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/collapse_feature_table_aa.py'

use rule collapse_feature_table_clonotype as collapse_feature_table_cdr3 with:
	output:
		feature_table = 'intermediate/{run}/cdr3/feature_table.biom'
	params:
		mapping_cols=['aaSVID','CDR3ID']
	log: 'results/logs/{run}/collapse_feature_table_cdr3.log'


rule concat_feature_tables:
	input: expand('intermediate/{run}/aa_filter/feature_table.biom', run=get_runs())
	output: 'intermediate/aa/feature_table.biom'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/concat_feature_tables.py'

rule sum_feature_table_runs:
	input:
		feature_table='intermediate/aa/feature_table.biom',
		mapping='config/guids.tsv'
	output: 'results/tables/aa/feature_table.biom'
	params:
		mapping_cols=['guid','ID'],
		axis='sample'
	log: 'results/logs/sum_feature_table_runs.log'
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/collapse_feature_table_aa.py'

def input_concat_asvs(wc):
	return { library: expand(
		'intermediate/{run}/aa_filter/{library}_aa.csv',
			library=library, run=get_runs())
		for library in get_libraries() }

rule concat_asvs_aa:
	input: unpack(input_concat_asvs)
	output: 'results/tables/aa/asvs.csv'
	threads: workflow.cores/2
	conda: 'envs/biopython-pysam.yaml'
	script: 'scripts/concat_asvs.py'

# moved to common.smk
# --------------------------
rule extract_cdrs:
	input: 'results/tables/aa/asvs.csv'
	output: 'results/tables/aa/cdrs.csv'
	log: 'results/logs/aa/extract_cdrs.log'
	conda: 'envs/biopython-pysam.yaml',
	params:
		library_CDRs = {library: get_library_param(library, 'CDRs') for library in get_libraries() }
	script: 'scripts/extract_cdrs.py'

# moved to downstream.smk
# --------------------------
# rule link_asvs_for_mmseqs2:
# 	wildcard_constraints:
# 		space=list_to_regex(['na','aa','cdrs','cdr3'])
# 	input: 'results/tables/{space}/asvs.csv'
# 	output: 'intermediate/{space}/asvs.csv'
# 	shell: """
# 	cp "{input}" "{output}"
# 	"""

# rule asvs_to_mmseqs2db:
# 	input: 'intermediate/{space}/asvs.fasta'
# 	output: directory('intermediate/{space}/features_db')
# 	conda: 'envs/mmseqs2-vsearch.yaml'
# 	shell: """
# 	rm -r -f {output}
# 	mkdir -p {output}
# 	mmseqs createdb {input} {output}/features
# 	mmseqs createindex {output}/features {config[scratch]}
# 	"""

use rule collapse_feature_table_clonotype as collapse_feature_table_clonotype_all with:
	input:
		feature_table = 'results/tables/aa/feature_table.biom',
		mapping =       'results/tables/aa/cdrs.csv'
	output:
		feature_table = 'results/tables/cdrs/feature_table.biom'
	params:
		mapping_cols=['aaSVID','clonotypeID']
	log: 'results/logs/collapse_feature_table_clonotype.log'

	threads: workflow.cores//2

use rule collapse_feature_table_clonotype_all as collapse_feature_table_cdr3_all with:
	output:
		feature_table = 'results/tables/cdr3/feature_table.biom'
	params:
		mapping_cols=['aaSVID','CDR3ID']
	log: 'results/logs/collapse_feature_table_cdr3.log'



# use rule collapse_feature_table_clonotype as collapse_feature_table_cdr3_all with:
# 	input:
# 		feature_table = 'results/tables/aa/feature_table.biom',
# 		mapping =       'results/tables/aa/cdrs.csv'
# 	output:
# 		feature_table = 'results/tables/cdr3/feature_table.biom'
# 	params:
# 		mapping_cols=['aaSVID','CDR3ID']
# 	log: 'results/logs/collapse_feature_table_cdr3.log'
# 	threads: workflow.cores//2

rule copy_cdrs_clonotype:
	input: 'results/tables/aa/cdrs.csv'
	output: 'results/tables/{space}/asvs.csv'
	threads: workflow.cores//2
	script: 'scripts/copy_cdrs_clonotype.py'


rule dummy:
	input: 'results/tables/aa/feature_table.biom'
