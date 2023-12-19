
from pathlib import Path
import pandas as pd
from scripts.common import *

configfile: 'config/config.yaml'


def get_library_params(library):
	if library in config['libraries']:
		return config['libraries'][library]
	else:
		raise Exception(f"Unknown library {library} (is it listed in config.yaml under `libraries`?)")

def get_library_param(library, param):
	return get_library_params(library)[param]

def get_libraries():
	return config['libraries'].keys()

def get_sample_library_params(run, sample):
	library = SAMPLES.loc[SAMPLES['guid'] == sample, 'phage_library'].iat[0].lower()
	return get_library_params(library)

def get_sample_library_param(run, sample, param):
	return get_sample_library_params(run, sample)[param]

def get_samples_per_run(run):
	return SAMPLES.loc[SAMPLES['run_id'] == run,'guid']

def get_runs():
	return RUNS.ID.astype(str)

def get_metadata_values(columns):
	return SAMPLES.loc[:,columns].drop_duplicates()




RUNS = pd.read_csv('config/runs.tsv', sep="\t")
SAMPLES = expand_metadata_by_runs(pd.read_csv('config/samples.tsv', sep="\t"), RUNS)
ANTIGENS = pd.read_csv('config/phenotypes.csv').query("type == 'antigen'")



rule summarize_feature_table:
	wildcard_constraints:
		kind=list_to_regex(['nochimeras', 'aa', 'aa_cluster', 'aa_filter', 'cdrs', 'cdr3']),
		run=list_to_regex(get_runs())
	input:
		feature_table = 'intermediate/{run}/{kind}/feature_table.biom'
	output:
		reads          = 'intermediate/{run}/{kind}/summary_reads.txt',
		features       = 'intermediate/{run}/{kind}/summary_features.txt',
		total_features = 'intermediate/{run}/{kind}/summary_total_features.txt'
	params:
		header=lambda wc: wc.kind,
		total_row=lambda wc: wc.run
	log: 'results/logs/{run}/{kind}_summarize_feature_table.log'
	conda: '../envs/biopython-pysam.yaml'
	script: '../scripts/summarize_feature_table.py'

rule summarize_feature_table_checkpoint:
	wildcard_constraints:
		kind=list_to_regex(['denoise', 'align']),
		run=list_to_regex(get_runs())
	input:
		feature_table = 'results/runs/{run}/{kind}/feature_table.biom'
	output:
		reads          = 'intermediate/{run}/{kind}/summary_reads.txt',
		features       = 'intermediate/{run}/{kind}/summary_features.txt',
		total_features = 'intermediate/{run}/{kind}/summary_total_features.txt'
	params:
		header=lambda wc: wc.kind,
		total_row=lambda wc: wc.run
	log: 'results/logs/{run}/{kind}_summarize_feature_table.log'
	conda: '../envs/biopython-pysam.yaml'
	script: '../scripts/summarize_feature_table.py'


rule asvs_csv_to_fasta:
	input:  'intermediate/{subdirs}/asvs.csv'
	output: 'intermediate/{subdirs}/asvs.fasta'
	conda: '../envs/biopython-pysam.yaml'
	script: '../scripts/csv_to_fasta.py'
