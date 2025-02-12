from pathlib import Path
import yaml
import textwrap

import pandas as pd
from snakemake.utils import validate


include: 'rules/common.smk'

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_good(*args):
	text = "\n".join(args)
	text = textwrap.wrap(text, width=70, *, initial_indent=u'\u2713 ', subsequent_indent='  ')
	print(bcolors.OKGREEN + text + bcolors.ENDC)

def print_warn(*args):
	text = "\n".join(args)
	text = textwrap.wrap(text, width=70, *, initial_indent=u'\u2713 ', subsequent_indent='  ')
	print(bcolors.WARNING + text + bcolors.ENDC)


# load and validate configuration
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")

print_good('config.yaml is valid')

# validate samples.tsv
validate(SAMPLES, "schemas/samples.schema.yaml")

# phage_library corresponds to config['libraries']
assert SAMPLES["phage_library"].str.lower().isin(config['libraries'].keys()).all(), "each entry in column `phage_library` of `samples.tsv` must correspond to a key in `libraries` within `config.yaml`"

print_good('samples.tsv is valid', f"{len(SAMPLES)} samples: {list(SAMPLES.ID)}")


metadata = pd.read_csv('config/metadata.csv')
validate(metadata, "schemas/metadata.schema.yaml")
metadata_schema = yaml.load("schemas/metadata.schema.yaml")


# expt corresponds between samples.tsv and metadata.csv
samples_expts = set(SAMPLES['expt'].unique())
metadata_expts = set(metadata['expt'].unique())
if (samples_expts != metadata_expts):
	samples_only_expts = samples_expts.difference(metadata_expts)
	metadata_only_expts = metadata_expts.difference(samples_expts)

	raise Exception("Values of `expt` in `samples.tsv` and `metadata.csv` do not match:"
			"`expt`s only in `samples.tsv`: {samples_only_expts}"
			"`expt`s only in `metadata.csv`: {metadata_only_expts}"
			)

	# TODO: check that, for each sample in samples.tsv, expt matches that given in metdata.csv

# metadata.selection corresponds to samples.sample
samples_selections = set(SAMPLES['sample'].unique())
metadata_selections = set(metadata['selection'].unique())
if (samples_selections != metadata_selections):
	if len(samples_selections - metadata_selections) > 0:
		bad_selections = samples_selections - metadata_selections
		bad_samples = SAMPLES.loc[SAMPLES['sample'].isin(bad_selections),:]

		raise Exception(
			"There are entries in column `sample` in `samples.tsv` that do not have a corresponding entry in column `selection` in `metadata.csv`. Each sample defined in `samples.tsv` must have a corresponding selection defined in `metadata.csv`."
			"Bad samples (from `samples.tsv`):\n"
			f"{bad_samples}")
	elif len(metadata_selections - samples_selections) > 0:
		bad_selections = metadata_selections - samples_selections
		bad_metadata = metadata.loc[metadata['selection'].isin(bad_selections),:]
		raise Exception(
			"There are entries in column `selection` in `metadata.csv` that do not have a corresponding entry in column `sample` in `samples.tsv`. Each selection defined in `metadata.csv` must be represented by at least one sample in `samples.tsv`."
			"Bad metadata entries (from `metadata.csv`):\n"
			f"{bad_metadata}")

print_good('metadata.csv is valid', f"{len(metadata)} selections: {list(metadata.selections)}")


PHENOTYPES = pd.read_csv('config/phenotypes.csv')
validate(PHENOTYPES, "schemas/samples.schema.yaml")

# check if all extra metadata columns correspond to phenotypes
metadata_predefined_columns = metadata_schema['properties'].keys()
metadata_phenotype_columns = set(metadata.columns) - set(metadata_predefined_columns)
metadata_unexpected_columns = set(metadata_phenotype_columns) - set(PHENOTYPES.name)
if len(metadata_unexpected_columns) > 0:
	print_warn(
		"There are additional columns in `metadata.csv` that are not defined in phenotypes.csv.",
		f"- Phenotypes (from phenotyes.csv): {list(phenotypes.names)}",
		f"- Unexpected metadata columns: {metadata_unexpected_columns}",
		"This is not an error, but it is recommended that phenotype columns be defined in phenotypes.csv"
		)

print_good('phenotypes.csv is valid', f"{len(PHENOTYPES)} phenotypes: {list(phenotypes.name)}")

# runs.tsv has ID column
assert 'ID' in RUNS, "`runs.tsv` must contain a non-empty column `ID`"

# run IDs are all non-empty
assert not (RUNS['ID'].isna().any()),  "Column `ID` in runs.tsv must not be empty"
assert not ((RUNS['ID'] == '').any()), "Column `ID` in runs.tsv must not be empty"

print_good('runs.tsv is valid', f"{len(RUNS)} phenotypes: {list(RUNS.ID)}")


# check that raw_sequence_dir exists
raw_sequence_dir = Path(config['raw_sequence_dir'])
assert raw_sequence_dir.exists(), f"config['raw_sequence_dir'] = {raw_sequence_dir} does not exist"
assert raw_sequence_dir.is_dir(), f"config['raw_sequence_dir'] = {raw_sequence_dir} is a file, not a directory"

# check that there is a run data directory for each file
for run in RUNS['ID']:
	run_sequence_dir = raw_sequence_dir / run

	assert raw_sequence_dir.exists(), f"for run {run}, expected directory {run_sequence_dir} does not exist"
	assert raw_sequence_dir.is_dir(), f"for run {run}, expected directory {run_sequence_dir} is a file, not a directory"

	# check that sample subdirectories exist
	for sample_id in SAMPLES.ID:
		sample_dir = run_sequence_dir / f"Sample_{sample_id}"

		# TODO: currently checks that a sample directory exists for each run; I believe this is currently required but technically need not be the case 
		assert sample_dir.exists(), f"for run {run} and sample {sample_id}, expected directory {sample_dir} does not exist"
		assert sample_dir.is_dir(), f"for run {run} and sample {sample_id}, expected directory {sample_dir} is a file, not a directory"

		# TODO: glob to check that raw sample .fastq.gz files actually exist

print_good(f'raw sequences were found in {raw_sequence_dir}')



guids = guids_to_samples(samples, runs)


rule dummy:
	run: 
		pass