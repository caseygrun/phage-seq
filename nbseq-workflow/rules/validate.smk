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

def print_color(color, *args, **kwargs):
	# text = "\n".join(args)
	import itertools
	text = "\n".join(itertools.chain(*(textwrap.wrap(arg, width=80, subsequent_indent='  ', **kwargs) for arg in args)))
	print(color + text + bcolors.ENDC)

def print_info(*args):
	print_color(bcolors.OKCYAN, *args, initial_indent=u'  ',)

def print_good(*args):
	print_color(bcolors.OKGREEN, *args, initial_indent=u'\u2713 ',)

def print_warn(*args):
	print_color(bcolors.WARNING, *args, initial_indent=u'! ',)


# load and validate configuration
# (config is actually loaded from rules/common.smk)
# configfile: 'config/config.yaml'
print_info('Validating `config/config.yaml`...')
validate(config, "schemas/config.schema.yaml")

# TODO: check that all properties in `min_CDR_length` are in `CDRs`
for name, library in config['libraries'].items():
	cdrs = set(library['CDRs'].keys())
	for cdr in library['min_CDR_length']:
		assert cdr in cdrs,(f"In config.yaml, in library '{library}', a minimum length ('min_CDR_length') "
				f"is given for domain '{cdr}', but the boundaries of this domain are not given in "
				f"'CDRs' (domains defined: {list(cdrs)}).")

print_good('config.yaml is valid')


# validate samples.tsv
print_info('Validating `config/samples.tsv`...')

# load samples.tsv separately, since in common.smk it is augmented with run, guid
SAMPLES_TSV = pd.read_csv('config/samples.tsv', sep="\t")

# check that there are no empty (all null) rows, as this is a common source of error
rows_of_all_na = SAMPLES_TSV.isnull().all(1)
assert sum(rows_of_all_na) == 0, f"There the following rows in samples.tsv are empty: {SAMPLES_TSV.index[rows_of_all_na]}. Microsoft Excel can add empty rows to a CSV file; try opening in a plain text editor to remove these empty rows."

validate(SAMPLES_TSV, "schemas/samples.schema.yaml")

# check that each phage_library corresponds to config['libraries']
assert SAMPLES["phage_library"].str.lower().isin(config['libraries'].keys()).all(), "each entry in column `phage_library` of `samples.tsv` must correspond to a key in `libraries` within `config.yaml`"

print_good('samples.tsv is valid;', f"{len(SAMPLES)} samples: {list(SAMPLES.ID)}")



print_info('Validating `config/metadata.csv`...')
metadata = pd.read_csv('config/metadata.csv')

# check that there are no empty (all null) rows, as this is a common source of error
rows_of_all_na = metadata.isnull().all(1)
assert sum(rows_of_all_na) == 0, f"There the following rows in metadata.csv are empty: {metadata.index[rows_of_all_na]}. Microsoft Excel can add empty rows to a CSV file; try opening in a plain text editor to remove these empty rows."


validate(metadata, "schemas/metadata.schema.yaml")
with open("workflow/schemas/metadata.schema.yaml", 'r') as f:
	metadata_schema = yaml.safe_load(f)


# check that expts correspond between samples.tsv and metadata.csv
samples_expts = set(SAMPLES['expt'].unique())
metadata_expts = set(metadata['expt'].unique())
if (samples_expts != metadata_expts):
	samples_only_expts = samples_expts.difference(metadata_expts)
	metadata_only_expts = metadata_expts.difference(samples_expts)

	warning = ["Warning: values of `expt` in `samples.tsv` and `metadata.csv` do not match:"]
	if len(samples_only_expts) > 0:
		warning = warning + [
			f"`expt`s only in `samples.tsv`: {samples_only_expts}",
			(f"This is not an error, but samples with `expt` in {samples_only_expts} will not gain any associated metadata from metadata.csv. "
				"Please check that this is what you intend.")
		]

	if len(metadata_only_expts) > 0:
		warning = warning + [
			f"`expt`s only in `metadata.csv`: {metadata_only_expts}",
			(f"This is not an error, but no samples from selections in `expt`(s) {metadata_only_expts} have been processed."
				"Please check that this is what you intend.")
		]
	print_warn(*warning)

	# TODO: check that, for each sample in samples.tsv, expt matches that given in metdata.csv

# check that metadata.selection corresponds to samples.sample
samples_selections = set(SAMPLES['sample'].unique())
metadata_selections = set(metadata['selection'].unique())
if (samples_selections != metadata_selections):
	if len(samples_selections - metadata_selections) > 0:
		bad_selections = samples_selections - metadata_selections
		bad_samples = SAMPLES.loc[SAMPLES['sample'].isin(bad_selections),:]

		print_warn(
			"Warning: there are entries in column `sample` in `config/samples.tsv` that do not have a corresponding entry in column `selection` in `config/metadata.csv`. ",
			"Each sample defined in `config/samples.tsv` should have a corresponding selection defined in `config/metadata.csv`.",
			"This is not an error, but the samples below will not gain any metadata from `metadata.csv`",
			"Samples with no selection metadata (from `config/samples.tsv`):"
			f"{bad_samples}")
	elif len(metadata_selections - samples_selections) > 0:
		bad_selections = metadata_selections - samples_selections
		bad_metadata = metadata.loc[metadata['selection'].isin(bad_selections),:]
		print_warn(
			"There are entries in column `selection` in `config/metadata.csv` that do not have a corresponding entry in column `sample` in `config/samples.tsv`. "
			"Each selection defined in `config/metadata.csv` must be represented by at least one sample in `config/samples.tsv`."
			"This is not an error, but the selections below do not have any associated samples."
			"Selections (from `config/metadata.csv`) with no associated samples:"
			f"{bad_metadata}")

print_good('metadata.csv is valid', f"{len(metadata)} selections: {list(metadata.selection)}")


print_info("Validating `config/phenotypes`.csv...")
PHENOTYPES = pd.read_csv('config/phenotypes.csv')
validate(PHENOTYPES, "schemas/phenotypes.schema.yaml")

# check if all extra metadata columns correspond to phenotypes
metadata_predefined_columns = metadata_schema['properties'].keys()
metadata_phenotype_columns = set(metadata.columns) - set(metadata_predefined_columns)
metadata_unexpected_columns = set(metadata_phenotype_columns) - set(PHENOTYPES.name)
if len(metadata_unexpected_columns) > 0:
	print_warn(
		"There are additional columns in `metadata.csv` that are not defined in phenotypes.csv.",
		f"- Phenotypes (from phenotyes.csv): {list(PHENOTYPES.name)}",
		f"- Unexpected metadata columns: {metadata_unexpected_columns}",
		"This is not an error, but it is strongly recommended that phenotype columns be defined in phenotypes.csv"
		)

print_good('phenotypes.csv is valid', f"{len(PHENOTYPES)} phenotypes: {list(phenotypes.name)}")

print_info("Validating `config/runs.tsv`...")

# runs.tsv has ID column
assert 'ID' in RUNS, "`runs.tsv` must contain a non-empty column `ID`"

# run IDs are all non-empty
assert not (RUNS['ID'].isna().any()),  "Column `ID` in runs.tsv must not be empty"
assert not ((RUNS['ID'] == '').any()), "Column `ID` in runs.tsv must not be empty"

print_good('runs.tsv is valid', f"{len(RUNS)} phenotypes: {list(RUNS.ID)}")



print_info("Checking input files...")

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


print_good(f'configuration and input data are valid!')

rule dummy:
	run: 
		pass