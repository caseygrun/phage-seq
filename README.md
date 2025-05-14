# Code and Data for "Bacterial cell surface characterization by phage display coupled to high-throughput sequencing."

[![DOI](https://zenodo.org/badge/733679522.svg)](https://zenodo.org/doi/10.5281/zenodo.12863463)

Code used to process the high-throughput sequencing data in our manuscript "Bacterial cell surface characterization by phage display coupled to high-throughput sequencing." 

This repository contains several components:

- `nbseq-workflow` contains shared code for a Snakemake-based workflow for processing raw Phage-seq sequencing data. The other directories include symbolic links to the `scripts`, conda environments (`envs`), and `rules` defined in this directory. However, each specific experiment contains its own configuration, resources, and Snakemake workflow definitions (`Snakefile` and `*.smk`). Additionally, each experiment contains a set of Jupyter notebooks in `workflow/notebooks` which were used to generate figures in the manuscript.
	- `panning-small` contains code used to process the data from solid-phase and small-scale cell-based  phage display selection campaigns reported in Fig. 2
	- `panning-massive` contains code used to process the data from high-throughput cell-based phage display selection campaigns reported in Fig. 3D
	- `panning-extended` contains code used to process the data from the extended rounds of high-throughput cell-based phage display selection reported in Fig. 3E
	- `panning-minimal` is a self-contained, runnable example containing a small subset of the data from `panning-extended`. The code is otherwise the same as `panning-extended`.
- `alpaca-library` contains a Snakebake-based workflow for processing data where the entire input library was sequenced using longer reads in order to characterize its diversity.
- `other-figures` contains raw data and code used to generate other figures in the paper that do not rely on sequencing data (e.g. ELISAs and phagemid titers).

The [`nbseq`](https://github.com/caseygrun/nbseq) library is a close companion to this repository. 

## Installation and usage

1. Clone this entire repository:

	[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system or cluster, into the place where you want to perform the data analysis:

	```bash
	git clone https://github.com/caseygrun/phage-seq.git
	cd phage-seq
	```

2. Install the [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html#mamba-install) (or [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)) package manager if it is not already installed on your machine or cluster. I recommend using the [`miniforge`](https://github.com/conda-forge/miniforge) distribution. 
3. Install Snakemake using `conda` (or `mamba`):

	```bash
	conda create -c conda-forge -c bioconda -n snakemake snakemake mamba
	```

	This will create a conda environment called `snakemake`. All execution of the workflows will be done in this environment. For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

4. Prepare to execute the workflow

	Activate the conda environment:

	```bash
	conda activate snakemake
	```

	Change to the relevant directory, e.g. for `panning-minimal`:

	```bash
	cd panning-minimal
	```

5. Execute the workflow

	You may want to test your configuration first by performing a dry-run via:

	```bash
	snakemake --use-conda -n --cores all
	```

	This will print the execution plan but not run any of the code.

	Each workflow can then be executed by running this command from the relevant directory:

	```bash
	snakemake --use-conda --cores all
	```

	Alternatively, some experiments are broken into multiple Snakemake workflows (specified in distinct `.smk` files)
	which can be invoked individually. The workflows are arragned this way because the `preprocess` step creates lots 
	of large temporary files which you may not want to preserve on cluster storage; once this step is completed, the 
	`intermediate/` directory can be safely deleted and subsequent workflows invoked separately. The sub-workflows
	can be invoked using the `-s` flag to specify a path to a specific workflow file, e.g.:

	```bash
	snakemake --use-conda --cores all -s workflow/preprocess.smk
	snakemake --use-conda --cores all -s workflow/feature_table.smk
	snakemake --use-conda --cores all -s workflow/downstream.smk
	```

	See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for further details. In particular, you likely want to adjust `--cores $N` 
	to specify the number of CPU cores available on your machine or cluster instance.

## Included demonstrations

Several demonstrations of the workflow code and notebooks are available:

- **`panning-minimal`, a subset of the raw data from `panning-extended` is embedded in the repository** so that you can run the entire analysis pipeline, then interactively explore the results using the `nbseq` library. 

- **The full processed datasets for `panning-small`, `panning-massive`, `panning-extended`, and `alpaca-library` are available for download from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11246657).** Using these data packages, you do not need to run the workflow(s) but can explore the results interactively using the `nbseq` library. These datasets are suitable to generate the figures in the paper using the code in the  `workflow/notebooks` subdirectories. See detailed instructions below.

### `panning-minimal` self-contained demonstration
	
The `panning-minimal` dataset consists of six selections (three biological replicates each of two different conditions), plus several samples where the raw input library was sequenced without panning. Each sample has been arbitrarily downsamples to 7500 reads for the sake of file space and processing time.

To run the complete workflow, follow the steps above. The analysis requires 15--30 minutes to download and install conda packages, then ~10--15 min to run on a 16-core workstation.

After running the Snakemake workflow, launch a Jupyter Lab server as described in the [`nbseq` repository](http://github.com/caseygrun/nbseq), and navigate to the Jupyter notebook `panning-minimal/workflow/notebooks/analysis.ipynb`. 

Importantly, the first time you open this notebook, you will be prompted to choose the "notebook kernel;" this will dictate in which conda environment the Python process runs. Assuming you installed `nbseq` into an environment called `nbseq`, choose the entry titled "Python [conda env:nbseq]" and click "Select."

### `panning-extended` processed dataset
  
Download and extract the processed dataset from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11246657):

```bash
cd panning-extended
wget https://zenodo.org/records/12825488/files/panning-extended-results.tar.gz
tar vzxf panning-extended-results.tar.gz
```

Then, launch a Jupyter Lab server as described in the [`nbseq` repository](http://github.com/caseygrun/nbseq), and navigate to the Jupyter notebook [`panning-extended/workflow/notebooks/analysis.ipynb`](panning-extended/workflow/notebooks/analysis.ipynb). Importantly, the first time you open this notebook, you will be prompted to choose the "notebook kernel;" this will dictate in which conda environment the Python process runs. Assuming you installed `nbseq` into an environment called `nbseq`, choose the entry titled "Python [conda env:nbseq]" and click "Select." There are several other Jupyter notebooks within [`panning-extended/workflow/notebooks/`](panning-extended/workflow/notebooks/) which produce figures that appear in the manuscript and can be executed similarly (i.e. using the `nbseq` conda environment).

The data package will create `results/` and `intermediate/` subdirectories, which will contain the full amino acid and CDR3 feature tables, feature sequences, and a SQLite database for interactive exploration. All data necessary to execute the demonstration notebook is included.

Note that several additional files such as various transformed feature tables and beta-diversity calculations are omitted for the sake of simplicity and file size, 
but they can be regenerated by running the included snakemake workflow, `workflow/downstream.smk`. Likewise, the large `mmseqs2` database (which is needed to search the dataset for VHHs with similar sequences) is not included; it can be re-generated on-demand by running `snakemake --use-conda --cores all -s workflow/downstream.smk -- intermediate/cdr3/features_db/` for the CDR3 feature space or `snakemake --use-conda --cores all -s workflow/downstream.smk -- intermediate/aa/features_db/` for the amino acid feature space.

### `panning-massive` processed dataset
  
Download and extract the processed dataset from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11246657):

```bash
cd panning-massive
wget https://zenodo.org/records/12825488/files/panning-massive-results.tar.gz
tar vzxf panning-massive-results.tar.gz
```
The data package will create the `results/` subdirectory, containing all data necessary to execute the Jupyter notebooks in [`panning-massive/workflow/notebooks/`](panning-massive/workflow/notebooks/); these notebooks generate the figures that appear in the manuscript and can be executed as described above. 

Note that the notebook [`panning-massive/workflow/notebooks/fig-suppl-learning.ipynb`](panning-massive/workflow/notebooks/fig-suppl-learning.ipynb) requires creating a different conda environment, `nbseq-xgb`; this environment can be created by running:

```bash
cd panning-massive
conda env create -f envs/nbseq-xgb.yaml
```

### `panning-small` processed dataset
  
Download and extract the processed dataset from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11246657):

```bash
cd panning-small
wget https://zenodo.org/records/12825488/files/panning-small-results.tar.gz
tar vzxf panning-small-results.tar.gz
```

Then explore the notebooks in [`panning-small/workflow/notebooks/`](panning-small/workflow/notebooks/)

### `alpaca-library` processed dataset
  
Download and extract the processed dataset from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11246657):

```bash
cd alpaca-library
wget https://zenodo.org/records/12825488/files/alpaca-library-results.tar.gz
tar vzxf alpaca-library-results.tar.gz
```

Then explore the notebooks in [`alpaca-library/workflow/notebooks/`](alpaca-library/workflow/notebooks/).

### `other-figures`

`other-figures/data` includes the data necessary to run the notebooks within `other-figures/code`. These notebooks produce other figures in the paper.

## Adapting the workflow to future experiments

### Workflow data model

These workflows are designed to accommodate several realities of our experiments:

- Each experiment may involve multiple selection campaigns
- Each phage display experiment involves observing at least one population of VHHs multiple times over rounds of biopanning
- Each sample may be sequenced across multiple sequencing runs
- There may be technical replicate observations at the level of library preparation
- Different selections in the experiment may use entirely different starting phage-displayed VHH libraries with different consensus reference sequences

Key elements of the data model:

- A **selection** refers to a distinct population of phage-displayed VHHs which is enriched via multiple rounds of panning against a particular pair of antigen conditions (e.g. counter-selection and selection bacterial cells, antigen-negative and antigen-positive protein-coated wells, etc.); For example, a selection for flagella using a _∆fliC_ mutant as the antigen-negative cells and a wild-type strain as the antigen-positive cells occurring over 4 rounds. 
	- Each selection is assigned a single **"phage library"**; this mainly dictates to which _reference sequence_ reads should be aligned. For example, if selections in the experiment involved two different starting libraries---one from an alpaca immunization and one created synthetically---these would be considered different "phage libraries".
	- Multiple _biological replicates_ with the same  should generally be treated as _separate_ selections
- A **sample** is a single technical replicate observation of a population of phage/VHHs. For example: the post-amplification ("input") phage before round 2 of the above selection is one sample. Each sample must be prepared separately for sequencing (e.g. with a distinct barcode sequence of barcode pair) and demultiplexed outside of this workflow into a distinct pair of .fastq.gz files.
	- One _selection_ will generally have multiple _samples_; there will be at least one sample for each round of selection. If both the pre-amplification (output) and post-amplification (input) phage populations are sequenced at a given round, these will be separate samples. If phage are eluted from the antigen-negative cells and sequenced, this will be a separate sample. Additionally, technical replicates of the same population may be prepared. 
- Multiple **sequencing runs** may be performed. Each sequencing run will be **denoised** separately, then replicate observations of the same _sample_ will be summed in the `feature_table` step.


### Configuring the workflow for future experiments

`panning-minimal` can be duplicated and used as a template for future experiments. Follow the installation instructions above. Duplicate the `panning-minimal` directory and rename it to your project of interest. The following files need to be modified to suit the workflow. Note that all file paths are relative to this directory (e.g. whatever you rename your copy of `panning-minimal`):

1. Edit `config.yaml` to configure the experiment and specify features of the starting VHH library/libraries, namely: 

	- `raw_sequence_dir`: Path to folder containing raw input sequences. This folder must contain one subdirectory for each sequencing run, named according to `ID` in `runs.tsv`.

	- `scratch`: Path to a temporary or scratch directory

	- `input_query`: Pandas query <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html> to identify which rows in `metadata.csv` correspond to an un-panned "input" phage library. This is used downstream to build a null enrichment distribution and calculate enrichment probabilities. One way to accomplish this is to give the input library a distinct value for `expt` (e.g. set the `expt` column to `'input.library'` then set `input_query` to `expt == 'input.library'` ).

	- `libraries`: Describes properties of each distinct input VHH library included in this experiment. Input libraries are distinct if they originated from a different biological or synthetic source and thus have a different reference sequence. For example, libraries derived from different organisms or different synthetic methodolodies will be distinct. This object should have keys which are library names (corresponding to values of `phage_library` in `samples.tsv`, converted to lowercase) and values which are objects with the following attributes.

		- `primer_fwd` (default ""):
		Forward primer sequence (5' to 3'); used in the `cutadapt` preprocessing step to identify reads corresponding to properly-prepared amplicons
		
		- `primer_rev` (default ""):
		Reverse primer sequence (5' to 3'); used in the `cutadapt` preprocessing step to identify reads corresponding to properly-prepared amplicons

		- `reference` (default None):
		Path to the reference sequence, in FASTA format

		- `reference_frame_start_nt` (default 0):
		Nucleic acid position (0-based) indicating the first base of the first codon of the reference sequence
		
		- `reference_length_aa` (default 0):
		Length of the reference sequence in amino acids; if the reference sequence is longer than (`reference_frame_start_nt` + (`reference_length_aa` * 3)) nt, it will be trimmed.
		
		- `min_fwd_end` (default 0): 3' (distal) end of the forward read must align to this NA position or 
		later; position 0 is the first base of the reference sequence, irrespective of `reference_frame_start_nt`.

		- `max_rev_start` (default 0):
		3' (distal) end of reverse read must align to this NA position or 
		earlier; position 0 is the first base of the reference sequence, irrespective of `reference_frame_start_nt`.

		- `min_aa_length` (default 69): Reads where the aligned amino acid sequence (excluding gap characters) is shorter than this length will be dropped
		
		- `CDRs` (default {}): Position of CDR and FR regions within the reference sequence, in amino 
		acid space. In this object, position 0 refers to the amino acid 
		corresponding to  `reference_frame_start_nt`. This object should
		be a dict mapping domain names (e.g. 'CDR1', 'FR2', etc.) to `[start, end]`
		positions, where `start` is inclusive and `end` is exclusive (e.g. 
		half-open) intervals, following the Python convention.

			Example:

				'CDRs': {
					'FR1':  [0,  21],
					'CDR1': [21, 29],
					'FR2':  [29, 46],
					'CDR2': [46, 54],
					'FR3':  [54, 92],
					'CDR3': [92, 116],
					'FR4':  [116,127],
				}

		- `min_CDR_length` (default {}):
		Reads with domains (CDR or FR regions) shorter than this length (in amino acids)
		will be dropped. Dict where keys are domain names (e.g. 'CDR1', 'FR2', etc.
		and should correspond to domains defined in `CDRs`) and values are 
		minimum lengths (in amino acids)


2. Edit `metadata.csv` to specify parameters of each **selection**. Each row represents a single selection. Columns describe parameters of the selection, including the bacterial strains used for selection and counter-selection. The following columns are required. Additional columns can be added to identify phenotypes associated with the selection. 

	|     Column    |                                                               Description                                                               |
	|---------------|-----------------------------------------------------------------------------------------------------------------------------------------|
	| expt          | Sub-experiment name/number. Must correspond to the same column in `samples.tsv`                                                         |
	| selection     | Name for the selection; typically this is in the form `{plate}.{well}`, corresponding to the location within the biopanning microplate. |
	| antigen       | Which antigen is targeted by this selection                                                                                             |
	| background_CS | Genetic background for the counter-selection strain                                                                                     |
	| genotype_CS   | Genotype of the counter-selection strain, relative to the genetic background                                                            |
	| strain_CS     | Strain number or identifier of the counter-selection strain, if applicable                                                              |
	| cond_CS       | Growth condition of the counter-selection cells                                                                                         |
	| background_S  | Genetic background for the selection strain                                                                                             |
	| genotype_S    | Genotype of the selection strain, relative to the genetic background                                                                    |
	| strain_S      | Strain number or identifier of the selection strain, if applicable                                                                      |
	| cond_S        | Growth condition of the selection cells                                                                                                 |

3. Edit `phenotypes.csv` to provide information about the phenotypes, e.g. biological covariates of each selection denoted by additional columns of `metadata.csv`. Importantly, you may want to mark certain phenotypes as "antigens" for the sake of classifier training routines.

	|    Column   |                                  Description                                   |
	|-------------|--------------------------------------------------------------------------------|
	| name        | name of the phenotype; must correspond to one of the columns in `metadata.csv` |
	| type        | "antigen" or "phenotype"                                                       |
	| category    | (optional) the category of antigen or phenotype, e.g. "motility"               |
	| locus       | (optional) locus tag                                                           |
	| description | (optional) description of the antigen or phenotype                             |

4. Edit `samples.tsv` to list each **technical replicate**:

	
	|     Column    |                                                                                           Description                                                                                            |
	|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
	| ID            | sample identifier, generally `{plate}-{well}`; must correspond to a pair of raw sequence files described below                                                                                   |
	| expt          | sub-experiment number                                                                                                                                                                            |
	| sample        | **name of the selection**; must correspond to **selection** column in `metadata.csv` above                                                                                                       |
	| round         | selection round, in format `R{n}{io}` where `{n}` is replaced by the round number and `{io}` is replaced with `i` for input (e.g. post-amplification) or `o` for output (e.g. pre-amplification) |
	| phage_library | name of the starting library, for the sake of determining the applicable reference sequence and related parameters. Must correspond to a key within the `libraries` entry in `config.yaml`       |
	| plate         | (optional) HTS library preparation plate number                                                                                                                                                  |
	| well          | (optional) HTS library preparation well number                                                                                                                                                   |
	| depth         | (optional) expected relative sequencing depth                                                                                                                                                    |
	| notes         | (optional) notes about the selection                                                                                                                                                             |

5. Edit `runs.tsv` to specify details about each sequencing run. The only required column is `ID`; each run must have a corresponding directory in `${raw_sequence_dir}/${ID}` containing the demultiplexed sequencing data as described below. Additional columns may be used to record other data about the run, e.g. instrument, flow cell ID, date, etc.

6. Add the raw sequencing data to the input directory (the path specified by `raw_sequence_dir` in `config.yaml`; by default, this is a directory called `input`): 
	- There must be one subdirectory for each sequencing run, i.e.: `${raw_sequence_dir}/${sequencing_run_ID}`
	- Within each run subdirectory, there must be one folder for each sample, i.e. `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}`, where `${sequencing_run_ID}` corresponds to the column `ID` in `runs.tsv` above and `${sample_id}` corresponds to the `ID` column in `samples.tsv`.
	- Within that sample directory `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}`, there should be two files named:
		- `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}/*_R1_*.fastq.gz` containing the forward reads, and
		- `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}/*_R2_*.fastq.gz` containing the reverse reads 

7. Remove `guids.tsv`, `metadata_full.csv`, and `metadata_full.tsv` from the `config/` directory; these files will be regenerated by the workflow:

	```bash
	rm config/guids.tsv config/metadata_full.csv config/metadata_full.tsv
	```

Once you have prepared the configuration and input data, you can use the validation workflow to check your work:

```bash
snakemake --snakefile workflow/rules/validate.smk --directory .
```

This will give you a comprehensive report and describe how to correct any errors. 

### FAQ

0. What is the general organization of the workflow code?

	The analysis pipelines are written using [Snakemake](https://snakemake.github.io) workflows. A Snakemake workflow is composed of **rules** which describe how to convert input files into output files. The first "rule" of each workflow describes the target of the workflow---which output file(s) the workflow should try to create. The [`snakemake`](https://snakemake.github.io) executable reads the workflow and the configuration files, examines the available input (and output) file(s), and figures out how to create the output from the input by running rules. This process is designed to be **reproducible** (i.e. the workflow can be run again on the same input to yield the same output), **portable** (i.e. all code and dependencies are specified by the workflow repository), and **robust** (i.e. the workflow can be interrupted and can resume execution with minimal loss of time/progress).

	Each experiment from the paper has an analysis workflow in a separate subdirectory of this repository (`panning-small`, `panning-extended`, `panning-minimal`). Each of these subdirectories contains Snakemake workflows, configuration, and associated code:

		- Several Snakemake workflow definitions (`workflow/Snakefile` and `workflow/*.smk`). These are mostly the same between different experiments.
		- Configuration specific to that workflow in the `config/` subdirectory; files within specify which samples are included in the experiment, the structure of the reference sequence(s)
		- The `resources/` subdirectory, which contains the VHH reference sequence(s) to which reads are aligned
		- Jupyter Notebooks in `workflow/notebooks` which make plots that appear in the paper
		- Symbolic links within `workflow` to `scripts`, conda environment definitions (`envs`), and Snakemake `rules` in the top-level `nbseq-workflow` directory; this code is shared across multiple experiments. 

	Note that the input data is stored outside this repository, except for `panning-minimal` (see "[Included demonstrations](#included-demonstrations)" above).

	The Snakemake workflows are broken into several parts:

		- `workflow/preprocess.smk`: identifies and trims the primer sequences from reads, dereplicates and filters reads by quality, and denoises reads to remove sequencing errors using Dada2. This is the most computationally-intensive part of the workflow. The output is a **feature table** (a matrix of samples by VHH nucleic acid sequences) for each sequencing run and a list of VHH nucleic acid sequences. 
		- `workflow/feature_table.smk`: aligns these nucleic acid sequences to a reference sequence, translates them to amino acids, sums reads corresponding to identical or overlapping amino acid sequences, and similarly constructs feature tables which combine reads with identical clonotypes (CDR1--3 sequence) or CDR3 sequence
		- `workflow/downstream.smk`: performs various transformations on the feature tables, including calculating alpha diversity, pairwise sample distance/beta diversity, ordination (dimensionality reduction), multiple sequence alignment and phylogeny calculation of VHHs, and machine learning classification of selections.
		- `workflow/Snakefile` invokes each of the above sub-workflows. 
		
	It is recommended to run `workflow/Snakefile` directly to execute the entire workflow top-to-bottom (this is what happens by default when you run `snakemake` from within an experiment subdirectory). However, the sub-workflows are separated because it may be useful to archive or delete intermediate files after completing the `preprocess` or `feature_table` steps, in order to save storage space. The `feature_table.smk` or `downstream.smk` workflows can then be executed independently. .

1. **Where should the input data for the workflows be placed?**

	The input data should be placed in a subdirectory called `input` within the experiment directory (for example, `panning-extended/input`). Alternatively, the input data can be placed elsewhere and symbolic linked from `raw_sequence_dir`. 

	See section "[Included demonstrations](#included-demonstrations)" to obtain processed _output_ data (`results` and `intermediate`). 

2. **Where is the output of the workflow stored?**

	Output of the workflow is stored in three directories:

	- `config/`: some processing of the provided metadata (`samples.tsv`, `metadata.csv`, `runs.tsv`, `phenotypes.csv`) is done and stored in the files `metadata_full.csv`, `metadata_full.tsv` (these two have the same content), and `guids.tsv`. 
	- `intermediate/`: contains files resulting from intermediate steps in the workflow that are either not important for final interpretation *or* are fast/straightforward to re-generate from the file contents of `results/`
	- `results/`: contains the substantive output of the workflow
		- `logs/`: contains log files recording the output of each rule from the snakemake workflow. Generally, these are named the same as the rule name in the workflow. It's gonna be useful for debugging if a rule exits due to an error.
		- `runs/`: contains one subdirectory for each sequencing run. The denoising is first performed separately for each sequencing run, because the error characteristics/batch effects may differ across sequencing runs; later, feature tables from each sequencing run are summed. 
		- `plots/`: contains plots generated by the workflow and/or the Jupyter notebooks in `workflow/notebooks`
		- `tables/`: contains the feature tables and feature sequences themselves. Contains sub directories for each of 3 "feature spaces;" the "feature space" describes what the columns of the feature table mean (the rows are generally samples):
			- `aa`: columns represent distinct VHH amino acid sequences; reads from any columns in the original feature table which translate to the same amino acid sequence are summed.
			- `cdrs`: columns represent VHHs with distinct CDR1, 2, and 3 sequence; columns in the `aa` feature table which differ only in the framework sequences (FR1, 2, 3, 4) are summed.
			- `cdr3`: columns represent VHHs with distinct CDR3 sequences

3. **I ran `snakemake` non-interactively  (e.g. as a batch job) and I can't tell whether it completed successfuly.**

	First, check whether your batch job is still running (e.g. if using slurm, you can use `squeue`)! If it is not...

	Second, check the log files created by `snakemake`; these essentially have the same content as is printed to stdout while `snakemake` is running. Within the working directory for some experiment (e.g. `panning-minimal`), these files are in `.snakemake/log` (e.g. `panning-minimal/.snakemake/log`; note the `.` in `.snakemake`!). Scroll to the bottom of that file to give you some clue whether snakemake completed successfully or exited due to an error. 

	Third, whether or not there was obviously an error, check the rule-specific log(s) for the last rule or rules to run. Each rule generally has a log in `results/logs/{rule}.log` or similar; this path should also be given in the main snakemake log.

	Third, try running `snakemake --dry-run` again. If the output looks reasonable, then run `snakemake`; if it continues successfully, problem solved. If you get a message about the working directory being locked (this happens sometimes if your job is killed by your HPC cluster scheduler), you may need to run `snakemake --unlock`. If `snakemake` died in the middle of a task, you may get a message about how the output files are incomplete; run `snakemake --rerun-incomplete`. 
4. I've run part of the workflow, but now `snakemake` wants to re-run a bunch of rules, perhaps even the whole workflow.

	Snakemake will, by default, re-run a rule _and all rules that depend on that rule's output_ if any of the following have changed since the rule was last run:

	- Input files for that rule (modification time or `mtime` specifically has changed, even if file content has not)
	- Code for that rule
	- Any parameters (`params`) to the rule
	- The `conda` environment
	
	This may have happened because you changed any of the above. Moving files usually doesn't trigger this, but moving `conda` environments can. You can approach this by changing `--rerun-triggers mtime`


