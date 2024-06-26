Code used to process the high-throughput sequencing data in our manuscript "Bacterial cell surface characterization by phage display coupled to high-throughput sequencing." 

This repository contains five components:

- `nbseq-workflow` contains shared code for a Snakemake-based workflow for processing raw Phage-seq sequencing data. The other directories include symbolic links to the `scripts`, conda environments (`envs`), and `rules` defined in this directory. However, each specific experiment contains its own configuration, resources, and Snakemake workflow definitions (`Snakefile` and `*.smk`).
- `panning-small` contains code used to process the data from solid-phase and small-scale cell-based  phage display selection campaigns reported in Fig. 2
- `panning-massive` contains code used to process the data from high-throughput cell-based phage display selection campaigns reported in Fig. 3D
- `panning-extended` contains code used to process the data from the extended rounds of high-throughput cell-based phage display selection reported in Fig. 3E
- `panning-minimal` is a runnable example containing a small subset of the data from `panning-extended`. The code is otherwise the same as `panning-extended`.


## Installation and usage

1. Clone this entire repository

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

	Each workflow can be executed by running this command from the relevant directory:

	```bash
	snakemake --use-conda --cores all
	```

	This will print the execution plan but not run any of the code.

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

Two demonstrations of the workflow code are included:

- **`panning-minimal`: a subset of the raw data from `panning-extended` is embedded in the repository** so that you can run the entire analysis pipeline, then interactively explore the results using the `nbseq` library. 
	
	The `panning-minimal` dataset consists of six selections (three biological replicates each of two different conditions), plus several samples where the raw input library was sequenced without panning. Each sample has been arbitrarily downsamples to 7500 reads for the sake of file space and processing time.

	To run the complete workflow, follow the steps above. The analysis requires 15--30 minutes to download and install conda packages, then ~10--15 min to run on a 16-core workstation.

	After running the Snakemake workflow, launch a Jupyter Lab server as described in the [`nbseq` repository](http://github.com/caseygrun/nbseq), and navigate to the Jupyter notebook `panning-minimal/workflow/notebooks/analysis.ipynb`. 

	Importantly, the first time you open this notebook, you will be prompted to choose the "notebook kernel;" this will dictate in which conda environment the Python process runs. Assuming you installed `nbseq` into an environment called `nbseq`, choose the entry titled "Python [conda env:nbseq]" and click "Select."

- **The full processed dataset for `panning-extended` is available for exploration.** You do not need to run the workflow but can explore the results interactively using the `nbseq` library. 
  
	Download and extract the processed dataset from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11246657):

	```bash
	cd panning-extended
	wget https://zenodo.org/records/11246658/files/panning-extended-results.tar.gz
	tar vzxf panning-extended-results.tar.gz
	```

	Then, launch a Jupyter Lab server as described in the [`nbseq` repository](http://github.com/caseygrun/nbseq), and navigate to the Jupyter notebook `panning-minimal/workflow/notebooks/analysis.ipynb`. Importantly, the first time you open this notebook, you will be prompted to choose the "notebook kernel;" this will dictate in which conda environment the Python process runs. Assuming you installed `nbseq` into an environment called `nbseq`, choose the entry titled "Python [conda env:nbseq]" and click "Select."

	The data package will create `results/` and `intermediate/` subdirectories, which will contain the full amino acid and CDR3 feature tables, feature sequences, and a SQLite database for interactive exploration. All data necessary to execute the demonstration notebook is included.

	Note that several additional files such as various transformed feature tables and beta-diversity calculations are omitted for the sake of simplicity and file size, 
	but they can be regenerated by running the included snakemake workflow, `workflow/downstream.smk`. Likewise, the large `mmseqs2` database (which is needed to search the dataset for VHHs with similar sequences) is not included; it can be re-generated on-demand by running `snakemake --use-conda --cores all -s workflow/downstream.smk -- intermediate/cdr3/features_db/` for the CDR3 feature space or `snakemake --use-conda --cores all -s workflow/downstream.smk -- intermediate/aa/features_db/` for the amino acid feature space.


## Adapting the workflow to future experiments

### Workflow data model

These workflows are designed to accommodate several realities of our experiments:

- Each experiment may involve multiple selection campaigns
- Each phage display experiment involves observing at least one population of VHHs multiple times over rounds of biopanning
- Each sample may be sequenced across multiple sequencing runs
- There may be technical replicate observations at the level of library preparation
- Different selections in the experiment may use entirely different starting phage-displayed VHH libraries with different consensus reference sequences

Key elements of the data mode

- A **selection** refers to enrichment of a distinct population of phage-displayed VHHs via multiple rounds of panning against a particular pair of antigen conditions (e.g. counter-selection and selection bacterial cells, antigen-negative and antigen-positive protein-coated wells, etc.); For example, a selection for flagella using a _∆fliC_ mutant as the counter-selection cells and a wild-type strain as the selection cells occurring over 4 rounds of selection. 
	- Each selection is assigned a single **"phage library"**; this mainly dictates to which _reference sequence_ reads should be aligned. For example, if selections in the experiment involved two different starting libraries---one from an alpaca immunization and one created synthetically---these would be considered different "phage libraries".
- A **sample** is a single technical replicate observation of a population of phage/VHHs, for example: the post-amplification ("input") phage before round 2 of the above selection, prepared for sequencing with a distinct pair of barcode sequences. 
- Multiple **sequencing runs** may be performed. Each sequencing run will be **denoised** separately, then replicate observations of the same _sample_ will be summed in the `feature_table` step.


### Configuring the workflow for future experiments

`panning-extended` can be duplicated and used as a template for future experiments. Follow the installation instructions above. Duplicate the `panning-extended` directory and rename it. The following files need to be modified to suit the workflow:

1. Edit `config.yaml` to configure the experiment and specify features of the starting VHH library/libraries, namely: 

	- `raw_sequence_dir`: Path to folder containing raw input sequences

	- `scratch`: Path to scratch directory

	- `libraries`: Phage display (VHH) libraries included in this experiment. Dict where keys are
	library names and values are objects with the following attributes:


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

	|     Column    |                                                              Description                                                               |
	|---------------|----------------------------------------------------------------------------------------------------------------------------------------|
	| expt          | Sub-experiment                                                                                                                         |
	| selection     | Name for the selection; typically this is in the form `{plate}.{well}`, corresponding to the location within the biopanning microplate |
	| antigen       | Which antigen is targeted by this selection                                                                                            |
	| background_CS | Genetic background for the counter-selection strain                                                                                    |
	| genotype_CS   | Genotype of the counter-selection strain, relative to the genetic background                                                           |
	| strain_CS     | Strain number or identifier of the counter-selection strain, if applicable                                                             |
	| cond_CS       | Growth condition of the counter-selection cells                                                                                        |
	| background_S  | Genetic background for the selection strain                                                                                            |
	| genotype_S    | Genotype of the selection strain, relative to the genetic background                                                                   |
	| strain_S      | Strain number or identifier of the selection strain, if applicable                                                                     |
	| cond_S        | Growth condition of the selection cells                                                                                                |

3. Edit `phenotypes.csv` to provide information about the phenotypes, e.g. biological covariates of each selection denoted by additional columns of `metadata.csv`. Importantly, you may want to mark certain phenotypes as "antigens" for the sake of classifier training routines.

	|    Column   |                                   Description                                    |
	|-------------|----------------------------------------------------------------------------------|
	| name        | name of the phenotype; should correspond to one of the columns in `metadata.csv` |
	| type        | "antigen" or "phenotype"                                                         |
	| category    | (optional) the category of antigen or phenotype, e.g. "motility"                 |
	| locus       | (optional) locus tag                                                             |
	| description | (optional) description of the antigen or phenotype                               |

4. Edit `samples.tsv` to list each **technical replicate**:

	
	|     Column    |                                                                                           Description                                                                                            |
	|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
	| ID            | sample identifier, generally `{plate}-{well}`; should correspond to a pair of raw sequence files described below                                                                                 |
	| expt          | sub-experiment number                                                                                                                                                                            |
	| sample        | **name of the selection**, corresponds to **selection** column in `metadata.csv` above                                                                                                           |
	| round         | selection round, in format `R{n}{io}` where `{n}` is replaced by the round number and `{io}` is replaced with `i` for input (e.g. post-amplification) or `o` for output (e.g. pre-amplification) |
	| phage_library | name of the starting library, for the sake of determining the applicable reference sequence and related parameters. Must correspond to a key within the `libraries` entry in `config.yaml`       |
	| plate         | (optional) HTS library preparation plate number                                                                                                                                                  |
	| well          | (optional) HTS library preparation well number                                                                                                                                                   |
	| depth         | (optional) expected relative sequencing depth                                                                                                                                                    |
	| notes         | (optional) notes about the selection                                                                                                                                                             |

5. Edit `runs.tsv` to specify details about each sequencing run. The only required column is `ID`; each run must have a corresponding directory in `${raw_sequence_dir}/${ID}` containing the demultiplexed sequencing data as described below. Additional columns may be used to record other data about the run, e.g. instrument, flow cell ID, date, etc.

6. Add the raw data to the path `raw_sequence_dir` specified in `config.yaml` (by default, this is a directory called `input`): 
	- There should be a subdirectory for each sequencing run, i.e.: `${raw_sequence_dir}/${sequencing_run_ID}`
	- There should be a folder for each sample, i.e. `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}`, where `${sequencing_run_ID}` corresponds to the column `ID` in `runs.tsv` above
	- Within that sample, there should be two files named `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}/*_R1_*.fastq.gz` and `${raw_sequence_dir}/${sequencing_run_ID}/Sample_${sample_id}/*_R2_*.fastq.gz` containing the forward and reverse reads respectively. `${sample_id}` corresponds to the column `ID` in `samples.tsv` above



