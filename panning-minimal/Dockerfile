FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="044d187d03ecd9b7df7ad659a76b2e37019c92c1f8e3ca374de196f2da5c7ee4"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/beta-diversity.yaml
#   prefix: /conda-envs/12118d3115c8565e50c5fda969d0a5bf
#   # code/environments/exact/cc085c6e7cf9387b6cc924cb88c5b0b9_.yaml
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - conda-forge::deicode=0.2.4
#   - scikit-bio=0.5.8
#   - unifrac=1.2
#   - biom-format=2.1.14
#   - biopython=1.80
#   - joblib=1.2.0
#   - scipy=1.10.0
#   - pandas=1.5.2
#   - openpyxl=3.0.10
#   - anndata=0.8.0
#   - scikit-learn=1.2.0
#   - conda-forge::simplejson=3.18.1
#   - pip
#   - pip:
#     - "git+https://github.com/caseygrun/nbseq.git#egg=nbseq"
RUN mkdir -p /conda-envs/12118d3115c8565e50c5fda969d0a5bf
COPY workflow/envs/beta-diversity.yaml /conda-envs/12118d3115c8565e50c5fda969d0a5bf/environment.yaml

# Conda environment:
#   source: workflow/envs/biopython-pysam.yaml
#   prefix: /conda-envs/3a8e66a7ad949aaf6a4d9cb01a91dea1
#   # code/environments/exact/8a0b04e2ac4c67c083cb820309af2500_.yaml
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - biopython=1.80
#   - pandas=1.5.2
#   - pysam=0.20.0
#   - scipy=1.9.3
#   - biom-format=2.1.13
#   - anndata=0.8.0
#   - joblib=1.2.0
#   - sqlite=3.40.0
#   - sqlalchemy=1.4.45
#   - pip
#   - pip:
#     - "git+https://github.com/caseygrun/nbseq.git#egg=nbseq"
RUN mkdir -p /conda-envs/3a8e66a7ad949aaf6a4d9cb01a91dea1
COPY workflow/envs/biopython-pysam.yaml /conda-envs/3a8e66a7ad949aaf6a4d9cb01a91dea1/environment.yaml

# Conda environment:
#   source: workflow/envs/bowtie2.yaml
#   prefix: /conda-envs/6a5eb90c6ff9a883e31a798150db79a0
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   # - bowtie2=2.3.4.3=py36h5c9b4e4_1
#   # - tbb=2019.3=h04f5b5a_1000
#   # - zlib=1.2.11=h1de35cc_1004
#   # - perl=5.32.0=hbcb3906_0
#   - bowtie2=2.3.4.3
#   - tbb=2019.3
#   - zlib=1.2.11
#   - perl=5.32.0
#   - samtools=1.15.1
#   - biopython=1.79
#   - pandas=1.3.5
#   - python>=3.7
#   - pip
#   - pip:
#     - "git+https://github.com/caseygrun/nbseq.git#egg=nbseq"
RUN mkdir -p /conda-envs/6a5eb90c6ff9a883e31a798150db79a0
COPY workflow/envs/bowtie2.yaml /conda-envs/6a5eb90c6ff9a883e31a798150db79a0/environment.yaml

# Conda environment:
#   source: workflow/envs/breakaway.yaml
#   prefix: /conda-envs/8448d2cbdccdc2abdf724bf1dcbe400f
#   # a124dc4f299b4b25632a0b82d522639b.yaml
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - bioconductor-rhdf5=2.38.0
#   - r-breakaway=3.0
#   - r-ggplot2=3.3.6
#   - r-lobstr=1.1.2
#   - r-matrix=1.4_1
#   - r-pbapply=1.5_0
#   - r-tidyverse=1.3.2
RUN mkdir -p /conda-envs/8448d2cbdccdc2abdf724bf1dcbe400f
COPY workflow/envs/breakaway.yaml /conda-envs/8448d2cbdccdc2abdf724bf1dcbe400f/environment.yaml

# Conda environment:
#   source: workflow/envs/cutadapt.yaml
#   prefix: /conda-envs/dd444f859e73c5e562ebcf217729c29c
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - cutadapt=2.3
#   - pandas=1.3.5
RUN mkdir -p /conda-envs/dd444f859e73c5e562ebcf217729c29c
COPY workflow/envs/cutadapt.yaml /conda-envs/dd444f859e73c5e562ebcf217729c29c/environment.yaml

# Conda environment:
#   source: workflow/envs/dada2.yaml
#   prefix: /conda-envs/96b8555f5e904d4efea1f53cab1174ba
#   # code/environments/exact/270c1d6d608fcc294d8576b354642db3_.yaml
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - r-matrix=1.3-2 # https://github.com/benjjneb/dada2/issues/212
#   - bioconductor-dada2=1.16.0
#   - r-tidyverse=1.3.2
#   - r-ggplot2=3.3.6
#   - bioconductor-rhdf5=2.34.0
#   - pigz=2.6
RUN mkdir -p /conda-envs/96b8555f5e904d4efea1f53cab1174ba
COPY workflow/envs/dada2.yaml /conda-envs/96b8555f5e904d4efea1f53cab1174ba/environment.yaml

# Conda environment:
#   source: workflow/envs/mmseqs2-vsearch.yaml
#   prefix: /conda-envs/0d7c8fda1d890a7c7d1a17363ffc35f0
#   # code/environments/exact/ebca16bd2afa552b3e747a9b50274acd_.yaml
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - anndata=0.8.0
#   - biom-format=2.1.14
#   - biopython=1.80
#   - mmseqs2=14.7e284
#   - pandas=1.5.2
#   - python=3.11.0
#   - samtools=1.16.1
#   - vsearch=2.22.1
#   - zlib=1.2.13
#   - pip
#   - pip:
#     - "git+https://github.com/caseygrun/nbseq.git#egg=nbseq"
RUN mkdir -p /conda-envs/0d7c8fda1d890a7c7d1a17363ffc35f0
COPY workflow/envs/mmseqs2-vsearch.yaml /conda-envs/0d7c8fda1d890a7c7d1a17363ffc35f0/environment.yaml

# Conda environment:
#   source: workflow/envs/msa.yaml
#   prefix: /conda-envs/f91e78c37d6079c10d247f1d2fb7464a
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - biopython=1.79
#   - pandas=1.4.3
#   - python=3.10.5
#   - joblib=1.1.0
#   - fasttree=2.1.11
#   - mafft=7.525
#   - trimal=1.4.1
#   - clustalo
#   - muscle
RUN mkdir -p /conda-envs/f91e78c37d6079c10d247f1d2fb7464a
COPY workflow/envs/msa.yaml /conda-envs/f91e78c37d6079c10d247f1d2fb7464a/environment.yaml

# Conda environment:
#   source: workflow/envs/scran-deseq2.yaml
#   prefix: /conda-envs/56d494dab6ae07c607c6629258d8c91c
#   # code/environments/exact/e63c021a43157a12bd145e585e5e1fb8.yaml
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - anndata=0.8.0
#     - bioconductor-ancombc=1.4.0
#     - bioconductor-biomformat=1.22.0
#     - bioconductor-deseq2=1.34.0
#     - bioconductor-metagenomeseq=1.36.0
#     - bioconductor-scater=1.22.0
#     - bioconductor-scran=1.22.1
#     - bioconductor-wrench=1.12.0
#     - biom-format=2.1.12
#     - numpy=1.23.1
#     - python=3.10.5
#     - r=4.1
#     - r-cowplot=1.1.1
#     - r-ggplot2=3.3.6
#     - r-ggrepel=0.9.1
#     - r-gridextra=2.3
#     - r-irkernel=1.3
#     - r-reticulate=1.25
#     - r-tidyverse=1.3.2
#     - scipy=1.8.1
#     - pip
#     - pip:
#       - "git+https://github.com/caseygrun/nbseq.git#egg=nbseq"
RUN mkdir -p /conda-envs/56d494dab6ae07c607c6629258d8c91c
COPY workflow/envs/scran-deseq2.yaml /conda-envs/56d494dab6ae07c607c6629258d8c91c/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/12118d3115c8565e50c5fda969d0a5bf --file /conda-envs/12118d3115c8565e50c5fda969d0a5bf/environment.yaml && \
    mamba env create --prefix /conda-envs/3a8e66a7ad949aaf6a4d9cb01a91dea1 --file /conda-envs/3a8e66a7ad949aaf6a4d9cb01a91dea1/environment.yaml && \
    mamba env create --prefix /conda-envs/6a5eb90c6ff9a883e31a798150db79a0 --file /conda-envs/6a5eb90c6ff9a883e31a798150db79a0/environment.yaml && \
    mamba env create --prefix /conda-envs/8448d2cbdccdc2abdf724bf1dcbe400f --file /conda-envs/8448d2cbdccdc2abdf724bf1dcbe400f/environment.yaml && \
    mamba env create --prefix /conda-envs/dd444f859e73c5e562ebcf217729c29c --file /conda-envs/dd444f859e73c5e562ebcf217729c29c/environment.yaml && \
    mamba env create --prefix /conda-envs/96b8555f5e904d4efea1f53cab1174ba --file /conda-envs/96b8555f5e904d4efea1f53cab1174ba/environment.yaml && \
    mamba env create --prefix /conda-envs/0d7c8fda1d890a7c7d1a17363ffc35f0 --file /conda-envs/0d7c8fda1d890a7c7d1a17363ffc35f0/environment.yaml && \
    mamba env create --prefix /conda-envs/f91e78c37d6079c10d247f1d2fb7464a --file /conda-envs/f91e78c37d6079c10d247f1d2fb7464a/environment.yaml && \
    mamba env create --prefix /conda-envs/56d494dab6ae07c607c6629258d8c91c --file /conda-envs/56d494dab6ae07c607c6629258d8c91c/environment.yaml && \
    mamba clean --all -y
