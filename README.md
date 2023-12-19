Code used to process the high-throughput sequencing data in our manuscript "Bacterial cell surface characterization by phage display coupled to high-throughput sequencing." 


This repository contains four components:

- `nbseq-workflow` is a template Snakemake-based workflow for processing raw Phage-seq sequencing data. The other directories include symbolic links to the `scripts`, conda environments (`envs`), and `rules` defined in this directory. However, each specific experiment contains its own configuration, resources, and Snakemake workflow definitions (`Snakefile` and `*.smk`).
- `panning-small` contains code used to process the data from solid-phase and small-scale cell-based  phage display selection campaigns reported in Fig. 2
- `panning-massive` contains code used to process the data from high-throughput cell-based phage display selection campaigns reported in Fig. 3D
- `panning-massive` contains code used to process the data from the extended rounds of high-throughput cell-based phage display selection reported in Fig. 3E

