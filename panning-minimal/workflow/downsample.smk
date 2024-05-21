# include: 'rules/common.smk'

import os
import shutil
import pandas as pd

configfile: 'config/config.yaml'

src_seq_dir = '/home/cng2/palmer_scratch/input/2022-12-13-CG026n/'
unsampled_seq_dir = '/home/cng2/palmer_scratch/input/panning-minimal/'
runs = ['20221208']
run = runs[0]
# glob = src_seq_dir + '/{sample}/'
# print(glob)
# samples, = glob_wildcards(glob)

samples = [ f.name for f in os.scandir(src_seq_dir) if f.is_dir() ]
samples = pd.read_csv('config/samples.tsv', sep="\t").ID.values

print(f"Copying the following samples from '{src_seq_dir}' to '{unsampled_seq_dir}', ")
print(f"then downsampling to '{str(Path(config['raw_sequence_dir']) / run)}'")
print(samples)

rule all:
	input: expand(str(Path(config['raw_sequence_dir']) / '{run}/Sample_{sample}/'), run=runs, sample=samples)




# src_seq_dir       = '/home/cng2/palmer_scratch/input/2022-12-13-CG026n/'
# unsampled_seq_dir = '/home/cng2/palmer_scratch/input/panning-minimal/'
# pattern = 'Sample_{ID}'

# os.makedirs(unsampled_seq_dir, exist_ok=True)

# for ID in samples['ID']:
# 	sample_folder = pattern.format(ID=ID)

# 	src, dest = Path(src_seq_dir) / sample_folder, Path(dst_seq_dir) / sample_folder
# 	print(f"{src} -> {dest}")
# 	shutil.copytree(src, dest, dirs_exist_ok=True)

rule copy_seqs:
	input:  str(Path(src_seq_dir) / '{sample}/')
	output: directory(str(Path(unsampled_seq_dir) / '{sample}/'))
	run:
		for i,o in zip(input, output):
			print(f"{i} -> {o}")
			shutil.copytree(i, o, dirs_exist_ok=True)		

rule downsample:
	conda: 'etc/downsample.yaml'
	input: str(Path(unsampled_seq_dir) / '{sample}/')
	output: directory(str(Path(config['raw_sequence_dir']) / run / "{sample}/"))
	params:
		reads=7500
	shell:"""
	shopt -s nullglob
	mkdir -p {output}
	
	for f in "{input}"/*.fastq.gz; do
		seqtk sample "$f" {params.reads} | gzip > "{output}/$(basename $f)"
	done
	"""