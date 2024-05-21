import shutil
from pathlib import Path
import os
import pandas as pd
import yaml

samples = pd.read_csv('samples.tsv', sep="\t")


src_seq_dir = '/home/cng2/palmer_scratch/input/2022-12-13-CG026n/'
dst_seq_dir = '/home/cng2/palmer_scratch/input/panning-minimal/'
pattern = 'Sample_{ID}'

os.makedirs(dst_seq_dir, exist_ok=True)

for ID in samples['ID']:
	sample_folder = pattern.format(ID=ID)

	src, dest = Path(src_seq_dir) / sample_folder, Path(dst_seq_dir) / sample_folder
	print(f"{src} -> {dest}")
	shutil.copytree(src, dest, dirs_exist_ok=True)



# with open("config.yaml") as stream:
#     try:
#         config = yaml.safe_load(stream)
#         dst_seq_dir = config['raw_sequence_dir']

#         for ID in samples['ID']:
#         	sample_folder = pattern.format(ID=ID)

#         src, dest = Path(src_seq_dir) / sample_folder, Path(dst_seq_dir) / sample_folder
#         print(f"{src} -> {dest}")
#         shutil.copy(src, dest)

#     except yaml.YAMLError as exc:
#         print(exc)