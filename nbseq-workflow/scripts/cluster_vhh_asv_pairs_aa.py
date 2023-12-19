from pathlib import Path
import tempfile
import subprocess

import pandas as pd
import nbseq.utils
from nbseq.asvs.cluster import cluster_asvs

from common import *

# with open(snakemake.log[0],'w') as logfile:
with snakemake_log(snakemake) as logfile:

	# contains one row per aaSV
	# aaSVID,seq
	input = pd.read_csv(snakemake.input[0])

	output = cluster_asvs(input,
		asvid_col='aaSVID',
		seq_col='translated',
		read_id_col=None,
		cluster_path=snakemake.output['clusters'],
		threads=snakemake.threads)
	output.to_csv(snakemake.output['aa'], index=False)
