from pathlib import Path
import tempfile
import subprocess

import pandas as pd
import nbseq.utils
from nbseq.asvs.cluster import cluster_asvs
from common import *

# with open(snakemake.log[0],'w') as logfile:
with snakemake_log(snakemake) as logfile:

	# contains one row per read pair:
	# ASVID,seq,pair_id,fwd,int,rev
	# (Must have ASVID, seq, pair_id; other columns optional)
	# ASVID is the MD5 hash of the full length ungapped `seq`
	input = pd.read_csv(snakemake.input[0])

	output = cluster_asvs(input, cluster_path=snakemake.output['clusters'], threads=snakemake.threads)
	output.to_csv(snakemake.output['na'], index=False)
