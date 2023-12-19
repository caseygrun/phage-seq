import numpy as np
from common import *

with snakemake_log(snakemake):
	def norm_scran(input, output):
		from nbseq.norm import scran
		scran(ft_path=input, 
			out_path=output,
			verbose=True)


	norm_functions = {
		'scran': norm_scran
	}

	norm_f = norm_functions[snakemake.wildcards.norm]
	norm_f(snakemake.input[0], snakemake.output[0])