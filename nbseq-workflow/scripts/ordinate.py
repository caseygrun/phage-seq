import pickle
import nbseq.ft
from nbseq.ordination import *

from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	ft = nbseq.ft.read_feature_table_biom(snakemake.input[0])
	method = snakemake.wildcards['method']
	params = snakemake.params

	if method == 'TSVD':
		kwargs = { 'n_components':100, **params }
	elif method == 'TSVD-TSNE':
		kwargs = { 'decomp_kwargs':{'n_components':100}, **params }
	else:
		kwargs = params

	print(f"Ordinating with method `{method}`: {kwargs}")
	ord_skl, ord_skbio = ordinate(ft, method, **kwargs)

	print(f"Writing skbio output to {snakemake.output['skbio']}")
	ord_skbio.write(snakemake.output['skbio'])

	print(f"Writing pickled sklearn output to {snakemake.output['sklearn']}")
	with open(snakemake.output['sklearn'],'wb') as f:
		pickle.dump(ord_skl, f)
