import pickle
import nbseq.ft
from nbseq.ordination import *

from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	ft = nbseq.ft.read_feature_table_biom(snakemake.input[0])
	method = snakemake.wildcards['method']
	params = snakemake.params

	# n_components = min([n_samples, n_features, 100])
	n_components = min(list(ft.shape) + [100])
	if method == 'TSVD':
		kwargs = { 'n_components':n_components, **params }
	elif method == 'TSVD-TSNE':
		perplexity = min([float(ft.shape[1])/2, 30.0])
		kwargs = { 'decomp_kwargs':{'n_components':n_components}, 'perplexity':perplexity, **params }
	else:
		kwargs = params

	print(f"Ordinating with method `{method}`: {kwargs}")
	ord_skl, ord_skbio = ordinate(ft, method, **kwargs)

	print(f"Writing skbio output to {snakemake.output['skbio']}")
	ord_skbio.write(snakemake.output['skbio'])

	print(f"Writing pickled sklearn output to {snakemake.output['sklearn']}")
	with open(snakemake.output['sklearn'],'wb') as f:
		pickle.dump(ord_skl, f)
