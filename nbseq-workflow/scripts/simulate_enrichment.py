import pickle

import numpy as np

import nbseq.utils
import nbseq.ft
import nbseq.select

from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	metadata = nbseq.utils.read_delim_auto(snakemake.input['metadata']).set_index(snakemake.params['id_col'] if 'id_col' in snakemake.params.keys() else 'ID')
	ft = nbseq.ft.read(snakemake.input['feature_table'], obs=metadata)
	print('Input feature:')
	print(ft)
	print()

	# "expt == '027j.lib' & name != '027j.lib.P2.8'"
	query = snakemake.params['query']
	ft_lib = nbseq.ft.query(ft, query, axis='obs')
	print('Feature table used for simulating enrichments, after "{query}"')
	print(ft_lib)
	print()

	print("Simulating enrichments...")
	df_enr_sim = nbseq.select.simulate_enrichments(ft_lib)
	df_enr_sim['log_enrichment'] = np.log10(df_enr_sim.enrichment).replace([float('inf'), float('-inf')],0)
	print(f"Simulated {len(df_enr_sim)} enrichment observations")

	print("Fitting ECDF...")
	cond_ecdf = nbseq.select.enrichment_abundance_ecdf(df_enr_sim, plot=False, verbose=True)

	cond_ecdf_path = snakemake.output[0]
	print(f"Writing enrichment model to '{cond_ecdf_path}'")
	nbseq.utils.mkdirp_file(cond_ecdf_path)
	with open(cond_ecdf_path, 'wb') as f:
	    pickle.dump(cond_ecdf, f)