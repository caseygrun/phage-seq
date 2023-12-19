from common import *
import pandas as pd
import nbseq.ft
import numpy as np

metadata = read_delim_auto(snakemake.input['metadata'])

param = snakemake.params['param']
axis = snakemake.params['axis'] if 'axis' in snakemake.params.keys() else 'sample'
metadata_id_col = snakemake.params['id_col'] if 'id_col' in snakemake.params.keys() else 'ID'

print(f"Filtering feature table by metadata column {param} == {snakemake.wildcards[param]}...")
metadata_filtered = (metadata.loc[
	metadata[param] == snakemake.wildcards[param], :
])
print(f"- {len(metadata_filtered)} {axis}s after filtering...")

ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
print(f"- Feature table before filtering: {repr(ft)}")

# ids_to_keep = np.intersect1d(metadata_filtered.iloc[:,0].values,
# 	ft.ids(axis),
# 	assume_unique = True)
ft, ids_to_keep = nbseq.ft.filter_intersect(ft,
	metadata_filtered.loc[:,metadata_id_col].values,
	axis=axis, remove_empty=True)

print(f"- {len(ids_to_keep)} {axis}s after intersection with feature table...")


if 'metadata' in snakemake.output.keys():
	(metadata_filtered
		.set_index(metadata_id_col)
		.loc[ids_to_keep,:]
		.to_csv(snakemake.output['metadata'], index=True))

	del metadata_filtered

# keep only feature IDs from the correct library
ft.filter(ids_to_keep, axis=axis, inplace=True)

# drop {other_axis}s that have no {axis}s after filtering
nbseq.ft.drop_empty(ft,
	axis={'sample':'observation', 'observation':'sample'}[axis],
	inplace=True)
print(f"- Feature table after filtering: {repr(ft)}")
nbseq.ft.write_feature_table_biom(ft, snakemake.output['feature_table'])


if ('feature_data' in snakemake.output.keys()):
	feature_data = read_delim_auto(snakemake.input['feature_data'])
	ft, ids_to_keep = nbseq.ft.filter_intersect(ft,
		feature_data.iloc[:,0].values,
		axis='observation', remove_empty=True)
		
	# indexing into pd.DataFrame with set is deprecated
	ids_to_keep = list(ids_to_keep)

	(feature_data
		.set_index(list(feature_data.columns)[0])
		.loc[ids_to_keep,:]
		.to_csv(snakemake.output['feature_data'], index=True))
