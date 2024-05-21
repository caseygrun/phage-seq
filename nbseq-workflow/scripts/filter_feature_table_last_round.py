from common import *
import pandas as pd
import nbseq.ft
import numpy as np

metadata = read_delim_auto(snakemake.input['metadata'])

param = snakemake.params['param'] if 'axis' in snakemake.params.keys() else 'round'
axis = 'sample'
metadata_id_col = snakemake.params['id_col'] if 'id_col' in snakemake.params.keys() else 'ID'


# feature table may already contain a subset of samples (e.g. one expt); take intersection of
# feature table and metadata to start
ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
print(f"- Feature table before filtering: {repr(ft)}")

ft, ids_to_keep = nbseq.ft.filter_intersect(ft,
	metadata.loc[:,metadata_id_col].values,
	axis=axis, remove_empty=True)

print(f"- {len(ids_to_keep)} {axis}s after intersection between metadata and feature table...")
metadata_filtered = metadata.loc[:,metadata[metadata_id_col].isin(ids_to_keep)]


# find the last round
rounds = sorted(metadata[param].unique())
print(f"Rounds: {rounds}")
print(f"Filtering feature table to {param} = '{rounds[-1]}'")

# filter metadata and intersect with feature table again
metadata_filtered = (metadata_filtered.loc[
	metadata_filtered[param] == rounds[-1], :
])
print(f"- {len(metadata_filtered)} samples after filtering...")
ft, ids_to_keep = nbseq.ft.filter_intersect(ft,
	metadata_filtered.loc[:,metadata_id_col].values,
	axis=axis, remove_empty=True)


# if indicated, write filtered metadata then delete to save memory
if 'metadata' in snakemake.output.keys():
	(metadata_filtered
		.set_index(metadata_id_col)
		.loc[ids_to_keep,:]
		.to_csv(snakemake.output['metadata'], index=True))

	del metadata_filtered

# keep only feature IDs from the correct round
ft.filter(ids_to_keep, axis=axis, inplace=True)

# drop {other_axis}s that have no {axis}s after filtering
nbseq.ft.drop_empty(ft,
	axis={'sample':'observation', 'observation':'sample'}[axis],
	inplace=True)
print(f"- Feature table after filtering: {repr(ft)}")
nbseq.ft.write_feature_table_biom(ft, snakemake.output['feature_table'])

# if indicated, write filtered feature data (e.g. asvs.csv)
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
