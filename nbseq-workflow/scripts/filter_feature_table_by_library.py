from common import *
import pandas as pd
import nbseq.ft
import numpy as np

feature_data = read_delim_auto(snakemake.input['feature_data'])

print(f"Filtering feature table to library {snakemake.wildcards['library']}...")
feature_data_filtered = (feature_data.loc[
	feature_data['library'] == snakemake.wildcards['library'], :
])
print(f"- {len(feature_data_filtered)} features after filtering...")

ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
print(f"- Feature table before filtering: {repr(ft)}")

# ids_to_keep = np.intersect1d(feature_data_filtered.iloc[:,0].values,
# 	ft.ids('observation'),
# 	assume_unique = True)
ft, ids_to_keep = nbseq.ft.filter_intersect(ft,
	feature_data_filtered.iloc[:,0].values,
	axis='observation', remove_empty=True)

print(f"- {len(ids_to_keep)} features after intersection with feature table...")

(feature_data_filtered
	.set_index(list(feature_data_filtered.columns)[0])
	.loc[list(ids_to_keep),:]
	.to_csv(snakemake.output['feature_data'], index=True))

del feature_data_filtered

# keep only feature IDs from the correct library
ft.filter(ids_to_keep, axis='observation', inplace=True)

# drop samples that have no features after filtering
nbseq.ft.drop_empty(ft,axis='sample', inplace=True)
print(f"- Feature table after filtering: {repr(ft)}")
nbseq.ft.write_feature_table_biom(ft, snakemake.output['feature_table'])
