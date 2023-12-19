from common import *
import pandas as pd

with snakemake_log(snakemake) as logfile:

	# we need a unique mapping of member ASVIDs to cluster ASVIDs, otherwise
	# group_feature_table breaks. This script handles the case where two libraries
	# assigned the same ASVID to different clusters (usually one of them to a
	# cluster and the other to a singleton cluster).

	# these are uncommon edge cases and probably just represent contaminants or
	# sequencing errors anyway

	infiles = list(snakemake.input)
	print(f"Resolving clusters from {infiles}...")
	mapping = pd.concat([read_delim_auto(f) for f in infiles])

	# same member ID, same cluster ID
	mapping = mapping.drop_duplicates() #.drop_duplicates(['ASVID','aaSVID'])

	if 'mapping_cols' in snakemake.params.keys():
		mapping_cols = list(snakemake.params['mapping_cols'])
		print(f"Using custom mapping columns: {mapping_cols}")
	else:
		mapping_cols = list(mapping.columns)
		print(f"Using first two columns from CSV/TSV file: {mapping_cols}")
	print(f"Mapping feature table IDs from {mapping_cols[0]} to {mapping_cols[1]}:")
	print(mapping)

	# some ASVIDs may have been clustered in one library but not clustered in the
	# other. Since the clustering is extremely conserative, these clusters are
	# probably real. We will therefor prefer the assignment where
	# member_ASVID != cluster_ASVID , since that implies the ASVID is not in a
	# cluster by itself. In the case this is true for both, the tie is broken
	# arbitrarily.
	mapping["_member_equals_cluster"] = (mapping[mapping_cols[0]] == mapping[mapping_cols[1]])
	mapping = mapping.sort_values("_member_equals_cluster").groupby(mapping_cols[0]).first()

	# make sure we got rid of any duplicates
	assert mapping.index.is_unique
	mapping.to_csv(snakemake.output[0], index=True)
