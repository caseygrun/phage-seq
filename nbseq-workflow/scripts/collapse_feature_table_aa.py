from common import *
import pandas as pd
import nbseq.ft

with snakemake_log(snakemake) as logfile:

	assert any(x in snakemake.input.keys() for x in ['features','mapping']), f"Error: 'features' and/or 'mapping' file(s) required. Input keys: {snakemake.input.keys()}"

	# mapping of NA ASVID to AA aaSVID
	mapping = None
	if 'mapping' in snakemake.input.keys():
		mapping = pd.concat([read_delim_auto(f) for f in listify(snakemake.input['mapping'])])
		mapping = mapping.drop_duplicates() #.drop_duplicates(['ASVID','aaSVID'])

		if 'mapping_cols' in snakemake.params.keys():
			mapping_cols = list(snakemake.params['mapping_cols'])
			print(f"Using custom mapping columns: {mapping_cols}")
		else:
			mapping_cols = list(mapping.columns)
			print(f"Using first two columns from CSV/TSV file: {mapping_cols}")
		print(f"Mapping feature table IDs from {mapping_cols[0]} to {mapping_cols[1]}:")
		print(mapping[mapping_cols])


	if 'axis' in snakemake.params.keys():
		axis = snakemake.params['axis']
	else: axis = 'observation'

	print(f"Using axis: {axis}")

	# features (aaSVs) on which to filter
	ids = None
	if 'features' in snakemake.input.keys():
		ids = pd.concat([pd.read_csv(f) for f in snakemake.input['features']])

		if 'filter_col' in snakemake.params.keys(): id_col = snakemake.params['filter_col']
		else: id_col = 'aaSVID'
		ids = ids.drop_duplicates(id_col)
		print(f"Filtering feature table to {len(ids)} {axis}s from column {id_col}...")

		if mapping is not None:
			mapping = mapping.loc[mapping.loc[:,mapping_cols[0]].isin(ids[id_col]),:]


	ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
	print(f"Input feature table: {repr(ft)}")

	if mapping is not None:

		# make sure any ASVs that aren't in the mapping are dropped
		ft = ft.filter(mapping[mapping_cols[0]], axis=axis)
		print(f"feature table after filtering: {repr(ft)}")

		# group feature table, mapping NA ASVID to AA aaSVID
		# mapping = mapping.set_index('ASVID')['aaSVID']
		# set index as the 0th column, then get a series with values in the 1st column
		mapping = mapping.set_index(mapping_cols[0]).loc[:,mapping_cols[1]]
		assert mapping.index.is_unique

		ft = nbseq.ft.group_feature_table(ft, mapping, axis=axis)
		print(f"feature table after grouping: {repr(ft)}")
	elif ids is not None:

		ft = nbseq.ft.filter_feature_table(ft, ids[id_col])
		print(f"feature table after filtering: {repr(ft)}")

	outfile = snakemake.output[0]
	print(f"Writing output to {outfile}")
	nbseq.ft.write_feature_table_biom(ft, outfile)
