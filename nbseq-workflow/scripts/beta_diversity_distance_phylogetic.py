from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	import nbseq.ft
	import nbseq.distance
	import skbio.tree

	# ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
	phylogeny = skbio.tree.TreeNode.read(snakemake.input['phylogeny'])

	# distance = nbseq.distance.pairwise_distance(
	# 	ft,
	# 	metric=snakemake.wildcards['metric'],
	# 	tree=phylogeny
	# )
	distance = nbseq.distance.pairwise_distance_phylo(
		snakemake.input['feature_table'],
		tree_path = snakemake.input['phylogeny'],
		metric=snakemake.wildcards['metric'],
		outfile=snakemake.output[0],
		threads=snakemake.threads,
		intersection=True)

	# nbseq.distance.save_distance_matrix(distance, snakemake.output[0])
