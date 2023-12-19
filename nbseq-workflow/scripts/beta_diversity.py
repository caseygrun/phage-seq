import nbseq.ft
import nbseq.distance

from common import snakemake_log
with snakemake_log(snakemake) as logfile:

	ft = nbseq.ft.read_feature_table_biom(snakemake.input['feature_table'])
	distance = nbseq.distance.pairwise_distance(
		ft,
		metric=snakemake.wildcards['metric'], threads=snakemake.threads)

	nbseq.distance.save_distance_matrix(distance, snakemake.output[0])

