snakemake@source('common.R')
snakemake_log_start(snakemake)
options(mc.cores = snakemake@threads)

library(dada2)
library(parallel)
library(Matrix)


ft <- read_feature_table_sparse(snakemake@input[["feature_table"]])
aligned_asvs <- read.csv(snakemake@input[["aligned_asvs"]],
	stringsAsFactors=FALSE, check.names=FALSE)

ft_dir <- snakemake@output[["feature_table"]]
ft_aligned <- ft[,aligned_asvs$ASVID]
message("Writing feature table for aligned ASVs to ", ft_dir)
write_feature_table_sparse(ft_aligned, dir = ft_dir)

summary_path <- snakemake@output[["summary"]]
message("Writing summary of aligned reads per sample to ", summary_path)
write_tsv(
	summarize_samples(ft_aligned, name = 'aligned'),
	summary_path)

snakemake_log_end()
