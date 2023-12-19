snakemake@source('common.R')
library(dada2)

# get list of summaries (one per sample), then stack them into one data.frame
filtered.summary <- data.frame(
	do.call(rbind,
		lapply(snakemake@input[["filter_summaries"]], readRDS)
	)
)
colnames(filtered.summary) <- c('prefilter','filtered')
filtered.summary <- tibble::rownames_to_column(filtered.summary, var="sample")
write_tsv(filtered.summary, snakemake@output[["summary_filter"]])
