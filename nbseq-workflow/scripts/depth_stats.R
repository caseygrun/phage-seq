snakemake@source('common.R')
snakemake_log_start(snakemake)

library(parallel)
library(tidyverse)
library(Matrix)
# library(rhdf5)


snakemake_use_threads(snakemake)


sample_to_freq_table <- function(ft, sample, colnames=c('reads_per_feature','n_features','ID')) {
	tab <- tabulate(ft[sample,])
	df <- data.frame(
	    reads_per_feature = seq_along(tab),
	    n_features = tab,
        'ID' = sample
	)
    if (!is.null(colnames)) {
        names(df) <- colnames
    }
    df
}


feature_stats <- function(freq_table) {
    freq_table %>%
        arrange(desc(reads_per_feature)) %>%
        filter(n_features > 0) %>%
        mutate(
            cum_n_features = cumsum(n_features),
            frac_n_features = cum_n_features/sum(n_features),
            cum_total_reads = cumsum(n_features * reads_per_feature),
            frac_total_reads = cum_total_reads / max(cum_total_reads))
}


calculate_sample_feature_stats <- function(ft, samples=NULL) {
    if (is.null(samples)) {
        samples <- rownames(ft)
    }
    mclapply(samples, function(s) {
        sample_to_freq_table(ft, s) %>% feature_stats
    }) %>% {do.call(rbind, .)}
}


ft <- read_feature_table_hd5(snakemake@input[[1]])

saveRDS(
	calculate_sample_feature_stats(ft),
	snakemake@output[[1]]
)

snakemake_log_end()