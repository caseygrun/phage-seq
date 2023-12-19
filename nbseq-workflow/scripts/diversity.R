library(tidyverse)
library(magrittr)
library(Matrix)
library(breakaway)
library(parallel)

snakemake@source('common.R')
snakemake_log_start(snakemake)
snakemake_use_threads(snakemake)
options(error = rlang::entrace)

message("Reading feature table...")
ft <- read_feature_table_hd5(snakemake@input[['feature_table']]) #read_feature_table_sparse(snakemake@input[['feature_table']])


mclapply_retry <- function(data,f, n_retries=3, expect_class=NULL) {
	library(parallel)

    # first try
    out <- mclapply(data, f)

    # examine results
    out_classes <- lapply(out,class)
    if ( any(out_classes == "try-error") || (!is.null(expect_class) && any(out_classes != expect_class)) ) {
        if ( any(out_classes == "try-error") ) {
            tmessage("BAD: Errors in one or more cores: ", sum(out_classes == "try-error"), " items affected:")
            print(as.character(out[out_classes == "try-error"]))
        }
        if (!is.null(expect_class) && any(out_classes != expect_class)) {
            unexpected_classes <- (out_classes != expect_class)
            tmessage("BAD: Unexpected results in ", sum(unexpected_classes), " items:")
            print(head(as.character(out[unexpected_classes])))
        }

        # second try, only on failed results
        tmessage("Retrying...")
        bad_classes <- (out_classes == "try-error")
        retry <- mclapply(data[bad_classes], f)

        out_classes <- lapply(retry, class)
        if ( any(out_classes == "try-error") || (!is.null(expect_class) && any(out_classes != expect_class)) ) {
            if ( any(out_classes == "try-error") ) {
                tmessage("BAD: Errors in one or more cores: ", sum(out_classes == "try-error"), " items affected:")
                print(as.character(retry[out_classes == "try-error"]))
            }
            if (!is.null(expect_class) && any(out_classes != expect_class)) {
                unexpected_classes <- (out_classes != expect_class)
                tmessage("BAD: Unexpected results in ", sum(unexpected_classes), " items:")
                print(head(as.character(retry[unexpected_classes])))
            }
            tmessage("BAD: Unable to obtain correct results even with retry. Exiting...")
            stop()
        }
        # if successful, overwrite the bad cases
        out[bad_classes] <- retry
    }

    tmessage("OK: completed")
    return(out)
}

# ft should have samples in rows x observations in columns
get_alpha_diversity_breakaway <- function(ft) {
	ftt <- t(ft)
	n_samples <- ncol(ftt)
    calculate_sample_alpha_diversity <- function(i) {
        withCallingHandlers({
			message("- tabulating sample ", i," of ", n_samples)
			tab <- tabulate(ftt[,i])
			freq_table <- data.frame(
				index = seq_along(tab),
				freq = tab
			)
            message("- breakaway  sample ", i," of ", n_samples)
			x <- breakaway(freq_table,
						   answers=TRUE,
						   print = FALSE, # print = TRUE,
						   plot = FALSE)
			unlist(c(
				x[c('est','seest')],
				uci=x[['ci']][[2]],
				lci=x[['ci']][[1]]
			))
		},
		error = function(cnd) {
			print(cnd)
			print(lobstr::cst())
		})
	}

	alphas <- mclapply_retry(seq_len(n_samples),calculate_sample_alpha_diversity)
	names(alphas) <- colnames(ftt)

	# some cores may have returned errors, which R expresses bizarrely; filter
	# those out
	good_alphas <- purrr::keep(alphas, ~{
		if (inherits(alpha, 'try-error')) {
			message("Error in ", alpha)
		}
		!inherits(alpha, 'try-error')
	})

	data.frame(
		do.call(rbind, good_alphas)
	)

}

# ft should have observations in rows x samples in columns
get_diversity <- function(ft, metadata = NULL) {
    message("Calculating alpha diversity on ", nrow(ft), " samples...")
	# ft <- t(ft)
	richness <- rowSums(ft != 0)
	sample_reads <- rowSums(ft)
	reads_per_asv <- sample_reads / richness
    message("Running breakaway...")
	alpha <- get_alpha_diversity_breakaway(ft) %>%
		magrittr::set_colnames(paste0('richness_',colnames(.)))


	diversity <- data.frame(
		sample_richness = richness,
		reads_per_asv = reads_per_asv,
		asvs_per_read = 1/reads_per_asv,
		sample_reads = sample_reads,
		row.names = rownames(ft)
	)

	# print(diversity)
	# print(alpha)

	# diversity <- (cbind(diversity, alpha[rownames(diversity),]) %>% rownames_to_column(var = 'ID'))

	# some breakaway models may not converge
	diversity <- (merge(
		diversity, alpha,
		by='row.names', all=TRUE
	) %>% rename('ID' = 'Row.names'))

	if (!is.null(metadata)) {
		return(inner_join(diversity, metadata, by = 'ID'))
	} else {
		return(diversity)
	}
}

# ft <- ft[1:100,]
# ft <- ft[800:980,]
system.time(
	diversity <- get_diversity(ft)
)

output_file <- snakemake@output[['alpha_diversity']]
message("Writing output to ", output_file, " ...")
readr::write_tsv(diversity, file=output_file)
message("Done.")

snakemake_log_end()
