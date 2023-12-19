snakemake@source('common.R')
snakemake_log_start(snakemake)

library(dada2)
library(parallel)
library(Matrix)

snakemake_use_threads(snakemake)

# MAKE FEATURE TABLE
# ==================

tmessage("Loading merged pairs...")
mergers <- readRDS(snakemake@input[["mergers"]])

tmessage("Size: mergers:              ", format(object.size(mergers),units='auto'))
tmessage("Merged pairs from ", length(mergers), " samples")
tmessage("Filtering merged sequences...")

# # now filter out the rejects before making a sequence table
# # we gave a list of inputs, so mergers is a list of `data.frame`s, each of which
# # has a column $accept = TRUE for merged sequences DADA2 wants to keep
# mergers <- lapply(mergers, function(sample_mergers) {
#     return(sample_mergers[sample_mergers$accept, , drop = FALSE])
# })

tmessage("Identifying unique sequences for each sample...")


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


# get a list of unique ASVs per sample, as well as a list of all unique ASVs
# and their MD5 hashes
uniq_seqs_per_sample <- mclapply_retry(seq_along(mergers), function(i) getUniques(mergers[[i]]) ) #mclapply(mergers, getUniques)
names(uniq_seqs_per_sample) <- names(mergers)

# out_classes <- lapply(uniq_seqs_per_sample,class)
# if (any(out_classes == "try-error")) {
#     print(as.character(uniq_seqs_per_sample[out_classes == "try-error"]))
#     stop()
# }

rm(mergers)
tmessage("Size: uniq_seqs_per_sample: ", format(object.size(uniq_seqs_per_sample),units='auto'))
gc()

tmessage("Merging unique sequences across all samples...")
all_uniq_seqs <- unique(do.call(c, mclapply(uniq_seqs_per_sample, names)))

tmessage(length(all_uniq_seqs), " unique sequences identified.")
tmessage("Size: all_uniq_seqs:        ", format(object.size(all_uniq_seqs),units='auto'))

tmessage("Hashing unique sequences...")
all_uniq_seqs_to_hash <- unlist(
    mclapply(all_uniq_seqs, digest::digest, algo="md5", serialize=FALSE))
names(all_uniq_seqs_to_hash) <- all_uniq_seqs



makeSequenceTableFast <- function(samples, all_uniq_seqs, verbose=TRUE, parallel=TRUE) {
    gc()
    start_time <- proc.time()

    # Samples are rows, columns are sequences
    # dense matrix of samples x ASVs will be large:
    if (verbose) cat("size feature table as dense matrix: ", format_size(4*length(samples)*length(all_uniq_seqs)), "\n")

    # with a sparse matrix, we only need an entry for each non-zero element
    n_entries <- sum(unlist(mclapply(samples, length)))
    if (verbose) cat("expected size of feature table as sparse matrix (approx): ", format_size(4 * 3 * n_entries), "\n")

    if (parallel) {
        if (verbose) cat("using parallel processing... ")
        rval <- makeSequenceTableSparseParallel(samples, all_uniq_seqs, verbose)
    } else {
        rval <- makeSequenceTableSparse(samples, all_uniq_seqs, n_entries = n_entries, verbose = verbose)
    }
    if (verbose) cat("\n","actual size of feature table as sparse matrix: ", format(object.size(rval),units='auto'),"\n", sep="")
    dimnames(rval) <- list(
        names(samples),
        all_uniq_seqs
    )

    if (verbose) cat("\n","final size of feature table with annotations: ", format(object.size(rval),units='auto'),"\n", sep="")
    end_time <- proc.time()
    if (verbose) {
        cat("finished:\n")
        print(end_time - start_time)
    }

    return(rval)
}

makeSequenceTableSparseParallel <- function(samples, all_uniq_seqs, verbose=TRUE) {
    tmessage("Identifying index pairs per sample...")
    indexes = seq_along(samples)
    get_triples_per_sample <- function(i) {
        j_s <- match(names(samples[[i]]), all_uniq_seqs)
        sample_length <- length(j_s)
        i_s <- rep(i, sample_length)
        x_s <- samples[[i]]
        return(list(i = i_s, j = j_s, x = x_s))
    }
    triples_per_sample <- mclapply_retry(indexes, get_triples_per_sample, expect_class="list")

    tmessage("Merging pairs per sample into sparse matrix...")
    triples <- do.call(mapply, c(c, triples_per_sample, SIMPLIFY=FALSE))
    rval <- sparseMatrix(i = triples[['i']],
                         j = triples[['j']],
                         x = triples[['x']],
                         index1 = TRUE,
                         dims=c(length(samples),
                                length(all_uniq_seqs))
                         )
    return(rval)
}

# obsolete alternative implementation
makeSequenceTableSparse <- function(samples, all_uniq_seqs, n_entries, verbose=TRUE) {

    # establish sparse matrix triples
    i_s <- numeric(n_entries) # index of sample
    j_s <- numeric(n_entries) # index of ASV (in all_uniq_seqs)
    x_s <- numeric(n_entries) # abundance of ASV in sample i,j

    # index of next triple
    entry <- 1

    # for each sample
    for(i in seq_along(samples)) {

        # get ASVs in sample as indexes in all_uniq_seqs
        sample_js <- match(names(samples[[i]]), all_uniq_seqs)

        # figure out which triples the sample will occupy
        sample_entry_start <- entry
        sample_length <- length(sample_js)
        sample_entry_end <- entry + sample_length - 1

        # populate triples
        i_s[sample_entry_start:sample_entry_end] <- i
        j_s[sample_entry_start:sample_entry_end] <- sample_js
        x_s[sample_entry_start:sample_entry_end] <- samples[[i]]

        # incrememnt pointer
        entry <- entry + sample_length

        if (verbose) cat(names(samples)[[i]],", ")
    }

    rval <- sparseMatrix(i = i_s, j = j_s, x = x_s,
                         index1 = TRUE,
                         dims = c(length(samples),
                                  length(all_uniq_seqs))
                         )

    return(rval)
}

tmessage("Making feature table...")
ft <- makeSequenceTableFast(uniq_seqs_per_sample, all_uniq_seqs, parallel=TRUE)

output_folder <- snakemake@output[["feature_table"]]
tmessage("Writing feature table to ", output_folder)

# write feature table with colnames being hashes (ASVIDs)
write_feature_table_sparse(ft, dir = output_folder,
    features = all_uniq_seqs_to_hash,
    samples = names(uniq_seqs_per_sample))

snakemake_log_end()
