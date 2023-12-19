snakemake_log_start <- function(snakemake) {
	if(length(snakemake@log) == 0) { return() }
	log.file <- file(snakemake@log[[1]],open="wt")
	sink(log.file)
	sink(log.file,type="message")
	return(log.file)
}

snakemake_log_end <- function() {
	sink(type="message")
	sink()
}

snakemake_use_threads <- function(snakemake) {
	options('mc.cores' = snakemake@threads)
	message("Using parallel processing with ", getOption('mc.cores'), " cores...")
}

mkdirp <- function(...) {
	fn <- file.path(...)
	dir.create(fn, recursive = TRUE, showWarnings = FALSE)
	return(fn)
}

tmessage <- function(...) {
	return(message(format(Sys.time()), " ", ...))
}

write_feature_table_sparse <- function(ft, dir, samples = NULL, features = NULL) {
    dir.create(dir, recursive=TRUE, showWarnings = FALSE)
    writeMM(ft, file = file.path(dir, 'table.mtx'))
    if (is.null(samples)) samples <- rownames(ft)
    write.table(samples, file.path(dir, 'samples.tsv'), sep="\t", row.names = FALSE, col.names = FALSE)

    if (is.null(features)) features <- colnames(ft)
    write.table(features, file.path(dir, 'features.tsv'), sep="\t", row.names = FALSE, col.names = FALSE)
}

read_feature_table_sparse <- function(dir) {
    ft <- readMM(file = file.path(dir, 'table.mtx'))
    samples <- read.table(file.path(dir, 'samples.tsv'),
                          col.names = c('sample'),
                          check.names = FALSE, stringsAsFactors = FALSE)
    features <- read.table(file.path(dir, 'features.tsv'),
                          col.names = c('feature'),
                          check.names = FALSE, stringsAsFactors = FALSE)
    dimnames(ft) <- list(
        samples$sample,
        features$feature
    )
    return(ft)
}

read_biom_hd5_matrix <- function(biom_file) {
	x = rhdf5::h5read(biom_file,"/",read.attributes = TRUE)
	library(Matrix)

	generate_matrix <- function(x){

	    nrow = length(x$observation$ids)
	    ncol = length(x$sample$ids)


	    # x$sample$matrix stores data in compressed sparse column (CSC) format. csc_matrix in scipy, dgCMatrix in R Matrix
	    # x$observation$matrix stores data in compressed sparse row (CSR) format. csr_matrix in scipy, dgRMatrix in R Matrix

	    # for details, see
	    # class? CsparseMatrix
	    # class? dgCMatrix
	    dx = new("dgCMatrix",
	             # sample/matrix/indices      : <int32> A (nnz,) dataset containing the row indices (e.g., maps into observation/ids)

	             # i: Object of class "integer" of length nnzero (number of non-zero elements).
	             # These are the 0-based row numbers for each non-zero element in the matrix,
	             # i.e., i must be in 0:(nrow(.)-1).
	             i=as.integer(x$sample$matrix$indices),

	             # sample/matrix/indptr       : <int32> A (N+1,) dataset containing the compressed column offsets

	             # p: integer vector for providing pointers, one for each column, to the initial
	             # (zero-based) index of elements in the column. .@p is of length ncol(.) + 1,
	             # with p[1] == 0 and p[length(p)] == nnzero, such that in fact, diff(.@p) are
	             # the number of non-zero elements for each column. In other words, m@p[1:ncol(m)]
	             # contains the indices of those elements in m@x that are the first elements in the
	             # respective column of m.
	             p=as.integer(x$sample$matrix$indptr),

	             # sample/matrix/data         : <float64> A (nnz,) dataset containing the actual matrix data
	             # x: Object of class "numeric" - the non-zero elements of the matrix.
	             x=as.numeric(x$sample$matrix$data),

	             Dim=c(nrow, ncol),
	             Dimnames = list(x$observation$ids, x$sample$ids)

	            )
	    dx
	}
	generate_matrix(x)
}

read_feature_table_hd5 <- function(biom_file) {
	return(t(read_biom_hd5_matrix(biom_file)))
}

summarize_samples <- function(ft, name = 'sum') {
    sample_sums <- colSums(ft)
    res <- data.frame(
		'sample' = names(sample_sums),
        'sum' = sample_sums)
    colnames(res)[2] <- name
    return(res)
}

sample_to_freq_table <- function(fti) {
	tab <- tabulate(fti)
	data.frame(
	    reads = seq_along(tab),
	    freq = tab
	)
}

write_tsv <- readr::write_tsv

# write_tsv <- function(df, file, sep="\t", quote=FALSE, row.names=FALSE, ...) {
# 	return(write.table(df, file, sep=sep, quote=quote, row.names=row.names, ...))
# }

read_tsv <- readr::read_tsv

# read_tsv <- function(file, sep="\t", row.names=FALSE, ...) {
# 	return(read.table(file=file, sep=sep, row.names=row.names, ...))
# }

natsorted_factor <- function(...) {
    return(factor(..., levels=stringr::str_sort(unique(...),numeric=TRUE)))
}

figsize <- function(x, y) {
    if (length(x) == 2) {
        y = x[2]
        x = x[1]
    }
    options(repr.plot.width=x, repr.plot.height=y)
}

plot_quiet <- function(gg) {
    return(suppressMessages(suppressWarnings(print(gg))))
}


library(R6)
library(progress)

# override progress_bar to prevent it from inserting extra \r characters
progress_bar_file <- R6Class("progress_bar_file",inherit=progress_bar,
	private = list(
		clear_line = function(...) {},
		cursor_to_start = function(...) {}
	)
)


augment_metadata <- function(metadata, guids=NULL, sample_metadata=NULL) {
    metadata %<>% mutate(r = gsub('[Rio]','',round),
                     io = stringr::str_sub(round,-1),
                     # expt = ifelse(grepl('_',sample),paste0(expt,'_pooled'),expt),
                     sample = gsub("\\+",'',sample),
                     kind = ifelse(grepl('-',sample),'-','+'),
                     selection = gsub("[\\+\\-]",'',sample)
                    ) %>%
                group_by(expt, round, sample) %>%
                mutate(replicate=row_number(ID)) %>% ungroup() %>%
				mutate(name_full = paste(expt,sample,replicate,round,sep='.'),
					name = paste(expt,sample,replicate,sep='.')
						)

    if (!is.null(guids)) {
        metadata %<>% full_join(guids, by = 'ID')
    }
    if (!is.null(sample_metadata)) {
		# I am trying to standardize "selection" as the name for the particular
		# selection condition, e.g. row of the sample_metadata
		if ('selection' %in% colnames(sample_metadata)) {
			.by = c('expt',
				# `metadata` may have sample = '1' or sample = '1-'; attach
				# same sample metadata regardless
				'selection'='selection')
		} else {
			.by = c('expt',
				'selection'='sample')
		}
		metadata %<>% left_join(sample_metadata, by=.by) %>%
			mutate(description = paste(
				ifelse(is.na(genotype_CS), background_CS, paste(background_CS, genotype_CS)),
				"/",
				ifelse(is.na(background_S),
					ifelse(category == "clinical isolate", "clinical isolate","?"),
					ifelse(is.na(genotype_S),
						background_S,
						paste(background_S, genotype_S)))
				)
			) %>%
			separate(genotype_pair,sep="\\s*/\\s*",remove=FALSE, into=c('gene_CS','gene_S'))
    }
    return(metadata)
}


load_metadata <- function(guids=NULL) {
    if (!is.null(guids)) {
        if (guids == TRUE) {
            # guids <- 'intermediate/guids.tsv'
            guids <- 'config/guids.tsv'
        }
        guids <- readr::read_tsv(guids, col_types = cols(
            ID = col_character(),
            run_id = col_character(),
            guid = col_character()
        ))
    }
    samples <- readr::read_tsv('config/samples.tsv', col_types = cols(
      plate = col_double(),
      well = col_character(),
      depth = col_double(),
      expt = col_character(),
      round = col_character(),
      sample = col_character(),
      phage_library = col_character(),
      notes = col_character(),
      ID = col_character()
    ))
    sample_metadata <- readr::read_csv('config/metadata.csv',
                                       col_types = cols(.default = "c")
                                      )
    augment_metadata(samples,
                     guids=guids,
                     sample_metadata=sample_metadata)
}


readRDS.gz <- function(file,threads=(parallel::detectCores())-1) {
    con <- pipe(paste0("pigz -d -c -p",threads," ",shQuote(file)))
    object <- base::readRDS(file = con)
    close(con)
    return(object)
}

saveRDS.gz <- function(object,file,threads=(parallel::detectCores())-1,compression_level=6) {
    con <- pipe(paste0("pigz -c",compression_level," -p",threads," > ",shQuote(file)),"wb")
    saveRDS(object, file = con)
    close(con)
}

format_size <- function(x) {
    return(format(structure(x, class="object_size"), units="auto"))
}

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
