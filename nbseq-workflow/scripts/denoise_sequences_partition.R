snakemake@source('common.R')
snakemake_log_start(snakemake)

library(dada2)
library(RcppParallel)
library(tools)

params <- list(
	multithread = snakemake@threads,
	# verbose = TRUE,
	verbose = TRUE,
	pool = FALSE #pool = "pseudo"
)

# take extra parameters, e.g. maxEE, truncQ, minLen, etc.
if (length(snakemake@params) > 0) {
	params <- modifyList(params, snakemake@params[ names(snakemake@params) != "" & names(snakemake@params) != "merge" ])
}

message("Denoising parameters: ")
print(params)

message("Loading error structure...")
err <- readRDS(snakemake@input[["error"]])

# functions to enforce monotonicity of error function; necessary for binned
# quality scores from NovaSeq data. see https://github.com/benjjneb/dada2/issues/791
# estimated error should increase as Q-score decreases, but we prefer estimates
# from the higher Q-scores, hence `monotonic_up(..., direction='reverse')`.
# getErrors(err) is 16 x 40 matrix, rows are transitions (A2A, A2C, etc.),
# columns are Q-scores.
monotonic <- function(lst, compare=max, direction='forward') {
	new <- numeric(length(lst))
	names(new) <- names(lst)

	order <- seq_along(lst)
	if (direction == 'reverse') {
		order <- rev(order)
	}

	bar <- lst[[order[[1]]]]
	for (ii in order) {
		bar <- compare(bar, lst[[ii]])
		new[[ii]] <- bar
	}
	return(new)
}
monotonic_up <- function(lst,...) { return(monotonic(lst, compare=max,...)) }
monotonic_down <- function(lst,...) { return(monotonic(lst, compare=min,...)) }
enforce_monotonicity <- function(err) {
	return( t(apply(getErrors(err), 1, monotonic_up, direction='reverse')) )
}

# partition_dereps <- function(dereps, n_partitions) {
#
# 	chunk_size = floor(length(dereps$map) / n_partitions)
#
# 	new_dereps <- map(seq_len(n_partitions), function(i) {
# 		derepMap = dereps$map[(i-1)*chunk_size+1 : i*chunk_size]
# 		derepO <- list(uniques=dereps$uniques, quals=dereps$quals, map=derepMap)
# 		derepO <- as(derepO, "derep")
# 	})
#
# 	return new_dereps
# }
#
# denoise_with_partition <- function (drp, err, priors, max_partition_size = 100000) {
# 	n_partitions <- ceiling(length(drp$map) / max_partition_size)
# 	drp_partitions <- partition_dereps(drp, n_partitions)
#
# }


cache_folders <- character(0)

dada_pseudo_pool <- function(derep_files, err, params,
	monotonic=TRUE, output_folder=NULL, invalidate=FALSE) {

	if (monotonic) {
		err <- enforce_monotonicity(err)
	}
	denoise <- function(derep_files, err, params, priors=character(0), cache_folder=NULL, invalidate=FALSE) {
		if (!is.null(cache_folder)) {
			if (!invalidate) {
				tmessage("-> intermediate results will be read from and saved to to ", cache_folder)
			} else {
				tmessage("-> intermediate results will be saved to ", cache_folder, " but existing data in this folder will be overwritten.")
			}
			cache_folders <- c(cache_folders, cache_folder)
		}

		dds <- vector("list", length(derep_files))
		names(dds) <- gsub('.rds','',basename(derep_files))

		pb <- progress_bar_file$new(total = length(dds),
			format = "[:bar] :current/:total (:percent), :tick_rate samples/sec :elapsed elapsed, eta: :eta \n",
			force = TRUE
		)
		pb$tick(0)

		for (i in seq_along(derep_files)) {
			tmessage(derep_files[[i]])

			# load from cache if available
			if (!is.null(cache_folder)) {
				cache_path <- file.path(cache_folder, paste0(names(dds)[[i]], '.rds'))

				if (!invalidate && file.exists(cache_path)) {
					dds_i <- readRDS(cache_path)
					tmessage('-> reading from cache: ', cache_path)
					if (class(dds_i) == 'list') {
						dds[[i]] <- dds_i
					} else {
						tmessage('---> wrapping single dada result in list')
						dds[[i]] <- list(dds_i)
					}
					pb$tick(1)
					next
				}
			}

			sample_start_time <- Sys.time()

			# read dereplicated sequences
			drp <- readRDS(derep_files[[i]])
			if (class(drp[[1]]) == 'list') {
				drp_partitions <- drp[[1]]
				n_partitions = length(drp_partitions)
				if (n_partitions > 1) {
					tmessage("-> divided into ", n_partitions, " partitions...")
				}
				partition_dds <- vector(mode="list", length=length(drp_partitions))

				partition_start_time <- Sys.time()
				for (j in seq_along(drp_partitions)) {
					if (n_partitions > 1) {
						tmessage("---> denoising partition ", j, " of ", n_partitions, "...")
					}
					partition_dds[[j]] <- do.call(dada2::dada, c(
						list(derep = drp_partitions[[j]], err = err, priors = priors),
						params
					))
				}
				dds[[i]] <- partition_dds

			} else {

				# tmessage("-> ", sum(drp$uniques), " reads in ", length(drp$uniques), " unique sequences.")
				dds[[i]] <- list(
					do.call(dada2::dada, c(
						list(derep = drp, err = err, priors=priors),
						params
					))
				)
			}
			if (!is.null(cache_folder)) {
				saveRDS(dds[[i]], cache_path)
				tmessage('-> writing to cache: ', cache_path)
			}

			sample_time <- difftime(Sys.time(), sample_start_time, "auto")
			tmessage("-> done (", format(sample_time, format = "auto", digits = 2), ")")

			pb$tick(1)
		}

		return(dds)
	}

	cache_folder <- NULL
	if (!is.null(output_folder)) {
		cache_folder_initial <- mkdirp(output_folder, 'initial')
		cache_folder_pseudo <- mkdirp(output_folder, 'pseudo')
		tmessage("Will attempt to cache intermediate results to: ", output_folder)
	}

	# perform first pass
	tmessage("Performing first pass denoising...")
	dds <- denoise(derep_files, err, params, cache_folder=cache_folder_initial)
	tmessage("... Done.\n")

	# # identify priors
	# tmessage("Identifying prior sequences...")
	#
	# # dds is a list of length (# samples)
	# st <- dada2::makeSequenceTable(unlist(dds,recursive=FALSE))
	# opts <- dada2::getDadaOpt()
	# pseudo_priors <- colnames(st)[colSums(st>0) >= opts$PSEUDO_PREVALENCE | colSums(st) >= opts$PSEUDO_ABUNDANCE]
	# rm(st)
	# tmessage("Identified ", length(pseudo_priors), " prior sequences.\n")


	# repeat denoising with priors
	# tmessage("Repeating denoising with prior information (pseudo-pooling)...")
	# dds <- denoise(derep_files, err, params, priors=pseudo_priors, cache_folder=cache_folder_pseudo)
	tmessage("... Done.")
	tmessage("Writing to output file...")
	return(dds)
}

if ("invalidate" %in% snakemake@params) {
	invalidate <- snakemake@params[["invalidate"]]
} else {
	invalidate <- FALSE
}

message("Denoising ", snakemake@wildcards[["direction"]], " reads...")
saveRDS(
	dada_pseudo_pool(snakemake@input[["dereps"]], err, params,
		output_folder=file_path_sans_ext(snakemake@output[[1]]),invalidate=invalidate),
	snakemake@output[[1]]
)
message("... Done.")

# message("Cleaning up cache folders: ", cache_folders)
# unlink(cache_folders, recursive=TRUE)

snakemake_log_end()
