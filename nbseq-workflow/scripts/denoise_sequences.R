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
    params <- modifyList(params, snakemake@params[ names(snakemake@params) != "" ])
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
				message("-> intermediate results will be read from and saved to to ", cache_folder)
			} else {
				message("-> intermediate results will be saved to ", cache_folder, " but existing data in this folder will be overwritten.")
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
			message(derep_files[[i]])

			# load from cache if available
			if (!is.null(cache_folder)) {
				cache_path <- file.path(cache_folder, paste0(names(dds)[[i]], '.rds'))

				if (!invalidate && file.exists(cache_path)) {
					dds[[i]] <- readRDS(cache_path)
					message('-> reading from cache: ', cache_path)
					pb$tick(1)
					next
				}
			}

			drp <- readRDS(derep_files[[i]])

            if (partition) {
                ss <- sample(1:n_partitions, size=length(drp[[1]]$map), replace=TRUE)
                maps <- split(drp[[1]]$map, ss)
                drps <- vector(mode='list', length=n_partitions)
                for (j in seq_along(drps)) {
                    drps[[j]] <- as(
                        list(uniques = drp[[1]]$uniques, qualities = drp[[1]]$qualities, map = maps[[j]]),
                        "derep")
                }
            }

			# message("-> ", sum(drp$uniques), " reads in ", length(drp$uniques), " unique sequences.")
			dds[[i]] <- do.call(dada2::dada, c(
				list(derep = drp, err = err, priors=priors),
				params
			))

			if (!is.null(cache_folder)) {
				saveRDS(dds[[i]], cache_path)
				message('-> writing to cache: ', cache_path)
			}
			pb$tick(1)
		}

		return(dds)
	}

	cache_folder <- NULL
	if (!is.null(output_folder)) {
		cache_folder_initial <- mkdirp(output_folder, 'initial')
		cache_folder_pseudo <- mkdirp(output_folder, 'pseudo')
		message("Will attempt to cache intermediate results to: ", output_folder)
	}

	# perform first pass
	message("Performing first pass denoising...")
	dds <- denoise(derep_files, err, params, cache_folder=cache_folder_initial)
	message("... Done.\n")

	# identify priors
	message("Identifying prior sequences...")
	st <- dada2::makeSequenceTable(dds)
	opts <- dada2::getDadaOpt()
	pseudo_priors <- colnames(st)[colSums(st>0) >= opts$PSEUDO_PREVALENCE | colSums(st) >= opts$PSEUDO_ABUNDANCE]
	rm(st)
	message("Identified ", length(pseudo_priors), " prior sequences.\n")


	# repeat denoising with priors
	message("Repeating denoising with prior information (pseudo-pooling)...")
	dds <- denoise(derep_files, err, params, priors=pseudo_priors, cache_folder=cache_folder_pseudo)
	message("... Done.")
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
message("Done.")

message("Cleaning up cache folders: ", cache_folders)
unlink(cache_folders, recursive=TRUE)

snakemake_log_end()
