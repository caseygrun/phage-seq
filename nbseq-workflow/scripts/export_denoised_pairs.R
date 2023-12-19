snakemake@source('common.R')
snakemake_log_start(snakemake)

library(dada2)
library(parallel)

snakemake_use_threads(snakemake)
message("Merging denoised paired-end reads and writing feature table with denoised sequences...")
message("Streaming dereplicated and denoised data...")


# MERGE
# =====

# ask dada2 for denoised, paired-end sequences. It is kind of confusing to
# extract these sequences directly from the dada-class and derep-class objects,
# so we'll have dada2 do it.
# the `justConcatenate` parameter will cause the fwd and rc(rev) sequences to
# be joined by "NNNNNNNNNN"; we will split them back out and undo the rc after.
message("Merging paired end sequences by concatenating")
params <- list(
	justConcatenate = TRUE,
	verbose = TRUE)
if ("merge" %in% names(snakemake@params) && length(snakemake@params[["merge"]]) > 0) {
	params <- modifyList(params, snakemake@params[["merge"]])
}
print(params)
# args <- c(list( dadaF = ddFs, derepF = drpFs, dadaR = ddRs, derepR = drpRs ), params)
# mergers <- do.call(mergePairs, args)

#mergers <- vector(mode="list", length=length(ddFs))

basenames <- tools::file_path_sans_ext(basename(snakemake@input[["dereps_fwd"]]))

filenames <- data.frame(
	dereps_fwd = snakemake@input[["dereps_fwd"]],
	dereps_rev = snakemake@input[["dereps_rev"]],

	# hack: snakemake gives intermediate/{run}/denoise/fwd.rds, we really look in
	# intermediate/{run}/denoise/fwd/initial/{name}.rds
	# doing this because the denoise step runs on all samples at once,
	# which seemed like a good idea at one time...
	dadas_fwd  = paste0(file.path(tools::file_path_sans_ext(snakemake@input[["fwd"]]),"initial",basenames),".rds"),
	dadas_rev  = paste0(file.path(tools::file_path_sans_ext(snakemake@input[["rev"]]),"initial",basenames),".rds"),


	name = basenames
)
message("Filenames: ")
options(width=1000)
print(filenames)

# print(filenames)
message("Merging sequences from ",nrow(filenames)," samples...")
message("Using parallel processing with ", snakemake@threads, " threads...")

merge_seqs <- function(i) {

	bn_drpF <- tools::file_path_sans_ext(basename(filenames[i, 'dereps_fwd']))
	bn_ddsF <- tools::file_path_sans_ext(basename(filenames[i, 'dadas_fwd']))
	stopifnot(bn_drpF == bn_ddsF)

	# drpF and dprR are 1-element lists; the single element is a list of
	# derep-class objects. Remove this outermost level of list
	tryCatch({

		drpF <- readRDS(filenames[i, 'dereps_fwd'])[[1]] #unlist(readRDS(filenames[i, 'dereps_fwd']), recursive = FALSE)
		drpR <- readRDS(filenames[i, 'dereps_rev'])[[1]] #unlist(readRDS(filenames[i, 'dereps_rev']), recursive = FALSE)

		# dadaF and dadaR are either `dada` objects or `list`s of `dada` objects
		dadaF <- readRDS(filenames[i, 'dadas_fwd'])
		dadaR <- readRDS(filenames[i, 'dadas_rev'])

	}, error = function() {
		sink(stderr())
		on.exit(sink(NULL))
		message(paste("BAD:  In sample", name, "unable to read dereps_fwd, dereps_rev, dadas_fwd, or dadas_rev: "))
		# traceback(3, max.lines = 1L)
		traceback()
		# q(status=1)
	})

	name <- filenames[i, 'name']

	# HACKS FOR 20220111:
	# if either dadaF or dadaR was denoised on partition but not both, it would be hard to put them back together, so I'll just make a list and re-denoise those together
	if ((class(dadaF) != class(dadaR)) || (length(dadaF) != length(dadaR))) {
		message("BAD:  In sample ", name , " , lengths of dadaF and dadaR do not agree: ", paste(class(dadaF), "=", length(dadaF), class(dadaR), "=", length(dadaR)))
	}

	# if sample was denoised on an unpartitioned sample and thus stored not as a list
	if (
		((class(dadaF) == "dada") || (class(dadaF) == "list" && length(dadaF) == 1)) || 
		((class(dadaR) == "dada") || (class(dadaR) == "list" && length(dadaR) == 1))) {

		# message("INFO: In sample ", name, " , loading unpartitioned dereps.")
		# tryCatch({

		# 	# HACK for 20220111:
		# 	# find the original unpartitioned filtered_trimmed_dereplicated sequence
		# 	# partitioned: intermediate/20220111/filtered_trimmed_dereplicated/{fwd,rev}/sample.rds
		# 	# unpartitioned: intermediate/20220111/filtered_trimmed_dereplicated/{fwd,rev}_unpartitioned/sample.rds
		# 	drpF <- readRDS(file.path(paste0(dirname(filenames[i, 'dereps_fwd']), "_unpartitioned"), paste0(name,".rds")))
		# 	drpR <- readRDS(file.path(paste0(dirname(filenames[i, 'dereps_rev']), "_unpartitioned"), paste0(name,".rds")))

		# }, error = function() {
		# 	sink(stderr())
		# 	on.exit(sink(NULL))
		# 	message(paste("BAD:  In sample", name, "unable to read dereps_fwd_unpartitioned or dereps_rev_unpartitioned: "))
		# 	traceback()
		# 	# q(status=1)
		# })

		if ((class(dadaF) == "dada")) {
			message("INFO: In sample ", name, " , wrapping dadaF in list.")
			dadaF <- list(dadaF)
		}
		if ((class(dadaR) == "dada")) {
			message("INFO: In sample ", name, " , wrapping dadaR in list.")
			dadaR <- list(dadaR)
		}
	}

	if (length(unique(sapply(list(dadaF, dadaR, drpF, drpR), length))) != 1) {
		stop(paste("BAD:  In sample", name, "lengths of dadaF, dadaR, drpF, drpR are not equal:", length(dadaF),  length(dadaR), length(drpF), length(drpR)))
	}

	mrgs <- do.call(dada2::mergePairs, c(
			list(
				dadaF = dadaF,
				derepF = drpF,
				dadaR = dadaR,
				derepR = drpR ),
			params)
		)

	# if the sample was partitioned (e.g. dadaF, derepF, etc. are lists,
	# this will return a list instead of a data.frame. Bind elements of the
	# list into a single data.frame , use dplyr::bind_rows because it's faster
	# than do.call(rbind, ...)
	if (!is.data.frame(mrgs) && is.list(mrgs)) {
		# saveRDS(mrgs,snakemake@output[["mergers"]])
		message("INFO: In sample ", name, " , merging pairs from ", length(mrgs), " partitions into one data.frame.")
		# print(mrgs)
		mrgs <- dplyr::bind_rows(mrgs)
	} else if (!is.data.frame(mrgs)) {
		message("BAD:  In sample ", name, " dada2::mergePairs yielded object of class ", class(mrgs), "; ", head(str(mrgs)))
		stop()
	}
	message("OK:   ", name)
	rm(drpF, drpR, dadaF, dadaR)
	return(mrgs)
}

gc()

# merge all sequences using multiple cores
mergers <- mclapply(seq_len(nrow(filenames)), merge_seqs)

mergers_classes <- lapply(mergers,class)
if (!all(mergers_classes == "data.frame")) {
	message(sum(mergers_classes != "data.frame"), " of ", length(mergers_classes), " samples did not result properly.")

	# print(mergers_classes)
	if (any(mergers_classes == "try-error")) {
		errors <- mergers[mergers_classes != "data.frame"]
		print(as.character(errors))
		nerrors <- min(c(10, length(errors)))
		for (i in 1:nerrors) {
			traceback(errors[[i]])
			cat()
		}

		stop()
	}
}

message("Merged sequences from ",length(mergers)," samples.")
# print(mergers)
# print(length(mergers))
names(mergers) <- filenames$name


out_path_mergers <- snakemake@output[["mergers"]]
# message("Merges: ")
# print(mergers)
message("Saving merged sequences to ",out_path_mergers, " ...")
saveRDS(mergers, out_path_mergers, compress = FALSE)
message("Saved merged sequences.")

message("")

snakemake_log_end()
