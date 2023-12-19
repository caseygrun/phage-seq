log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")


library(dada2)

params <- list(
    multithread = snakemake@threads,
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
errF <- readRDS(snakemake@input[["error_fwd"]])
errR <- readRDS(snakemake@input[["error_rev"]])

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

dada_pseudo_pool <- function(derep_files, err, params, monotonic=TRUE) {

	if (monotonic) {
		err <- enforce_monotonicity(err)
	}
	denoise <- function(derep_files, err, params, priors=character(0)) {
		dds <- vector("list", length(derep_files))
		names(dds) <- gsub('.rds','',basename(derep_files))

		for (i in seq_along(derep_files)) {
			message(derep_files[[i]])
			drp <- readRDS(derep_files[[i]])
			dds[[i]] <- do.call(dada2::dada, c(
				list(derep = drp, err = err, priors=priors),
				params
			))
		}

		return(dds)
	}

	# perform first pass
	dds <- denoise(derep_files, err, params)

	# identify priors
	st <- dada2::makeSequenceTable(dds)
	opts <- dada2::getDadaOpt()
	pseudo_priors <- colnames(st)[colSums(st>0) >= opts$PSEUDO_PREVALENCE | colSums(st) >= opts$PSEUDO_ABUNDANCE]
	rm(st)

	# repeat denoising with priors
	dds <- denoise(derep_files, err, params, priors=pseudo_priors)
	return(dds)
}

message("Denoising forward reads...")
saveRDS(
	dada_pseudo_pool(snakemake@input[["dereps_fwd"]], errF, params),
	snakemake@output[["fwd"]]
)

message("Denoising reverse reads...")
saveRDS(
	dada_pseudo_pool(snakemake@input[["dereps_rev"]], errR, params),
	snakemake@output[["rev"]]
)


# # read in each dereplicated sample from a separate RDS file, but all at once so
# # I can do pseudopooling.
# # each RDS contains a 1-element list, with name = sample name, and
# # value = derep-class.
# # TODO: may rewrite this, implementing pseudopooling explicitly, so
# # I can stream in the read sequences and reduce memory requirements
# message("Loading dereplicated sequences...")
# drpFs <- do.call(c, lapply(snakemake@input[["dereps_fwd"]], readRDS)) #readRDS(snakemake@input[["dereps_fwd"]])
#
# # perform denoising with pseudo-pooling
# message("Denoising forward reads...")
# saveRDS(
# 	do.call(dada2::dada, c(
# 		list(derep = drpFs, err = enforce_monotonicity(errF)),
# 		params
# 	)),
# 	snakemake@output[["fwd"]]
# )
# rm(drpFs)
# message("Done.")
#
#
# drpRs <- do.call(c, lapply(snakemake@input[["dereps_rev"]], readRDS)) #readRDS(snakemake@input[["dereps_rev"]])
# message("")
# message("Denoising reverse reads...")
# saveRDS(
# 	do.call(dada2::dada, c(
# 		list(derep = drpRs, err = enforce_monotonicity(errR)),
# 		params
# 	)),
# 	snakemake@output[["rev"]]
# )

message("Done.")
message("Finished denoising.")
sink(type="message")
sink()
