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

# message("Denoising parameters: ")
# print(params)
#
# message("Loading error structure...")
err <- readRDS(snakemake@input[["error"]])

if ("priors" %in% snakemake@input) {
	priors <- readRDS(snakemake@input[["priors"]])
	message("Using priors (pseudo-pooling)")
} else {
	priors <- NULL
}

message("Denoising: ", snakemake@input[["derep"]])
drp <- readRDS(snakemake@input[["derep"]])
saveRDS(do.call(dada2::dada, c(
	list(derep = drp, err = err, priors=priors),
	params
)),snakemake@output[["dada"]])
