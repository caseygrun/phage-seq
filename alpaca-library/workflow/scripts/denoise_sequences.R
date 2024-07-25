log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")


library(dada2)

params <- list(
    multithread = snakemake@threads,
	verbose = TRUE,
	pool = "pseudo"
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

# read in all dereplicated files from single RDS, so can do pseudopooling
# TODO: may rewrite this, implementing pseudopooling explicitly, so
# I can stream in the read sequences and reduce memory requirements
message("Loading dereplicated sequences...")
drpFs <- readRDS(snakemake@input[["dereps_fwd"]])
drpRs <- readRDS(snakemake@input[["dereps_rev"]])


# perform denoising with pseudo-pooling
message("Denoising forward reads...")
saveRDS(
	do.call(dada2::dada, c(
		list(derep = drpFs, err = errF),
		params
	)),
	snakemake@output[["fwd"]]
)

message("")
message("Denoising reverse reads...")
saveRDS(
	do.call(dada2::dada, c(
		list(derep = drpRs, err = errR),
		params
	)),
	snakemake@output[["rev"]]
)

message("Finished denoising.")
sink(type="message")
sink()


# for(j in seq(length(filtsF))) {
# 	drpF <- readRDS(derepsF[[j]])
# 	drpR <- readRDS(derepsR[[j]])
# 	ddF <- dada(drpF, err=errF, multithread=multithread, verbose=verbose)
# 	ddR <- dada(drpR, err=errR, multithread=multithread, verbose=verbose)
# }
