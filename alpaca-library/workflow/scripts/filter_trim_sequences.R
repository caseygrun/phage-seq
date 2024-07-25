log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")


library(methods)
library(dada2)
library(ggplot2)

#
# unfiltered_files_fwd <- snakemake@input[["fwd"]]
# unfiltered_files_rev <- snakemake@input[["rev"]]
#
# filtered_files_fwd <- snakemake@input[["fwd"]]
# filtered_files_rev <- snakemake@input[["rev"]]

params <- list(
    multithread = snakemake@threads,
    verbose = TRUE
)

# take extra parameters, e.g. maxEE, truncQ, minLen, etc.
if (length(snakemake@params) > 0) {
    params <- modifyList(params, snakemake@params[ names(snakemake@params) != "" ])
}

args <- c(
    list(
        fwd      = snakemake@input[["fwd"]],
        filt     = snakemake@output[["fwd"]],
        rev      = snakemake@input[["rev"]],
        filt.rev = snakemake@output[["rev"]]
    ),
    params
)

message("Filtering and trimming with dada2... Options: ")
message(args)

filtered_summary <- do.call(dada2::filterAndTrim, args)

# filtered.summary <- filterAndTrim(
#     unfiltered_files_fwd, filtered_files_fwd,
#     unfiltered_files_rev, filtered_files_rev,
#     # No truncLen for variable length sequences.
#     # truncLen=c(truncLenF, truncLenR),
#     trimLeft=c(trimLeftF, trimLeftR),
#     maxEE=maxEE, truncQ=truncQ, minLen=minLen,
#     rm.phix=TRUE,
#     multithread=multithread,
#     verbose=verbose
# )

saveRDS(filtered_summary,snakemake@output[["summary"]])


pquality <- plotQualityProfile(unlist(snakemake@input[["fwd"]]))
ggsave(snakemake@output[["plot_fwd"]], pquality, width = 4, height = 3, dpi = 300)
pquality <- plotQualityProfile(unlist(snakemake@input[["rev"]]))
ggsave(snakemake@output[["plot_rev"]], pquality, width = 4, height = 3, dpi = 300)


# cat(ifelse(file.exists(filtsF), ".", "x"), sep="")
# filtsF <- list.files(filtered.dirF, pattern=".fastq.gz$", full.names=TRUE)
# filtsR <- list.files(filtered.dirR, pattern=".fastq.gz$", full.names=TRUE)
# cat("\n")
# if(length(filtsF) == 0) { # All reads were filtered out
#   errQuit("No reads passed the filter (were truncLenF/R longer than the read lengths?)", status=2)
# }

sink(type="message")
sink()
