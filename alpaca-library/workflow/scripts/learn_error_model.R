log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")


library(dada2)
library(ggplot2)

# LEARN ERROR RATES
# -----------------

params <- list(
    qualityType = "Auto",
    multithread = snakemake@threads,
    verbose = TRUE,
    randomize = TRUE
)

# take extra parameters, e.g. randomize=True
if (length(snakemake@params) > 0) {
    params <- modifyList(params, snakemake@params[ names(snakemake@params) != "" ])
}

args <- c(
    list(fls = snakemake@input[["fwd"]]),
    params
)
print(args)
errF <- do.call(dada2::learnErrors, args)
saveRDS(errF, snakemake@output[["error_fwd"]])
err_plot <- dada2::plotErrors(errF, nominalQ = TRUE)
ggsave(snakemake@output[["plot_fwd"]], err_plot, width = 8, height = 8, dpi = 300)


errR <- do.call(dada2::learnErrors, c(
    list(fls = snakemake@input[["rev"]]),
    params
))
saveRDS(errR, snakemake@output[["error_rev"]])
err_plot <- dada2::plotErrors(errR, nominalQ = TRUE)
ggsave(snakemake@output[["plot_rev"]], err_plot, width = 8, height = 8, dpi = 300)


sink(type="message")
sink()
