log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")


library(dada2)

for (direction in c("fwd", "rev")) {
	message(paste("Dereplicating",direction,"reads..."))
	saveRDS(
		dada2::derepFastq(
			snakemake@input[[direction]], verbose=TRUE),
		snakemake@output[[direction]])
}

sink(type="message")
sink()

# for(j in seq(length(filtsF))) {
# 	drpF <- dada2::derepFastq(filtsF[[j]])
# 	saveRDS(drpF,derepsF[[j]])
# }
#
# for(j in seq(length(filtsR))) {
# 	drpR <- dada2::derepFastq(filtsR[[j]])
# 	saveRDS(drpR,derepsR[[j]])
# }
