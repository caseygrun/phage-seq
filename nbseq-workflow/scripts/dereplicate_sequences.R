log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")


library(dada2)

for (direction in c("fwd", "rev")) {
	message(paste("Dereplicating",direction,"reads..."))
	# save each file as a 1-item list, with the name given as the sample filename
	x = list(
			dada2::derepFastq(
				snakemake@input[[direction]], verbose=TRUE)
		)
	names(x) <- basename(snakemake@input[[direction]])
	saveRDS(x, file = snakemake@output[[direction]])
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
