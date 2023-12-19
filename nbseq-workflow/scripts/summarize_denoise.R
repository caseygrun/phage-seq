snakemake@source('common.R')
# log.file<-file(snakemake@log[[1]],open="wt")
# sink(log.file)
# sink(log.file,type="message")

library(dada2)
library(Matrix)

ft <- read_feature_table_sparse(snakemake@input[["feature_table"]])

summary_path <- snakemake@output[["summary"]]
message("Writing summary of denoised reads per sample to ", summary_path)
denoise_summary <- summarize_samples(ft, name = 'denoised')
write_tsv(denoise_summary, summary_path)


# # get list of summaries (one per sample), then stack them into one data.frame
# filtered.summary <- do.call(rbind, lapply(snakemake@input[["filter_summaries"]], readRDS))
# write.table(filtered.summary, snakemake@output[["summary_filter"]], sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#
#
# ddFs <- readRDS(snakemake@input[["dadas_fwd"]])
# ddRs <- readRDS(snakemake@input[["dadas_rev"]])
# seqtab <- t(read.table(snakemake@input[["feature_table_with_chimeras"]], header=TRUE, row.names=1))
# seqtab.nochim <- t(read.table(snakemake@input[["feature_table"]], header=TRUE, row.names=1))
#
# getN <- function(x) sum(getUniques(x))
#
# track <- cbind(
# 	# filtered.summary, # input, filtered
# 	sapply(ddFs, getN),
# 	sapply(ddRs, getN),
# 	rowSums(seqtab), #sapply(mergers, getN),
# 	rowSums(seqtab.nochim))
# colnames(track) <- c(
# 	#"input", "filtered",
# 	"denoisedF", "denoisedR", "merged", "nonchim"
# )
#
# write.table(track, snakemake@output[["summary_denoise"]], sep="\t", row.names=TRUE, col.names=NA,
# 	    quote=FALSE)
#
#
# 		# denoisedF <- rep(0, length(filtsF))
# 		# denoisedF[[j]] <- getN(ddF)
#
# 		# track <- cbind(filtered.summary, matrix(0, nrow=nrow(filtered.summary), ncol=3))
# 		# colnames(track) <- c("input", "filtered", "denoised", "merged", "non-chimeric")
# 		# passed.filtering <- track[,"filtered"] > 0
# 		# track[passed.filtering,"denoised"] <- denoisedF
# 		# track[passed.filtering,"merged"] <- rowSums(seqtab)
# 		# track[passed.filtering,"non-chimeric"] <- rowSums(seqtab.nochim)
#
# # sink(type="message")
# # sink()
