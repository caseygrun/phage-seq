log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")

library(dada2)
message("Merging denoised paired-end reads and writing feature table with denoised sequences...")

# read in each dereplicated file
message("Reading dereplicated sequences...")
drpFs <- readRDS(snakemake@input[["dereps_fwd"]])
drpRs <- readRDS(snakemake@input[["dereps_rev"]])

# read in denoised files
message("Reading denoised sequences (dadas)...")
ddFs <- readRDS(snakemake@input[["dadas_fwd"]])
ddRs <- readRDS(snakemake@input[["dadas_rev"]])


# MERGE
# =====

message("Merging paired end sequences by aligning overlaps...")
params <- list(
    # keep the rejects in the output at first so we can see why they were rejected
    returnRejects = TRUE,
    verbose = TRUE)
if ("merge" %in% names(snakemake@params) && length(snakemake@params[["merge"]]) > 0) {
    params <- modifyList(params, snakemake@params[["merge"]])
}
print(params)
args <- c(list( dadaF = ddFs, derepF = drpFs, dadaR = ddRs, derepR = drpRs ), params)
mergers <- do.call(mergePairs, args)
out_path_mergers <- snakemake@output[["mergers"]]

message("Merges: ")
print(mergers)
saveRDS(mergers, out_path_mergers)

message("")

# MAKE FEATURE TABLE
# ==================

# now filter out the rejects before making a sequence table
# we gave a list of inputs, so mergers is a list of `data.frame`s, each of which
# has a column $accept = TRUE for merged sequences DADA2 wants to keep
mergers <- lapply(mergers, function(sample_mergers) {
    return(sample_mergers[sample_mergers$accept, , drop = FALSE])
})
seqtab <- makeSequenceTable(mergers)

out_path_with_chimeras <- snakemake@output[["feature_table_with_chimeras"]]
message("Writing output to:")
message(out_path_with_chimeras)
col.names <- names(drpFs)
write.table(t(seqtab), out_path_with_chimeras, sep="\t", row.names=TRUE,
    col.names = col.names)

# REMOVE CHIMERAS
# ===============

message("Removing chimeras")
seqtab_nochim <- removeBimeraDenovo(seqtab, multithread=snakemake@threads)


# filtered.summary<- readRDS(file.path(scratch.path,"filtered/summary.rds"))
#
# track <- cbind(filtered.summary, matrix(0, nrow=nrow(filtered.summary), ncol=3))
# colnames(track) <- c("input", "filtered", "denoised", "merged", "non-chimeric")
# passed.filtering <- track[,"filtered"] > 0
# track[passed.filtering,"denoised"] <- denoisedF
# track[passed.filtering,"merged"] <- rowSums(seqtab)
# track[passed.filtering,"non-chimeric"] <- rowSums(seqtab_nochim)
# write.table(track, out.track, sep="\t", row.names=TRUE, col.names=NA,
# 	    quote=FALSE)

### WRITE OUTPUT AND QUIT ###
# Formatting as tsv plain-text sequence table table
out_path_nochim <- snakemake@output[["feature_table"]]
message("Writing output to:")
message(out_path_nochim)

seqtab_nochim <- t(seqtab_nochim) # QIIME has OTUs as rows
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab_nochim, out_path_nochim, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)

message("Finished writing output.")
sink(type="message")
sink()
