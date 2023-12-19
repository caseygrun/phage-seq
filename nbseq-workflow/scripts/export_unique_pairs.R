log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file,type="message")

library(Biostrings)
library(parallel)
library(digest)


message("Reading merged pairs...")

mergers <- readRDS(snakemake@input[["mergers"]])

message("Read ", length(mergers), " samples")

# if samples have zero merged pairs, their $sequence attribute will have the
# wrong type and cause and error in dplyr::bind_rows
good_samples <- unlist(lapply(mergers, function(m) length(m$sequence)!=0))
mergers <- mergers[good_samples]

message("Kept ", length(mergers), " samples with non-zero merged pairs.")
message("Extracting unique pairs...")

# mergers is a list of data.frames, one per sample
# extract the merged $sequence column from each sample, row-wise concatenate,
# and get unique elements
# uniq_pairs <- unique(do.call(rbind, lapply(mergers, "[", c("sequence"))))
# uniq_pairs <- dplyr::distinct(do.call(rbind, lapply(mergers, "[", c("sequence"))))
pairs <- dplyr::bind_rows(mergers)
rm(mergers)
uniq_pairs <- dplyr::distinct(pairs,sequence)

# assign row names as MD5 hash of concatenated sequence
rownames(uniq_pairs) <- unlist(mclapply(uniq_pairs$sequence, digest::digest, algo="md5", serialize=FALSE))

message("- separating fwd and rev sequences...")

# split each unique pair back into two sequences
uniq_pairs_sep <- tidyr::separate(uniq_pairs,
    col = sequence,
    into = c("fwd", "rev"),
    sep="NNNNNNNNNN",
    fill = "left")

# free memory
rm(pairs,uniq_pairs)


# This is very slow and seemingly not necessary, since you can change the
# expected orientation in bowtie and pysam seems to always report the read query
# sequence as if it were on the (+) strand...
#
# message("- reverse-complementing rev sequences...")
# reverse-complement the "rev" sequence
# uniq_pairs_sep$rev <- lapply(uniq_pairs_sep$rev, dada2::rc)
# uniq_pairs_sep$rev <- unlist(mclapply(uniq_pairs_sep$rev, function(x) as.character(reverseComplement(DNAString(x)))))

message("- writing to CSV file...")

# write the table of unique pairs to a CSV file
write.csv(uniq_pairs_sep, snakemake@output[["unique_pairs"]])

message("Done.")

message("Writing fasta files...")

sink(type="message")
sink()

# also write fwd and rev sequences to separate fasta files for input to bowtie2
write_fasta <- function(seqs, file = NULL, ids = NULL) {
    if(!is.null(file)) {
        sink(file)
    }
    if (is.null(ids)) {
        if (!is.null(names(seqs))) {
            ids <- names(seqs)
        } else {
            ids <- paste0("seq", seq_along(seqs))
        }
    }
    for (i in seq_along(seqs)) {
        cat(">",ids[[i]],"\n", sep="")
        cat(seqs[[i]],"\n",sep="")
    }
    if(!is.null(file)) {
        sink()
    }
}

# write_fasta(uniq_pairs_sep$fwd, snakemake@output[["unique_pairs_fwd"]])
# write_fasta(uniq_pairs_sep$rev, snakemake@output[["unique_pairs_rev"]])
