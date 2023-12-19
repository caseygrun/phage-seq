library(dada2)

# read serialized dada objects for each sample, from prior
dds <- lapply(unlist(snakemake@input), readRDS)

# identify sequences with abundance & prevalence matching threshold
st <- dada2::makeSequenceTable(dds)
opts <- dada2::getDadaOpt()
pseudo_priors <- colnames(st)[colSums(st>0) >= opts$PSEUDO_PREVALENCE | colSums(st) >= opts$PSEUDO_ABUNDANCE]
rm(st)

# save priors
saveRDS(pseudo_priors, snakemake@output[[1]])
