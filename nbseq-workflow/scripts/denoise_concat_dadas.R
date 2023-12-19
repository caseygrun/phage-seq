# snakemake@source('common.R')

library(dada2)

dds <- lapply(unlist(snakemake@input), readRDS)
names(dds) <- gsub('.rds','',basename(snakemake@input["dadas"]))

saveRDS(dds, snakemake@output[[1]])
