# snakemake@source('common.R')
library(readr)
library(magrittr)
library('dada2')


infiles <- unlist(snakemake@input)

dereps <- sapply(infiles, function(derep_path) {
	derep <- readRDS(derep_path)[[1]]
	unlist(length(derep$uniques))
})

df <- data.frame(
	sample = gsub('.rds','', basename(infiles)),
	derep = dereps
)
write_tsv(df, snakemake@output[[1]])
