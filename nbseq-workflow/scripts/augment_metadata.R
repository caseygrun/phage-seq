library(tidyverse)
library(magrittr)

snakemake@source('common.R')

metadata <- load_metadata()

ext = tools::file_ext(snakemake@output[[1]])
if (ext == 'csv') {
	sep = ","
} else {
	sep = "\t"
}

metadata %>% write_delim(snakemake@output[[1]], quote="all", delim=sep)
