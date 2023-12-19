feature_table_with_chimeras <- read.csv(snakemake@input[['feature_table_with_chimeras']], sep='\t',
	stringsAsFactors=FALSE, check.names=FALSE)
asv_seqs <- read.csv(snakemake@input[['asv_sequences']],
	stringsAsFactors=FALSE, check.names=FALSE)

# replace md5 hashes with sequences merged by reference
colnames(feature_table_with_chimeras) <- asv_seqs[colnames(feature_table_with_chimeras),'seq']

feature_table <- dada2::removeBimeraDenovo(feature_table_with_chimeras)

write.table(feature_table, snakemake@output[['feature_table']], sep="\t")
