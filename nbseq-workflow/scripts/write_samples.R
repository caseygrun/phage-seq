# temporary script to rewrite sample CSV file from 202201
readRDS.gz <- function(file,threads=(parallel::detectCores())-1) {
    con <- pipe(paste0("pigz -d -c -p",threads," ",shQuote(file)))
    object <- base::readRDS(file = con)
    close(con)
    return(object)
}

mergers <- readRDS.gz("intermediate/20220111/denoise/pairs.rds")
cat("\nRead merged sequences...\n")
print(names(mergers))
cat("\nWriting to feature table...\n")

write_feature_table_sparse <- function(dir, samples = NULL, features = NULL) {
    write.table(samples, file.path(dir, 'samples.tsv'), sep="\t", row.names = FALSE, col.names = FALSE)
}

write_feature_table_sparse(
	dir = "intermediate/20220111/denoise/feature_table/",
	samples = names(mergers)
)
