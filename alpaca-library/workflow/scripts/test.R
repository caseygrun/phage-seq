

print(snakemake@params)

cat("------\n")

print(names(snakemake@params))

cat("------\n")

print(snakemake@params[["merge"]])


cat("------\n")

params <- list(
    # keep the rejects in the output at first so we can see why they were rejected
    returnRejects = TRUE,
    verbose = TRUE
)
if ("merge" %in% names(snakemake@params) && length(snakemake@params[["merge"]]) > 0) {
    params <- modifyList(params, snakemake@params[["merge"]])
}

print(params)
