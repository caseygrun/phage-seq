snakemake@source('common.R')
snakemake_log_start(snakemake)


library(dada2)

max_partition_size <- as.numeric(snakemake@params[["max_partition_size"]])
# check <- TRUE
check <- FALSE

# n_partitions <- snakemake@params[["n_partitions"]]

# load and dereplicate FASTQ files
drp <- list("fwd" = NA, "rev" = NA)

for (direction in c("fwd", "rev")) {
	message("Loading ", direction, " reads from ", snakemake@input[[direction]], " ...")
	drp[[direction]] <- dada2::derepFastq(snakemake@input[[direction]], verbose=TRUE)
	print(drp[[direction]])
	message()
}

# if sample has more than `max_partition_size` reads, split it into partitions
# up to `max_partition_size` reads
drpF <- drp[["fwd"]]
n_reads <- length(drpF$map)
n_partitions <- ceiling(n_reads / max_partition_size)
message("Dividing ", n_reads, " into ", n_partitions, " partitions.")

# create containers to hold partitions
drp_partitions <- list(
	"fwd" = vector(mode='list', length=n_partitions),
	"rev" = vector(mode='list', length=n_partitions)
)

# drpF$map is an integer vector, with length = number of reads; values are
# indices into drpF$uniques.
# ss is an integer vector with length = number of reads, values = partition
# number 1..n_partitions
# mapsF[[j]] = values of drpF$map assigned to partition j
ss <- sample(1:n_partitions, size=length(drpF$map), replace=TRUE)
maps <- list(
	"fwd" = split(drp[["fwd"]]$map, ss),
	"rev" = split(drp[["rev"]]$map, ss)
)
message("Partition sizes: ", paste(lapply(maps[["fwd"]], length), collapse=", "))

# construct new derep-class objects containing the partitioned sequences
for (direction in c("fwd", "rev")) {
	message("Partitioning ", direction, " reads...")
	for (j in seq_len(n_partitions)) {

	    # re-tabulate how many reads were assigned to each unique sequence in
	    # the partition. replace names() with the original unique sequences
	    uniq_counts <- tabulate(maps[[direction]][[j]], nbins = length(drp[[direction]]$uniques))
	    names(uniq_counts) <- names(drp[[direction]]$uniques)

		# remove elements from uniq_counts that did not appear in this
		# partition, and re-number the indices in the map for this partition so
		# they correspond with the new values of uniq_counts.

		# logical vector, where uniq_in_partition[[x]] = TRUE if unique sequence
		# x is in this partition
	    uniq_in_partition <- (uniq_counts > 0)

		# integer vector mapping indices in drp[[direction]]$uniques to indices in
		# uniq_counts
	    uniq_map <- integer(length(uniq_counts))
	    uniq_map[uniq_in_partition] <- 1:sum(uniq_in_partition)

		# re-map indices in $map for this partition to the new indices
		map <- uniq_map[ maps[[direction]][[j]] ]

		# filter uniq_counts to include only those that appearedi in this partition
		uniq_counts <- uniq_counts[uniq_counts > 0]

		# filter quality score matrix to only include unique sequences that
		# actually appear in partition, and truncate to length of longest
		# sequence seen
		maxlen <- max(nchar(names(uniq_counts)))
		quals <- drp[[direction]]$quals[names(uniq_counts), 1:maxlen, drop=FALSE]

		drp_partitions[[direction]][[j]] <- as(
			list(uniques = uniq_counts, quals = quals, map = map),
			"derep"
		)

		# check that we got everything right
		if (check) {
			stopifnot(nrow(quals) == length(uniq_counts))
			stopifnot(ncol(quals) == maxlen)
			stopifnot(all(names(drp[[direction]]$uniques)[maps[[direction]][[j]]] == names(uniq_counts)[map]))
			message(" - OK")
		}
	}

	# make a single-element list with key = filename, value = derep-class
	drp_direction_out <- list( drp_partitions[[direction]] )
	names(drp_direction_out) <- basename(snakemake@input[[direction]])
	saveRDS(drp_direction_out, file = snakemake@output[[direction]])
}


# for (j in seq_len(n_partitions)) {
# 	drp_partitions_F[[j]] <- as(
# 		list(uniques = drpF$uniques, quals = drpF$quals, map = maps_F[[j]]),
# 		"derep"
# 	)
# 	drp_partitions_R[[j]] <- as(
# 		list(uniques = drpR$uniques, quals = drpR$quals, map = maps_R[[j]]),
# 		"derep"
# 	)
# }


#
# for (direction in c("fwd", "rev")) {
# 	message(paste("Dereplicating",direction,"reads..."))
# 	# save each file as a 1-item list, with the name given as the sample filename
# 	x = list(
# 			dada2::derepFastq(
# 				snakemake@input[[direction]], verbose=TRUE)
# 		)
# 	names(x) <- basename(snakemake@input[[direction]])
# 	saveRDS(x, file = snakemake@output[[direction]])
# }

snakemake_log_end()

# for(j in seq(length(filtsF))) {
# 	drpF <- dada2::derepFastq(filtsF[[j]])
# 	saveRDS(drpF,derepsF[[j]])
# }
#
# for(j in seq(length(filtsR))) {
# 	drpR <- dada2::derepFastq(filtsR[[j]])
# 	saveRDS(drpR,derepsR[[j]])
# }
