message("Loading error structure...")
err <- readRDS(snakemake@input[["error"]])

# functions to enforce monotonicity of error function; necessary for binned
# quality scores from NovaSeq data. see https://github.com/benjjneb/dada2/issues/791
# estimated error should increase as Q-score decreases, but we prefer estimates
# from the higher Q-scores, hence `monotonic_up(..., direction='reverse')`.
# getErrors(err) is 16 x 40 matrix, rows are transitions (A2A, A2C, etc.),
# columns are Q-scores.
monotonic <- function(lst, compare=max, direction='forward') {
    new <- numeric(length(lst))
    names(new) <- names(lst)

    order <- seq_along(lst)
    if (direction == 'reverse') {
        order <- rev(order)
    }

    bar <- lst[[order[[1]]]]
    for (ii in order) {
        bar <- compare(bar, lst[[ii]])
        new[[ii]] <- bar
    }
    return(new)
}
monotonic_up <- function(lst,...) { return(monotonic(lst, compare=max,...)) }
monotonic_down <- function(lst,...) { return(monotonic(lst, compare=min,...)) }
enforce_monotonicity <- function(err) {
    return( t(apply(getErrors(err), 1, monotonic_up, direction='reverse')) )
}

message("Saving monotonic error structure...")
saveRDS(enforce_monotonicity(err),snakemake@output[["error"]])
