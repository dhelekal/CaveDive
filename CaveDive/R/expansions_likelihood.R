#' @export
expansion_likelihood <- function(phylo.preprocessed, div.MRCA.nodes, div.times, diverging.rates, diverging.sizes, neutral.size, expansion.probs) {
    structured.log_lh <- structured_coal.likelihood(phylo.preprocessed, div.MRCA.nodes, div.times, diverging.rates, diverging.sizes, neutral.size)
    partition_counts <- structured.log_lh$partition_counts
    partition.log_lh <- sum(sapply(c(1:length(expansion.probs)), function (i) log(expansion.probs[i])*partition_counts[[i]]))

    return(structured.log_lh$log_lh + partition.log_lh)
}