para.initialiser <- function(pre, prior_r, prior_K){

    x_next <- vector(mode = "list", length = length(x_prev))

    rates <- prior_r()
    K <- prior_K()

    edges <- pre$edges.df
    total_len <- sum(edges$length)

    r <- runif(1, 0, total_len)
    i <- 0
    len <- 0
    
    while (r < len) {
        i <- i+1
        len <- len + edges$length[i] 
    }

    div.branch <- edges$id[i]
    div.times <- runif(1, pre$nodes.df$times[edges$node.parent[div.branch]], pre$nodes.df$times[edges$node.child[div.branch]])

    x_next[[1]] <- rates
    x_next[[2]] <- K
    x_next[[4]] <- div.times
    x_next[[5]] <- div.branch

    return(x_next)
}

para.log_lh <- function(x, prior_r, prior_K) {

    rates <- x[[1]]
    K <- x[[2]]

    out <- 0
    out <- out + prior_r(rates) + prior_K(K)

    return(out)
}

log_lh <- function(x, i, prior_i, prior_r, prior_k, prior_times, prior_branch, lh){

    N <- x[[1]]

    rates <- sapply(c(1:i), function(j) x[[j+1]][[1]])
    K <- sapply(c(1:i), function(j) x[[j+1]][[2]])
    div.branch <- sapply(c(1:i), function(j) x[[j+1]][[3]])
    div.times <- sapply(c(1:i), function(j) x[[j+1]][[4]])

    if (
      i >= 0 &&
      all(K > 0) &&
      all(rates > 0) &&
      all(N > 0) && 
      all(!is.na(div.branch))&&
      all(div.times > pre$nodes.df$times[pre$edges.df$node.parent[div.branch]]) &&
      all(div.times < pre$nodes.df$times[pre$edges.df$node.child[div.branch]])) 
    {
        MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)
        MRCAs <- c(MRCAs, root_MRCA)
        div.times <- c(div.times, root_div)
        if (all(!is.na(MRCAs))) {
            prior_br <- 0#sum(log(pre$edges.df$length[div.branch]/total_branch_len))
            lh <-structured_coal.likelihood(pre, MRCAs, div.times, rates, K, N, type="Sat")$log_lh
            prior <- prior_rates + prior_K + prior_N + prior_br
            lh <- lh + prior

        #if (pre$which_half[div.branch] == 2) {
        #  print(lh)
        #}
        } else {
            lh <- -Inf
        }
    } else {
        lh <- -Inf
    }
    return(lh)
}