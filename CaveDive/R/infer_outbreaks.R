#' @export
infer_outbreaks <- function(phy, prior_i, prior_N, prior_N.sample, prior_r, prior_r.sample, prior_K, prior_K.sample, n_it=1e6, thinning=1, debug=FALSE) {
    i_0 <- 0
    N_0 <- prior_N.sample()

    x_0 <- list()
    x_0[[1]] <- N_0 

    pre <- structured_coal.preprocess_phylo(tree)

    o <- rjmcmc(function(x, i) log_lh(x,
                                      i, 
                                      prior_i, 
                                      prior_r, 
                                      prior_K, 
                                      prior_times,
                                      prior_branch, 
                                      pre, 
                                      debug),
                function(x, i, x_given, i_given) prop.cond_lh(x,
                                                              i,
                                                              x_given, 
                                                              i_given, 
                                                              function(x_init) para.log_lh(x_init, prior_r, prior_K),
                                                              pre),
                function(x_prev, i_prev) prop.sampler(x_prev,
                                                      i_prev, 
                                                      pre, 
                                                      function() para.initialiser(prior_r.sample,
                                                                                  prior_K.sample, 
                                                                                  pre)
                                                      ),
                x_0, i_0, n_it, thinning)
    return(o)
}

para.initialiser <- function(prior_r, prior_K, pre){

    x_next <- vector(mode = "list", length = 4)

    rates <- prior_r()
    K <- prior_K()

    edges <- pre$edges.df

    inner_branches <- edges$id[which(edges$node.child>pre$n_tips)] 
    edges_subs <- edges[inner_branches,]
    total_len <- sum(edges_subs$length)

    r <- runif(1, 0, total_len)
    i <- 0 
    len <- 0
    
    while (r > len) {
        i <- i+1
        len <- len + edges_subs$length[i] 
    }

    div.branch <- edges_subs$id[i]
    div.times <- runif(1, pre$nodes.df$times[edges$node.parent[div.branch]],
                          pre$nodes.df$times[edges$node.child[div.branch]])

    x_next[[1]] <- rates
    x_next[[2]] <- K
    x_next[[3]] <- div.times
    x_next[[4]] <- div.branch

    return(x_next)
}

para.log_lh <- function(x, prior_r, prior_K) {

    rates <- x[[1]]
    K <- x[[2]]

    out <- 0
    out <- out + prior_r(rates) + prior_K(K)

    return(out)
}

log_lh <- function(x, i, prior_i, prior_r, prior_K, prior_times, prior_branch, pre, exclude_lh=FALSE){


    n_tips <- pre$n_tips
    total_branch_len <- sum(pre$edges.df$length)
    root_MRCA <- pre$phy$node.label[pre$nodes.df$id[which.min(pre$nodes.df$times)]-n_tips]
    root_div <- -Inf

    N <- x[[1]]
    print(x)
    print(i)

    if (i > 0) {
        rates <- sapply(c(1:i), function(j) x[[j+1]][[1]])
        K <- sapply(c(1:i), function(j) x[[j+1]][[2]])
        div.branch <- sapply(c(1:i), function(j) x[[j+1]][[3]])
        div.times <- sapply(c(1:i), function(j) x[[j+1]][[4]])

    } else {
        rates <- c()
        K <- c()
        div.branch <- c()
        div.times <- c()
    }

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
            prior <- 0
            if (i > 0) {
                prior_br <- sum(log(pre$edges.df$length[div.branch]/total_branch_len))
                prior <- prior_r(rates) + prior_K(K) + prior_N(N) + prior_br
            } 
            if (exclude_lh) lh <- 0 else lh <- structured_coal.likelihood(pre, MRCAs, div.times, rates, K, N, type="Sat")$log_lh
            lh <- lh + prior
        } else {
            lh <- -Inf
        }
    } else {
        lh <- -Inf
    }
    return(lh)
}#

