offset <- 2

#' @export
outbreaks_infer <- function(phy,
                            prior_i, 
                            prior_N, 
                            prior_N.sample, 
                            prior_r, 
                            prior_r.sample, 
                            prior_K, 
                            prior_K.sample, 
                            prior_t,
                            prior_t.sample,
                            concentration,
                            n_it=1e6, thinning=1, debug=FALSE) {

    if (debug) warning("Running in debug mode with only priors in use.")

    pre <- structured_coal.preprocess_phylo(phy)
    edges <- pre$edges.df
    nodes <- pre$nodes.df
    n_tips <- pre$n_tips

    inner_branches <- edges$id[which(edges$node.child>pre$n_tips)] 
    edges_subs <- edges[inner_branches,]
    total_branch_len <- sum(edges_subs$length)


    all_times <- extract_lineage_times(pre, pre$phy$node.label[pre$root_idx-pre$n_tips], -Inf)

    tree_height <- max(pre$nodes.df$times) - min(pre$nodes.df$times)

    const_log_lh <- function(n){
        if (n > 0){
            lh <- -coalescent_loglh(all_times$sam.times[[1]],
                        all_times$coal.times[[1]],
                        n,
                        0)
        } else {
            lh <- Inf
        }
        return(lh)
    }

    N_0 <- optim(1, const_log_lh, lower=1e-3, upper=1e5, method="Brent", control = list(maxit = 2000000))$par

    i_0 <- 0 
    x_0 <- list()
    x_0[[1]] <- N_0 
    x_0[[2]] <- c(1)

    prior_probs <- function(probs) {

        out <- -Inf
        if (abs(sum(probs)-1) < 1e-8 && all(probs > 0)) {
            out <- ddirichlet(t(matrix(probs)), alpha=rep(concentration, length(probs)), log=TRUE)
        }
        return(out) 
    }

    prior_edges <- function(div.edges) { ### this is actually the partition probability
        return(length(div.edges)*log(1/length(inner_branches)))
    }
 
    prop_branch_time <- function(times, div.branch) {
        return(log(1/length(inner_branches)) + log(1/tree_height))
    }

    prop_branch_time.sample <- function() {
        edges <- pre$edges.df
        nodes <- pre$nodes.df

        div.branch <- sample(inner_branches, 1)
        div.times <- runif(1, min(nodes$times), max(nodes$times))

        return(list(time=div.times, branch=div.branch))
    }

    o <- rjmcmc(function(x, i) log_lh(x,
                                      i, 
                                      pre, 
                                      debug),
                function(x, i) log_prior(x,
                                         i,
                                         prior_i, 
                                         prior_r, 
                                         prior_K, 
                                         prior_t, 
                                         prior_probs,
                                         prior_edges,
                                         pre),
                function(x_prev, i_prev) prop.sampler(x_prev,
                                                      i_prev, 
                                                      pre, 
                                                      function() para.initialiser(prior_r.sample,
                                                                                  prior_K.sample, 
                                                                                  prop_branch_time.sample),
                                                      function(x_init) para.log_lh(x_init,
                                                                                           prior_r, 
                                                                                           prior_K, 
                                                                                           prop_branch_time),
                                                      fn_log_J,
                                                      fn_log_J_inv,
                                                      pop_scale=(tree_height/2) ### good enough approximation of the population size
                                                      ),
                x_0, i_0, n_it, thinning)
    return(o)
}

fn_log_J <- function(i_prev, x_prev, x_next) {
    return(0)
}

fn_log_J_inv <- function(i_prev, x_prev, x_next, which_mdl_rm) {
    return(0)
}

para.initialiser <- function(prior_r, prior_K, prior_branch_time, pre){

    x_next <- vector(mode = "list", length = 4)

    rates <- prior_r()
    K <- prior_K()
    br_time <- prior_branch_time()

    x_next[[1]] <- rates
    x_next[[2]] <- K
    x_next[[3]] <- br_time$time
    x_next[[4]] <- br_time$branch

    return(x_next)
}

para.log_lh <- function(x, prior_r, prior_K, prior_branch_time) {

    rates <- x[[1]]
    K <- x[[2]]
    div.times <- x[[3]]
    div.branch <- x[[4]]

    out <- prior_branch_time(div.times, div.branch) + prior_r(rates) + prior_K(K)
    return(out)
}

log_prior <- function(x, i, prior_i, prior_r, prior_K, prior_t, prior_probs, prior_edges, pre) {

    n_tips <- pre$n_tips
    N <- x[[1]]
    probs <- x[[2]]

    if (i > 0) {
        rates <- sapply(c(1:i), function(j) x[[j+offset]][[1]])
        K <- sapply(c(1:i), function(j) x[[j+offset]][[2]])
        div.times <- sapply(c(1:i), function(j) x[[j+offset]][[3]])
        div.branch <- sapply(c(1:i), function(j) x[[j+offset]][[4]])
    } else {
        rates <- c()
        K <- c()
        div.branch <- c()
        div.times <- c()
    }

    if (
      (i >= 0) &&
      all(K > 0) &&
      all(rates > 0) &&
      all(N > 0) && 
      all(probs >= 0) &&
      (abs(sum(probs)-1) < 1e-8) &&
      all(!is.na(div.branch)) &&
      all(div.times > pre$nodes.df$times[pre$edges.df$node.parent[div.branch]]) && ## technically last two inequalities part of likelihood
      all(div.times < pre$nodes.df$times[pre$edges.df$node.child[div.branch]]))    ## but they assign 0 likelihood and easier to check here

    {
        prior <- prior_i(i) + prior_probs(probs) + prior_N(N)
        MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)

        if (all(!is.na(MRCAs))){
            if (i > 0) {
                prior <- prior + 
                         sum(prior_edges(div.branch)) +
                         sum(prior_r(rates)) +
                         sum(prior_K(K)) + 
                         sum(prior_t(div.times)) + lgamma(i+1)
            }
        } else {
            prior <- -Inf
        }
    } else {
        prior <- -Inf 
    }
    if (prior==Inf) print("Error-prior")
    return(prior)
}

log_lh <- function(x, i, pre, exclude_lh=FALSE){

    if (exclude_lh)
    {
        lh <- 0
    } else {
        n_tips <- pre$n_tips

        root_MRCA <- pre$phy$node.label[pre$nodes.df$id[which.min(pre$nodes.df$times)]-n_tips]
        root_div <- -Inf

        N <- x[[1]]
        probs <- x[[2]]

        if (i > 0) {
            rates <- sapply(c(1:i), function(j) x[[j+offset]][[1]])
            K <- sapply(c(1:i), function(j) x[[j+offset]][[2]])
            div.times <- sapply(c(1:i), function(j) x[[j+offset]][[3]])
            div.branch <- sapply(c(1:i), function(j) x[[j+offset]][[4]])
        } else {
            rates <- c()
            K <- c()
            div.branch <- c()
            div.times <- c()
        }

        MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)
        MRCAs <- c(MRCAs, root_MRCA)
        div.times_root <- c(div.times, root_div)

        lh <- outbreaks_likelihood(pre, MRCAs, div.times_root, rates, K, N, probs)
    }

    return(lh)
}
