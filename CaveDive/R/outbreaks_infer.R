offset <- 2

#' @export
outbreaks_infer <- function(phy,
                            prior_i, 
                            prior_N, 
                            prior_N.sample, 
                            prior_tmid_given_N, 
                            prior_tmid_given_N.sample, 
                            prior_K_given_N, 
                            prior_K_given_N.sample, 
                            prior_t_given_N,
                            prior_t_given_N.sample,
                            concentration,
                            n_it=1e6, thinning=1, init=NULL, debug=FALSE) {

    if (debug) warning("Running in debug mode with only priors in use.")

    pre <- structured_coal.preprocess_phylo(phy)
    edges <- pre$edges.df
    nodes <- pre$nodes.df
    n_tips <- pre$n_tips

    inner_branches <- edges$id[which(edges$node.child>pre$n_tips)] 
    edges_subs <- edges[inner_branches,]
    total_branch_len <- sum(edges_subs$length)


    all_times <- extract_lineage_times(pre, pre$phy$node.label[pre$root_idx-pre$n_tips], -Inf)

    tree_height <- max(nodes$times) - min(nodes$times)

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

    if (is.null(init)) {
    
        N_0 <- optim(1, const_log_lh, lower=(1e-2)*tree_height, upper=(1e2)*tree_height, method="Brent", control = list(maxit = 2000000))$par
    
        i_0 <- 0 
        x_0 <- list()
        x_0[[1]] <- N_0 
        x_0[[2]] <- c(1)
    } else {
        x_0 <- init$x_0
        i_0 <- init$i_0
    }

    prior_probs <- function(probs) {

        out <- -Inf
        if (abs(sum(probs)-1) < 1e-8 && all(probs > 0)) {
            out <- ddirichlet(t(matrix(probs)), alpha=rep(concentration, length(probs)), log=TRUE)
        }
        return(out) 
    }

    o <- rjmcmc(function(x, i) log_posterior(x,
                                         i,
                                         prior_i, 
                                         prior_N,
                                         prior_tmid_given_N, 
                                         prior_K_given_N, 
                                         prior_t_given_N, 
                                         prior_probs,
                                         pre),
                function(x_prev, i_prev) prop.sampler(x_prev,
                                                      i_prev, 
                                                      pre, 
                                                      function(N) para.initialiser(N, 
                                                                                   prior_tmid_given_N.sample,
                                                                                   prior_K_given_N.sample, 
                                                                                   prior_t_given_N.sample,
                                                                                   pre),
                                                      function(x_init, N) para.log_lh(x_init,
                                                                                      N,
                                                                                      prior_tmid_given_N, 
                                                                                      prior_K_given_N, 
                                                                                      prior_t_given_N,
                                                                                      pre),
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

para.initialiser <- function(N, prior_tmid_given_N, prior_K_given_N, prior_t_given_N, pre){
    edges <- pre$edges.df
    nodes <- pre$nodes.df
    x_next <- vector(mode = "list", length = 4)

    mid.times <- prior_tmid_given_N(N)
    K <- prior_K_given_N(N)
    div.times <- prior_t_given_N(N)

    ### which branches exist at time of divergence
    br_extant_before <- edges$id[which(nodes$times[edges$node.child] > div.times)]
    br_extant_after <- br_extant_before[which(nodes$times[edges$node.parent[br_extant_before]] < div.times)]

    ### filter out terminal branches as those have 0 prior mass
    extant_inner <- br_extant_after[which(edges$node.child[br_extant_after]>pre$n_tips)] 
    ### choose one at random

    if(length(extant_inner) < 1) extant_inner <- c(NA)


    div.branch <- extant_inner[sample.int(length(extant_inner),1)] 

    x_next[[1]] <- mid.times
    x_next[[2]] <- K
    x_next[[3]] <- div.times
    x_next[[4]] <- div.branch

    return(x_next)
}

para.log_lh <- function(x, N, prior_tmid_given_N, prior_K_given_N, prior_t_given_N, pre) {
    edges <- pre$edges.df
    nodes <- pre$nodes.df

    mid.times <- x[[1]]
    K <- x[[2]]
    div.times <- x[[3]]
    div.branch <- x[[4]]

    ### which branches exist at time of divergence
    br_extant_before <- edges$id[which(nodes$times[edges$node.child] > div.times)]
    br_extant_after <- br_extant_before[which(nodes$times[edges$node.parent[br_extant_before]] < div.times)]

    ### filter out terminal branches as those have 0 prior mass
    extant_inner <- br_extant_after[which(edges$node.child[br_extant_after]>pre$n_tips)] 
    if (length(extant_inner) < 1) extant_inner <- c(NA)

    out <- prior_tmid_given_N(mid.times,N) + prior_K_given_N(K,N) + prior_t_given_N(div.times,N) + log(1/length(extant_inner))
    return(out)
}

log_posterior <- function(x,
                            i, 
                            prior_i, 
                            prior_N, 
                            prior_tmid_given_N, 
                            prior_K_given_N, 
                            prior_t_given_N, 
                            prior_probs, 
                            pre) 
{
    prior <- 0
    lh <- 0

    edges <- pre$edges.df
    nodes <- pre$nodes.df

    n_tips <- pre$n_tips
    root_MRCA <- pre$phy$node.label[nodes$id[which.min(nodes$times)]-n_tips]
    root_div <- -Inf

    ### Extract values

    N <- x[[1]]
    probs <- x[[2]]
    if (i > 0) {
        div.times <- sapply(c(1:i), function(j) x[[j+offset]][[3]])

        div_ord <- order(-div.times)
        div.times <- div.times[div_ord]

        mid.times <- sapply(div_ord, function(j) x[[j+offset]][[1]])
        K <- sapply(div_ord, function(j) x[[j+offset]][[2]])
        div.branch <- sapply(div_ord, function(j) x[[j+offset]][[4]])
    } else {
        mid.times <- c()
        K <- c()
        div.branch <- c()
        div.times <- c()
    }

    ### Check that all values make sense and lie within support of prior and likelihood functions
    if (
      (i >= 0) &&
      all(K > 0) &&
      all(mid.times > 0) &&
      all(N > 0) && 
      all(probs >= 0) &&
      (abs(sum(probs)-1) < 1e-8) &&
      all(!is.na(div.branch)) &&
      all(div.times > nodes$times[edges$node.parent[div.branch]]) && 
      all(div.times < nodes$times[edges$node.child[div.branch]])) 
    {
        MRCAs <- sapply(edges$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)
        if (all(!is.na(MRCAs))) { ### Make sure no terminal branch is being proposed as that is a zero set. 
            
            prior <- prior_i(i) + prior_probs(probs) + prior_N(N) + 
                 lgamma(i+1) ### Correction for exchangeable RVs
            if (i > 0) {
                    prior <- prior + 
                             sum(prior_tmid_given_N(mid.times, N)) +
                             sum(prior_K_given_N(K,N)) + 
                             sum(prior_t_given_N(div.times,N)) -
                             lgamma(length(div.times)) ### prior on divergence events
            }

            MRCAs_root <- c(MRCAs, root_MRCA) ### add root for parent population
            div.times_root <- c(div.times, root_div) ### parent diverges at -Inf

            ### structured coal likelihood accepts rates, need to transform mid point to a rate
            rates <- sapply(mid.times, function (x) (1/x)**2)

            structured.log_lh <- structured_coal.likelihood(pre,
                                                            MRCAs_root, 
                                                            div.times_root, 
                                                            rates, 
                                                            K, 
                                                            N)

            partition_counts <- structured.log_lh$partition_counts
            partition_prior <- sum(sapply(c(1:length(probs)), function (i) log(probs[i])*partition_counts[[i]]))

            prior <- prior + partition_prior
            lh <- structured.log_lh$log_lh
        } else {
            prior <- -Inf
            lh <- -Inf
        }     
    } else {
        prior <- -Inf
        lh <- -Inf
    }
    return(list(prior=prior, lh=lh))
}