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
    i_0 <- 0
    N_0 <- prior_N.sample()

    x_0 <- list()
    x_0[[1]] <- N_0 
    x_0[[2]] <- c(1)

    pre <- structured_coal.preprocess_phylo(tree)
    edges <- pre$edges.df
    nodes <- pre$nodes.df
    n_tips <- pre$n_tips

    inner_branches <- edges$id[which(edges$node.child>pre$n_tips)] 
    edges_subs <- edges[inner_branches,]
    total_branch_len <- sum(edges_subs$length)

    prior_probs <- function(probs) {
        out <- log(ddirichlet(probs, alpha=rep(concentration, length(probs))))
        if(out == Inf) {
            print("err prob prior")
            print(probs)
        }
        return(out) 
    }


    prop_branch_time <- function(times, div.branch) log((times-nodes$times[edges$node.parent[div.branch]]) /
                                       pre$edges.df$length[div.branch]) + 
                                       log(pre$edges.df$length[div.branch]/total_branch_len)


    prop_branch_time.sample <- function() {
        out <- vector(mode = "list", length = 2)
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
                                         pre),
                function(x, i, x_given, i_given) prop.cond_lh(x,
                                                              i,
                                                              x_given, 
                                                              i_given, 
                                                              function(x_init) para.log_lh(x_init,
                                                                                           prior_r, 
                                                                                           prior_K, 
                                                                                           prop_branch_time),
                                                              pre),
                function(x_prev, i_prev) prop.sampler(x_prev,
                                                      i_prev, 
                                                      pre, 
                                                      function() para.initialiser(prior_r.sample,
                                                                                  prior_K.sample, 
                                                                                  prop_branch_time.sample),
                                                      prob.initialiser,
                                                      fn_log_J,
                                                      fn_log_J_inv
                                                      ),
                x_0, i_0, n_it, thinning)
    return(o)
}

prob.initialiser <- function(x_prev, i_prev) {
    new_exp_p <- 1/(i_prev+2)
    old_probs <- x_prev[[2]]
    old_probs <- old_probs*(1-new_exp_p)
    return(c(old_probs, new_exp_p))
}

fn_log_J <- function(i_prev, x_prev, x_next) {
    new_exp_p <- 1/(i_prev+2)
    return(log(1-new_exp_p)*(i_prev+1) + log(new_exp_p))
}

fn_log_J_inv <- function(i_prev, x_prev, x_next) {
    probs <- x_prev[[2]]
    old_p_inv <- 1/probs[(i_prev+1)]
    scale <- 1/sum(probs[-(i_prev+1)])
    return(log(scale)*(i_prev) + log(old_p_inv))
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

    out <- 0
    out <- out + prior_r(rates) + prior_K(K) + prior_branch_time(div.times, div.branch)

    return(out)
}

log_prior <- function(x, i, prior_i, prior_r, prior_K, prior_t, prior_probs, pre) {

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
      i >= 0 &&
      all(K > 0) &&
      all(rates > 0) &&
      all(N > 0) && 
      all(probs>0) &&
      abs(sum(probs)-1) < 1e-6 &&
      all(!is.na(div.branch)) &&
      all(div.times > pre$nodes.df$times[pre$edges.df$node.parent[div.branch]]) &&
      all(div.times < pre$nodes.df$times[pre$edges.df$node.child[div.branch]])) 
    {
        MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)

        if (all(!is.na(MRCAs))) {
            prior <- prior_i(i) + prior_probs(probs)
            if (i > 0) {
                prior <- prior + 
                         sum(prior_r(rates)) +
                         sum(prior_K(K)) + 
                         sum(prior_N(N)) + 
                         sum(prior_t(div.times))
            }  else {
                prior <- prior + sum(prior_N(N))
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
