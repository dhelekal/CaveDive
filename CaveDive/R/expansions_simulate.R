#' A phylogeny containing clonal expansions.
#' Note: The resulting phylogeny contains single child nodes denoting divergence events. use `collapse.singles` to remove these.
#' @param priors List of priors to simulate parameters from. See `standard priors` for details 
#' @param sampling_times Vector sampling times, entries must be negative or 0 
#' @param concentration Scalar concentration hyperparameter for dirichlet expansion model
#' @param given A list of variables with values given. Supported names: 'n_exp' - number of expansions, 'tip_colours' - tip expansion assignment, N' - Background population size, 'K' - Expansion carrying capacities, 't_mid' - Expansion times to midpoints, 'div_times' - Expansion divergence times, 'div_from' - parent lineages 
#' @param collapse_singles Whether to remove (collapse) divergence event nodes. Default: FALSE 
#' @return list of: `tree` - The simulated genealogy, `params` - the simulated parameters for the process
#' @export
simulate_expansion_phylo <- function(priors, sampling_times, concentration=2, given=list(), collapse_singles=FALSE) {
    sim <- expansions_simulate(priors, sampling_times, concentration=concentration, given=given)
    co <- sim$co
    params <- sim$params
    phy.div_nodes <- build_coal_tree.structured(sampling_times, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
    tree  <- read.tree(text = phy.div_nodes$full)
    if (collapse_singles){
        tree <- collapse.singles(tree)
    }
    return(list(tree=tree, params=params))
}


#' Simulate parameters and event times for a genealogy with clonal expansionss
#' @param priors List of priors to simulate parameters from. See `standard priors` for details 
#' @param sampling_times Vector sampling times, entries must be negative or 0 
#' @param concentration Scalar concentration hyperparameter for dirichlet expansion model
#' @param given A list of variables with values given. Supported names: 'n_exp' - number of expansions, 'tip_colours' - tip expansion assignment, N' - Background population size, 'K' - Expansion carrying capacities, 't_mid' - Expansion times to midpoints, 'div_times' - Expansion divergence times, 'div_from' - parent lineages  
#' @return list of: `co` - realisation of the prpcess, `params` - the simulated parameters for the process, `coal_log_lh` - the process likelihood, `param_log_lh` - the prior likelihood
#' @export
expansions_simulate <- function(priors, sampling_times, concentration, given=list()) {

    n_tips <- length(sampling_times)

    if(abs(max(sampling_times)) > 1e-6) {
        stop("maximum value of sampling_times must be 0")
    }

    if (is.null(given$n_exp)){
        n_exp <- priors$prior_i.sample()
    } else {
        n_exp <- given$n_exp
    }

    exp_probs <- rdirichlet(1, rep(concentration, (n_exp+1)))

    div_cols <- c(1:(n_exp+1))
    max_it <- 100
    
    if(is.null(given$tip_colours)){
        colouring <- c()
        it <- 0
        while((length(unique(colouring)) < (n_exp+1)) || (!all(clade_sizes>1) && n_exp > 0)) {
            if (it > max_it) {
                stop("Maximum sampling iterations exceeded. Invalid concentration or number of expansions. Unable to sample valid colouring")
            }
            colouring <- rmultinom(n_tips, 1, exp_probs)
            colouring <- sapply(c(1:n_tips), function (i) which(colouring[,i]>0))
            clade_sizes <- sapply(c(1:(n_exp+1)), function (i) length(which(colouring==i)))
            it <- it + 1
        }
    } else {
        colouring <- given$tip_colours 
    }

    exp_sizes <- sapply(div_cols, function(i) length(which(colouring==i)))
    if (is.null(given$N)) {
        N <- priors$prior_N.sample()
    } else {
        N <- given$N
    }

    if (is.null(given$K)) {
        K <- sapply(rep(N, n_exp), priors$prior_K_given_N.sample)
    } else {
        K <- given$K
    }

    if (is.null(given$t_mid)) {
        t_mid <- sapply(rep(N, n_exp), priors$prior_t_mid_given_N.sample)
    } else {
        t_mid <- given$t_mid
    }

    if(is.null(given$div_times)) {
        it <- 0
        div_times <- rep(Inf, (n_exp+1))
        while(!all(sapply(c(1:(n_exp+1)), function(i) all(div_times[i] < sampling_times[which(colouring==i)])))){
            if (it > max_it) {
                stop("Maximum sampling iterations exceeded.")
            }
            ## div time now
            div_times <- as.double(sapply(rep(N, n_exp), priors$prior_t_given_N.sample))
            div_times <- c(div_times,-Inf)

            ## Re-order everything so that divergence event numbering corresponds to their order of occurence
            div_times.ord <- order(-div_times) 
            div_times <- div_times[div_times.ord]

            colouring <- sapply(colouring, function(i) which(div_times.ord==i))
            it <- it+1
            exp_probs <- exp_probs[div_times.ord]
            exp_sizes <- exp_sizes[div_times.ord]
        }
    } else {
        div_times <- given$div_times
    }

    A <- sapply(t_mid, function(x) (1/x)**2)

    if (is.null(given$div_from)){
        if (n_exp > 0) {
            div_from <- sapply(c(1:n_exp), function(i) c((i+1):(n_exp+1))[runif(1,1,n_exp+2-i)])
        } else {
            div_from <- c()
        }
    } else {
        div_from <- given$div_from
    }

    clonal_co <- simulate_clonal_tree(n_exp, N, K, A, sampling_times, colouring, div_times, div_cols, div_from=div_from) 
    params <- list(n_exp=n_exp,
                    N=N, 
                    K=K, 
                    t_mid=t_mid, 
                    tip_colours=colouring, 
                    div_times=div_times, 
                    div_cols=div_cols, 
                    exp_probs=exp_probs)

    param_log_lh <- priors$prior_i(n_exp) +  
    priors$prior_N(N) +
    ddirichlet(t(matrix(exp_probs)), alpha=rep(concentration, length(exp_probs)), log=TRUE) +
    sum(sapply(c(1:length(exp_probs)), function (i) log(exp_probs[i])*exp_sizes[[i]])) -
    lgamma(n_exp+1) + ### Prior on parent populations
    lgamma(n_exp+1) ### Expansions are exchangeable

    if(n_exp > 0){
        param_log_lh <- param_log_lh +
        sum(priors$prior_K_given_N(K,N)) +
        sum(priors$prior_t_mid_given_N(t_mid,N)) +
        sum(priors$prior_t_given_N(div_times[-length(div_times)],N))  ## Do not compute likelihood for neutral population divergence time!!! 
    }

    return(list(co=clonal_co$co, params=params, coal_log_lh=clonal_co$log_lh, param_log_lh=param_log_lh))
}

#' Simulate an instance of coalescent process with local population structure 
#' @param n_exp number of expansions
#' @param N background population size
#' @param K carrying capacities
#' @param A growth rates
#' @param sampling_times sampling times
#' @param tip_colours tip clade assignment 
#' @param div_times clade divergence times
#' @param div_cols vector of colours ordered in sequence
#' @param div_from (optional) parent populations for individual clades. If NA will be randomised equiprobably.
#' @return list of: `co` - list containing event times, event colour assignment and parents of clonal expansions; log_lh: the log likelihood
#' @export
simulate_clonal_tree <- function(n_exp, N, K, A, sampling_times, tip_colours, div_times, div_cols, div_from=NA) {

    rates <- if(n_exp > 0) lapply(c(1:n_exp), function(i) return(function (s) sat.rate(s, K[i], A[i], div_times[i]))) else list()
    rates[[n_exp+1]] <- function (s) constant.rate(s, N)

    rate.ints <- if(n_exp > 0) lapply(c(1:n_exp), function(i) return(function (t, s) sat.rate.int(t, s, K[i], A[i], div_times[i]))) else list()
    rate.ints[[n_exp+1]] <- function(t, s) constant.rate.int(t, s, N)

    co <- structured_coal.simulate(sampling_times, tip_colours, div_times, div_cols, rates, rate.ints, div.from=div_from)

    return(list(co=co, log_lh=co$log_lh))
}