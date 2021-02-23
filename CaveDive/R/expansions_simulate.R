#' @export
expansions_simulate <- function(priors, sampling_times, concentration) {
    
    n_tips <- length(sampling_times)

    if(abs(max(sampling_times)) > 1e-6) {
        stop("maximum value of sampling_times must be 0")
    }

    n_exp <- priors$prior_i.sample()
    expansion_probs <- rdirichlet(1, rep(concentration, (n_exp+1)))
    div_cols <- c(1:(n_exp+1))
    colouring <- c()
    max_it <- 100
    it <- 0
    while((length(unique(colouring)) < (n_exp+1)) || (!all(clade_sizes>1) && n_exp > 0)) {
        if (it > max_it) {
            stop("Maximum sampling iterations exceeded. Invalid concentration or number of expansions. Unable to sample valid colouring")
        }
        exp_probs <- rdirichlet(1, concentration)
        colouring <- rmultinom(n_tips, 1, exp_probs)
        colouring <- sapply(c(1:n_tips), function (i) which(colouring[,i]>0))
        it <- it + 1
    }

    exp_sizes <- sapply(div_cols, function(i) length(which(colouring==i)))

    N <- priors$prior_N.sample()
    K <- sapply(rep(N, n_exp), priors$prior_K_given_N.sample)
    t_mid <- sapply(rep(N, n_exp), priors$prior_t_mid_given_N.sample)

    most_recent_sam <- sapply(div_cols, function (i) min(sampling_times[which(colouring==i)]))
    div_times <- rep(Inf, n_exp)

    it <- 0
    while (!all(div_times < most_recent_sam)) {
        if (it > max_it) {
            stop("Maximum sampling iterations exceeded. Invalid Sampling Times. Cannot simulate divergence times compatible with sampling times.")
        }
        div_times <- sapply(rep(N, n_exp), priors$prior_t_given_N.sample)
        div_times <- c(div_times, -Inf)
        it <- it + 1
    }

    A <- sapply(t_mid, function(x) (1/x)**2)

    co <- clonal_tree_process.simulate_tree(n_exp, N, K, A, sampling_times, colouring, div_times, div_cols) 
    params <- list(n_exp=n_exp, N=N, K=K, t_mid=t_mid, tip_colours=colouring, div_times=div_times, div_cols=div_cols, exp_probs=expansion_probs)
    param_log_lh <- priors$prior_i(n_exp) +  
                    priors$prior_N(N) +
                    sum(priors$prior_K_given_N(K,N)) +
                    sum(priors$prior_t_mid_given_N(t_mid,N)) + 
                    sum(priors$prior_t_given_N(div_times,N)) +
                    ddirichlet(t(matrix(exp_probs)), alpha=rep(concentration, length(exp_probs)), log=TRUE) +
                    sum(sapply(c(1:length(exp_probs)), function (i) log(exp_probs[i])*exp_sizes[[i]]))

    return(list(co=co, params=params, coal_log_lh=co$log_lh, param_log_lh=param_log_lh))
} 