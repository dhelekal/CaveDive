#' @export
simulate_outbreaks <- function(poi_rate, concentration, sampling_times, r_mean, r_sd, K_mean, K_sd, time_rate, time_shape) {
    
    n_tips <- length(sampling_times)

    sam <- sampling_times - max(sampling_times)
    sam <- sam[order(-sam)]

    n_exp <- rpois(1, poi_rate)+1
    expansion_probs <- rdirichlet(1, rep(concentration, n_exp))

    colouring <- rmultinom(n_tips, 1, expansion_probs)
    colouring <- sapply(c(1:n_tips), function (i) which(colouring[,i]>0))

    N <- rlnorm(1, meanlog = K_mean, sdlog = K_sd)
    K <- rlnorm(n_exp-1, meanlog = K_mean, sdlog = K_sd)
    A <- rlnorm(n_exp-1, meanlog = r_mean, sdlog = r_sd)

    div_times <- min(sam)-rgamma(n_exp-1, shape=time_shape, scale = 1/time_rate)
    div_times <- div_times[order(-div_times)]
    div_times <- c(div_times, -Inf)

    div_cols <- c(1:(n_exp))

    rates <- lapply(c(1:(n_exp-1)), function(i) return(function (s) sat.rate(s, K[i], A[i], div_times[i])))
    rates[[n_exp]] <- function (s) constant.rate(s, N)

    rate.ints <- lapply(c(1:(n_exp-1)), function(i) return(function (t, s) sat.rate.int(t, s, K[i], A[i], div_times[i])))
    rate.ints[[n_exp]] <- function(t, s) constant.rate.int(t, s, N)

    co <- structured_coal.simulate(sam, colouring, div_times, div_cols, rates, rate.ints)

    full_lh <- co$log_lh + 
               log(ddirichlet(expansion_probs, rep(concentration, n_exp))) +
               sum(sapply((div_cols), function(i) length(which(colouring==i)) * log(expansion_probs[i])))

    return(list(co=co, n_exp=n_exp, N=N, K=K, A=A, colours=colouring, div_times=div_times, div_cols=div_cols, exp_probs=expansion_probs, rates=rates, rate.ints=rate.ints, full_lh=full_lh))
}