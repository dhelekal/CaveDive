#' @export
simulate_outbreaks <- function(poi_rate, concentration, sampling_times, r_mean, r_sd, K_mean, K_sd, time_rate) {
    
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

    div_times <- min(sam)+rexp(n_exp-1, time_rate)
    div_times <- div_times[order(-div_times)]
}