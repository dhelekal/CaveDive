#' @export
clonal_tree_process.simulate_tree <- function(n_exp, N, K, A, sampling_times, tip_colours, div_times, div_cols) {

    rates <- if(n_exp > 1) lapply(c(1:(n_exp-1)), function(i) return(function (s) sat.rate(s, K[i], A[i], div_times[i]))) else list()
    rates[[n_exp]] <- function (s) constant.rate(s, N)

    rate.ints <- if(n_exp > 1) lapply(c(1:(n_exp-1)), function(i) return(function (t, s) sat.rate.int(t, s, K[i], A[i], div_times[i]))) else list()
    rate.ints[[n_exp]] <- function(t, s) constant.rate.int(t, s, N)

    co <- structured_coal.simulate(sampling_times, tip_colours, div_times, div_cols, rates, rate.ints)

    full_lh <- co$log_lh

    return(list(co=co, full_lh=full_lh))
}