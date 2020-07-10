#'Run metropolis-hastings mcmc. 
#' 
#' @param model.lh a function for the conditional model log-likelihood.
#' @param proposal.cond_lh a function for the conditional proposal log-likelihood \code{x_{t+1}} given \code{x_t}.
#' @param x0 initial value for the chain
#' @param n_it number of iterations
#' @return A list of mcmc jumps
#' @export
run_mcmc <- function(model.log_lh, proposal.cond_log_lh, proposal.sampler, x0, n_it){
	
	out <- vector(mode = "list", length = n_it)
	x_prev <- x0

	for (i in c(1:n_it)) {
		out[[i]] <- x_prev
		x_cand <- proposal.sampler(x_prev)
		a_prob <- log_metropolis_ratio(model.log_lh, proposal.cond_log_lh, x_cand, x_prev)

		u <- log(runif(1, 0, 1))
		if (u<a_prob) {
			x_prev <- x_cand
		}
	}

	return(out)
}

log_metropolis_ratio <- function(model.log_lh, proposal.cond_log_lh, x_cand, x_prev) {
	a <- proposal.cond_log_lh(x_prev, x_cand)+model.log_lh(x_cand)
	b <- proposal.cond_log_lh(x_cand, x_prev)+model.log_lh(x_prev)
	return(min(a-b, 0))
}

