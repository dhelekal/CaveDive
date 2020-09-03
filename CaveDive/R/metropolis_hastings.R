#'Run metropolis-hastings mcmc. 
#' 
#' @param model.lh a function for the conditional model log-likelihood.
#' @param proposal.cond_lh a function for the conditional proposal log-likelihood \code{x_{t+1}} given \code{x_t}.
#' @param x0 initial value for the chain
#' @param n_it number of iterations
#' @return A list of mcmc jumps
#' @export
run_mcmc <- function(model.log_lh, proposal.cond_log_lh, proposal.sampler, x0, n_it, verbose=FALSE){
	
	out <- vector(mode = "list", length = n_it)
	x_prev <- x0

	for (i in c(1:n_it)) {
		if(verbose && (i %% 10 == 0)) print(paste0("iteration: ",i))
		out[[i]] <- x_prev
		x_cand <- proposal.sampler(x_prev, i)
		a_prob <- log_metropolis_ratio(model.log_lh, proposal.cond_log_lh, x_cand, x_prev, i)

		u <- log(runif(1, 0, 1))
		if (u<a_prob) {
			x_prev <- x_cand
		}
	}

	return(out)
}

log_metropolis_ratio <- function(model.log_lh, proposal.cond_log_lh, x_cand, x_prev, it) {
	out <- -Inf
	cand_log_lh <- model.log_lh(x_cand)
	if (cand_log_lh > -Inf) {
		a <- proposal.cond_log_lh(x_prev, x_cand, it)+cand_log_lh
		b <- proposal.cond_log_lh(x_cand, x_prev, it)+model.log_lh(x_prev)
		out <- min(a-b, 0)
	}

	return(out)
}

