function <- mcmc(model.lh, proposal.cond_lh, proposal.sampler, x0, n_it){
	
	out <- rep(0,  n_it)
	x_prev <- x_0

	for (i in c(1:n_it)) {
		out[1] <- x_prev
		x_cand <- proposal.sampler(x_prev)
		a_prob <- metropolis_ratio(model.lh, proposal.cond_lh, x_cand, x_prev)

		u <- runif(1, 0, 1)
		if (u<a_prob) {
			x_prev <- x_cand
		}
	}

	return(out)
}

function <- metropolis_ratio(model.lh, proposal.cond_lh, x_cand, x_prev) {
	a <- proposal.cond_lh(x_prev, x_cand)*model.lh(x_cand)
	b <- proposal.cond_lh(x_cand, x_prev)*model.lh(x_prev)
	return(min(a/b ,1))
}

