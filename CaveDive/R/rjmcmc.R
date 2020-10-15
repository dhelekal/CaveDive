rjmcmc <- function(likelihood, prior, proposal.cond_log_lh, proposal.sampler, x0, i0, max_it, thinning=1) {

    pb <- utils::txtProgressBar(min=0, max=max_it, style = 3)
    out.para <- vector('list',max_it/thinning)
    out.i <- vector('list',max_it/thinning)
    out.lh <- vector('list',max_it/thinning)
    out.prior <- vector('list',max_it/thinning)

    x_prev <- x0
    i_prev <- i0
    prior_prev <- prior(x0, i0)
    lh_prev <- likelihood(x0, i0)

    for(it in c(1:max_it)) {
        if (it%%thinning == 0) {
            utils::setTxtProgressBar(pb, it)
            out.para[[it/thinning]] <- x_prev 
            out.i[[it/thinning]] <- i_prev
            out.lh[[it/thinning]] <- lh_prev
            out.prior[[it/thinning]] <- prior_prev
        } 

        prop <- proposal.sampler(x_prev, i_prev)
        x_prop <- prop$x
        i_prop <- prop$i 

        prior_prop <- prior(x_prop, i_prop)

        if (prior_prop > -Inf) {
            lh_prop <- likelihood(x_prop, i_prop)
        } else {
            lh_prop <- -Inf
        }

        if (lh_prop > -Inf){
            a <- proposal.cond_log_lh(x_prev, i_prev, x_prop, i_prop)+lh_prop+prior_prop
            b <- proposal.cond_log_lh(x_prop, i_prop, x_prev, i_prev)+lh_prev+prior_prev
            alpha <- min(a-b, 0)

            r <- log(runif(1))
            if (r < alpha) {
                lh_prev <- lh_prop
                prior_prev <- prior_prop
                x_prev <- x_prop
                i_prev <- i_prop
            }
        }
    }
    return(list(para=out.para, dims=out.i, log_lh=out.lh, log_prior=out.prior))
}