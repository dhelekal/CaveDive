rjmcmc <- function(likelihood, prior, proposal.sampler, x0, i0, max_it, thinning=1) {

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
        prop_Jacc <- prop$log_J 
        qr <- prop$qr

        prior_prop <- prior(x_prop, i_prop)

        if (prior_prop > -Inf) {
            lh_prop <- likelihood(x_prop, i_prop)
        } else {
            lh_prop <- -Inf
        }

        if (lh_prop > -Inf){
            a <- lh_prop+prior_prop
            b <- lh_prev+prior_prev

            mh <- qr + a - b

            alpha <- min(mh+prop_Jacc, 0)

            r <- log(runif(1))

            if(is.na(r<alpha)) print(paste0("qr: ", qr, " A: ", a, " B: ", b, " J: ", prop_Jacc))

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