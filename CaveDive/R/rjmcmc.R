rjmcmc <- function(likelihood, proposal.cond_log_lh, proposal.sampler, x0, i0, max_it, thinning=1) {

    pb <- utils::txtProgressBar(min=0, max=max_it, style = 3)
    out.para <- vector('list',max_it/thinning)
    out.i <- vector('list',max_it/thinning)

    x_prev <- x0
    i_prev <- i0
    prev_lh <- likelihood(x0, i0)

    for(it in c(1:max_it)) {
        if (it%%thinning == 0) {
            utils::setTxtProgressBar(pb, it)
            out.para[[it/thinning]] <- x_prev 
            out.i[[it/thinning]] <- i_prev
        } 

        prop <- proposal.sampler(x_prev, i_prev)
        x_prop <- prop$x
        i_prop <- prop$i 

        prop_lh <- likelihood(x_prop, i_prop)

        if (prop_lh > -Inf){
            a <- proposal.cond_log_lh(x_prev, i_prev, x_prop, i_prop)+prop_lh
            b <- proposal.cond_log_lh(x_prop, i_prop, x_prev, i_prev)+prev_lh
            alpha <- min(a-b, 0)

            r <- log(runif(1))
            if (alpha < r) {
                x_prev <- x_prop
                i_prev <- i_prop
            }
        }
    }
    return(list(para=out.para, dims=out.i))
}