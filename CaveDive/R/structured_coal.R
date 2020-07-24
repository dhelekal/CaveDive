structured_coal.simulate <- function(sampling_times, colours, div_times, div_events, Neg.rates,  Neg.rate.ints, Neg.rate.int_invs) {

    times_desc <- sampling_times[order(-sampling_times)]
    div_times_desc <- div_times[order(-div_times)]

    n_col <- length(unique(colours))

    extant_lineages <- rep(0, n_col)

    future_lineages <- sapply(unique(colours), function (x) sum(colours==x))
    coalescent_times <- rep(0, future_lineages)
    coalescent_cols <- rep(0, future_lineages)

    t <- 0
    t0 <- times_desc[1]

    div_idx <- 1
    coal_idx <- 1
    sam_idx <- 2

    log_lh <- 0

    while (max(future_lineages) > 0 || sum(extant_lineages) > 1) {
        if (max(extant_lineages) < 2) {
            ### If no lineage has more than 1 member continue to next divergence or sampling event
            if(sam_idx <= length(times_desc) && (div_idx > length(div_times_desc) || div_times_desc[div_idx] < times_desc[sam_idx])) {
                t <- t0 - times_desc[sam_idx]

                extant_lineages[colours[sam_idx]] <- extant_lineages[colours[sam_idx]] + 1    
                future_lineages[colours[sam_idx]] <- future_lineages[colours[sam_idx]] - 1
                sam_idx <- sam_idx + 1
            } else {
                which_div <- div_events[div_idx]
                if (extant_lineages[which_div] != 1) {
                    warning("Divergence event registered at this time but lineage MRCA has not been reached. This is likely a bug")
                    return(-1)
                }
                ### choose which lineage to diverge from
                extant_lineages[which_div] <- extant_lineages[which_div] - 1 
                i <- choose_reaction(extant_lineages)

                n_nodes <- sum(extant_lineages)
                extant_lineages[i] <- extant_lineages[i] + 1    
                div_idx <- div_idx + 1     

                log_lh <- log(extant_lineages[i])-log(n_nodes) 
            }
        } else {
            comb_ns <- sapply(extant_lineages, function (x) choose(x, 2))

            if (sam_idx > length(times_desc) && div_idx > length(div_times_desc)) {
                s <- Inf
            } else if(sam_idx > length(times_desc)){
                s <- t0 - t - div_times_desc[div_idx]
            } else if(div_idx > length(div_times_desc)){
                s <- t0 - t - times_desc[sam_idx]
            } else {
                s <- t0 - t - min(times_desc[sam_idx], div_times_desc[div_idx])
            }
            ### Compute rates for each lineage
            rates <- sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rate.ints[[x]](t, s))
            ### Decide if an coalescent event happens in this interval
            p_coal <- 1-exp(-sum(rates))

            r < -runif(1,0,1)

            if (r<=p_coal) { 
                ### choose which event happens and compute waiting time.
                i <- choose_reaction(rates) 
                log_lh + log(rates[i])-log(sum(rates))

                w_t <- inv_t_inhomogenous_exp_conditional(function(t, s)
                                                            comb_ns[i] * Neg.rate.ints[[i]](t, s),
                                                            Neg.rate.int_inv[[i]],
                                                            comb_ns[i],
                                                            t,
                                                            s)
                log_lh <- log_lh - log(comb_ns[i]) + inhomogenous_exp.loglh(function(s)
                                                            comb_ns[[i]] * Neg.rates[[i]](s),
                                                            function(t, s)
                                                            comb_ns[[i]] * Neg.rate.ints[[i]](t, s), t, w_t)
                                                        
        
                t <- t + w_t
                extant_lineages[i] <- extant_lineages[i] - 1
                coalescent_times[coal_idx] <- t0 - t
                coalescent_cols[coal_idx] <- i
                coal_idx <- coal_idx + 1
            } else {
                log_lh <- log_lh + log(1 - p_coal)
                ### determine whether sampling or divergence event comes next
                ### If no lineage has more than 1 member continue to next divergence or sampling event
                if(sam_idx <= length(times_desc) && (div_idx > length(div_times_desc) || div_times_desc[div_idx] < times_desc[sam_idx])) {
                    t <- t0 - times_desc[sam_idx]

                    extant_lineages[colours[sam_idx]] <- extant_lineages[colours[sam_idx]] + 1    
                    future_lineages[colours[sam_idx]] <- future_lineages[colours[sam_idx]] - 1
                    sam_idx <- sam_idx + 1
                } else {
                    which_div <- div_events[div_idx]
                    if (extant_lineages[which_div] != 1) {
                        warning("Divergence event registered at this time but lineage MRCA has not been reached. This is likely a bug")
                        return(-1)
                    }
                    ### choose which lineage to diverge from
                    extant_lineages[which_div] <- extant_lineages[which_div] - 1 
                    i <- choose_reaction(extant_lineages)

                    n_nodes <- sum(extant_lineages)
                    extant_lineages[i] <- extant_lineages[i] + 1    
                    div_idx <- div_idx + 1     

                    log_lh <- log_lh + log(extant_lineages[i])-log(n_nodes) 
                } 
            }
        }
    }
}

choose_reaction(rates) {
    total_rate <- sum(rates)
    r <- runif(1, 0, rates)

    i <- 0 
    rate_sum <- total_rate
    while (rate_sum > r) {
        i <- i+1
        rate_sum <- rate_sum-rates[i]
    }
    return(i)
}