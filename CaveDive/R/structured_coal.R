#' Simulate outbreak structured coalescent
#' Population size must be locally smooth around divergence time and must reach 0 at divergence time !!
#' Note: This simulation relies on inverse transform sampling and numerical inverses of Neg.rate.ints  
#'
#' @param sampling_times Nx1 times of leaves.
#' @param colours Nx1 colour assignment of leaves
#' @param div_times Kx1 divergence times for colours, neutral colour must diverge at -Inf. MUST BE IN DESCENDING ORDER
#' @param div_events Kx1 colour assignment for divergence times, i.e. which colour corresponds to which divergence time
#' @param Neg.rates Kx1 List of Functions (t)->(R+). 1/Neg(t) for each colour.
#' @param Neg.rate.ints Kx1 List of Functions (t, s)->(R+). \int_{t}^{s} 1/Neg(\tau)\,d\tau for each colour
#' @return A list consisting of the simulate coalescent times \code{coalescent_times} and the log-likelihood of the simulation \code{log_likelihood}.
#' @export
structured_coal.simulate <- function(sampling_times, colours, div_times, div_events, Neg.rates, Neg.rate.ints) {

    sam_ord <- order(-sampling_times)
    times_desc <- sampling_times[sam_ord]
    colours_desc <- colours[sam_ord]

    n_col <- length(unique(colours_desc))

    if (n_col != length(div_times)) {
        warning("Number of divergence events does not match number of unique colours")
        return(-1)
    }

    extant_lineages <- rep(0, n_col)

    future_lineages <- sapply(div_events, function (x) sum(colours_desc==x)) ###maybe wrong????
    future_lineages <- future_lineages[order(div_events)]

    coalescent_times <- rep(0, sum(future_lineages))
    coalescent_cols <- rep(0, sum(future_lineages))

    t <- 0
    t0 <- times_desc[1]

    div_idx <- 1
    coal_idx <- 1
    sam_idx <- 2

    extant_lineages[colours_desc[1]] <- 1
    future_lineages[colours_desc[1]] <- future_lineages[colours_desc[1]] - 1

    log_lh <- 0

    while (max(future_lineages) > 0 || sum(extant_lineages) > 1) {
        if (max(extant_lineages) < 2) {
            ### If no lineage has more than 1 member continue to next divergence or sampling event
            if(sam_idx <= length(times_desc) && (div_idx > length(div_times) || div_times[div_idx] < times_desc[sam_idx])) {
                t <- t0 - times_desc[sam_idx]

                which_colour <- colours_desc[sam_idx]

                extant_lineages[which_colour] <- extant_lineages[which_colour] + 1    
                future_lineages[which_colour] <- future_lineages[which_colour] - 1
                sam_idx <- sam_idx + 1

            } else {
                t <- t0 - div_times[div_idx]

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

            if (sam_idx > length(times_desc) && div_idx > (length(div_times)-1)) {
                s <- Inf
            } else if(sam_idx > length(times_desc)){
                s <- t0 - t - div_times[div_idx]
            } else if(div_idx > (length(div_times)-1)){
                s <- t0 - t - times_desc[sam_idx]
            } else {
                s <- t0 - t - max(times_desc[sam_idx], div_times[div_idx])
            }

            print(paste0("S: ", s, "T: ", t))
            ### Compute rates for each lineage
            rates <- sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rate.ints[[x]](t, s))
            print(paste0("rates: ", rates))
            ### Decide if an coalescent event happens in this interval
            p_coal <- 1-exp(-sum(rates))
            r <- runif(1,0,1)

            if (r<=p_coal) { 
                ### choose which event happens and compute waiting time.

                inv_int <- function(t, Fs) {
                    func <- function(s) rate_sum_int(Neg.rate.ints, comb_ns, n_col)(t,s)
                    return(inv_rates(func, Fs))
                }

                w_t <- inv_t_inhomogenous_exp_conditional(function(t, s) inv_int(t, s),
                                                            rate_sum_int(Neg.rate.ints, comb_ns, n_col),
                                                            sum(comb_ns),
                                                            t,
                                                            s)

                i <- choose_reaction(rates) 
                log_lh <- log_lh + log(rates[i])-log(sum(rates))

                
                log_lh <- log_lh - log(comb_ns[i]) + inhomogenous_exp.loglh(
                                                            function(s) sum(sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rates[[x]](s))),
                                                            rate_sum_int(Neg.rate.ints, comb_ns, n_col), t, w_t)
                                                        
        
                t <- t + w_t
                extant_lineages[i] <- extant_lineages[i] - 1
                coalescent_times[coal_idx] <- t0 - t
                coalescent_cols[coal_idx] <- i
                coal_idx <- coal_idx + 1
            } else {
                log_lh <- log_lh + log(1 - p_coal)
                ### determine whether sampling or divergence event comes next
                ### If no lineage has more than 1 member continue to next divergence or sampling event
                if(sam_idx <= length(times_desc) && (div_idx > (length(div_times)-1) || div_times[div_idx] < times_desc[sam_idx])) {
                    t <- t0 - times_desc[sam_idx]

                    which_colour <- colours_desc[sam_idx]

                    extant_lineages[which_colour] <- extant_lineages[which_colour] + 1    
                    future_lineages[which_colour] <- future_lineages[which_colour] - 1
                    sam_idx <- sam_idx + 1
                } else {
                    t <- t0 - div_times[div_idx]

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
    return(list(times=coalescent_times, colours=coalescent_cols, log_lh=log_lh))
}

choose_reaction <- function(rates) {
    total_rate <- sum(rates)
    r <- runif(1, 0, total_rate)

    i <- 0 
    rate_sum <- total_rate
    while (rate_sum > r) {
        i <- i+1
        rate_sum <- rate_sum-rates[i]
    }
    return(i)
}

inv_rates <- function(func, N){    
    f <- function(x) func(x) - N
    b <- bisect(f, 1e-6, 1e6, maxiter = 20)
    out <- b$root
    if (b$f.root > 1e-8) {
        n <- newtonRaphson(f,out)
        if (n$f.root >  1e-8) {
            warning(paste0("Function root suspiciously large, please revise. f.root: ", n$f.root))
        }
        out <-n$root
    }
    return(out)
}

rate_sum_int <- function(rate_ints, comb_ns, n_col) {
    return(function(t,s) sum(sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*rate_ints[[x]](t, s))))
}