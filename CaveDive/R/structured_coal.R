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

    coalescent_times <- rep(0, sum(future_lineages)-1)
    coalescent_cols <- rep(0, sum(future_lineages)-1)
    div_from <- rep(0, n_col-1)

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

                div_from[div_idx] <- i
 
                div_idx <- div_idx + 1     

                log_lh <- log(extant_lineages[i])-log(n_nodes) 
            }
        } else {
            comb_ns <- sapply(extant_lineages, function (x) choose(x, 2))
            print(paste0("Extant Lineages: ", extant_lineages))
            if (sam_idx > length(times_desc) && div_idx > (length(div_times)-1)) {
                s <- Inf
            } else if(sam_idx > length(times_desc)){
                s <- t0 - t - div_times[div_idx]
            } else if(div_idx > (length(div_times)-1)){
                s <- t0 - t - times_desc[sam_idx]
            } else {
                s <- t0 - t - max(times_desc[sam_idx], div_times[div_idx])
            }


            if (t > t+s) {
                print("next time step must be greater than current time")
                return(-1)
            }

            if (!all(comb_ns >= 0)){
                print("combination numbers must be positive")
                return(-1)
            }

            print(paste0("S: ", s, " T: ", t))
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
                    dfun <- function(s) sum(sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rates[[x]](t+s)))
                    return(inv_rates(func, Fs, dfun))
                }

                w_t <- inv_t_inhomogenous_exp_conditional(rate_sum_int(Neg.rate.ints, comb_ns, n_col),
                                                            function(t, s) inv_int(t, s),
                                                            1,
                                                            t,
                                                            s)

                i <- choose_reaction(rates)
                if (rates[i]!=Inf) { 
                    log_lh <- log_lh + log(rates[i])-log(sum(rates))
                }

                
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

                    div_from[div_idx] <- i

                    div_idx <- div_idx + 1     

                    log_lh <- log_lh + log(extant_lineages[i])-log(n_nodes) 
                } 
            }
        }
    }
    return(list(times=coalescent_times, colours=coalescent_cols, div_from=div_from, log_lh=log_lh))
}

structured_coal.preprocess_phylo <- function(phy){
    labs <- c(tree$node.label, tree$tip.label)
    nodes <- nodeid(phy, labs)
    is_tip <- c(rep(FALSE, length(tree$node.label)), rep(TRUE, length(tree$tip.label)))

    times <- node.depth.edgelength(phy)
    times <- times - max(times)
    times <- times[nodes]

    edges.parent<- phy$edge[1]
    edges.child <- phy$edge[2]
    edges.len <- phy$edge.length

    nodes.df <- data.frame(id=nodes, times=times, is_tip=is_tip)
    edges.df <- data.frame(parent=edges.parent, child=edges.child, length=edge.length)

    nodes.df <- nodes.df[order(nodes),]
    edges.df <- edges.df[order(parent),]

    clades.list <- lapply(nodes.df$nodes, function(x) extract.clade(phy, x))

    return(list(phy=phy, nodes.df = nodes.df, edges.df = edges.df, clades.list = clades.list))
}

structured_coal.likelihood <- function(phylo.preprocessed, div.MRCA.nodes, div.times,  Neg.rates, Neg.rate.ints){
    subtrees <- lapply(div.MRCA.nodes, function (x) phylo.preprocessed$clades.list[[x]]) 
    times.ord <- order(div.times)
    k_div <- length(div.times)

    lineage.trees <- lapply(c(1:k_div),
        function(x) drop.tip(subtrees[[times.ord[x]]], nodeid(subtrees[[times.ord[x]]], unlist(lapply(c(x:k_div),
            function (y) subtrees[[times.ord[y]]]$tip.label
            )))))

    
}

choose_reaction <- function(rates) {

    if (sum(rates==Inf) > 1) {
        warning("More than one entry equal to infinity. This should not happen")
    } else if (!all(rates != Inf)) {
        which_inf <- which(rates==Inf)
        return(which_inf[runif(1,1,length(which_inf+1))])
    } else {

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
}

inv_rates <- function(func, N, dfun=NULL){   
    print("Numerically inverting function")
    print(paste0("With N = ", N))
    f <- function(x) func(x) - N
    lo <- 0
    hi <- 1e0
    if (sign(f(lo)) == sign(f(hi))){

        sg <- sign(f(lo))

        min_hi <- 1
        while(sign(f(hi*(10^(min_hi)))) == sg && min_hi < 20){
            min_hi <- min_hi+1
        }

        if (sign(f(hi*(10^(min_hi)))) == sg) {
            warning(paste0("Cannot find initial bisection search interval. ", "min_lo: ", min_lo, " min_hi: ", min_hi))
        }


        hi <- hi*(10^(min_hi))
    }

    print("initiating BSS")
    b <- bisect(f, lo, hi, maxiter = 50)

    out <- b$root

    if (b$f.root > 1e-8) {
        print("initiating newtonRaphson")
        n <- newtonRaphson(f,out,dfun=dfun)
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