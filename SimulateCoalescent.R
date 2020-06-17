library("ape")

inv_t_conditional <- function(u,rate,delta_t){
        return((-1/rate)*log(1-u*(1-exp(-rate*delta_t))))
}

cond_likelihood <- function(rate, u, delta_t){
        return(rate*(1/(1-exp(-rate*delta_t)) - u))
}

simulate <- function(sampling_times, pop_size){
        
        lh <- 1 
                
        times_desc <- sampling_times[order(-sampling_times)]
        future_lineages <- length(sampling_times)-1
        extant_lineages <- 1 
        coalescent_times <- rep(0, future_lineages)
        t <- 0 
        idx <- 1
        coal_idx <- 1
        
        while (future_lineages > 0 || extant_lineages > 1 ) {
                if (extant_lineages < 2) {
                        #If one lineage continue to next sampling event
                        idx <- idx+1
                        extant_lineages <- extant_lineages+1
                        future_lineages <- future_lineages-1
                        t <- times_desc[idx]
                } else {
                        #Pick waiting time
                        if (idx+1 > length(sampling_times)) {
                                delta_t <- Inf
                        } else {
                                delta_t <- t-times_desc[idx+1]
                        }
                        rate <- choose(extant_lineages,2)/pop_size*extant_lineages
                        w_t <- (-1.0/rate)*log(1-runif(n=1, min = 0, max= 1))
                        print(delta_t)
                        if (w_t>delta_t){
                                #if waiting time longer than interval
                                #move to next sampling event, restart
                                idx <- idx+1
                                extant_lineages <- extant_lineages+1
                                future_lineages <- future_lineages-1
                                t <- times_desc[idx]
                        } else {
                                #if waiting time within interval, coalesce
                                t <- t + w_t
                                extant_lineages <- extant_lineages-1 
                                coalescent_times[coal_idx] <- t
                                coal_idx <- coal_idx+1
                        }
                }
        }
        ### add return likelihood
        return(coalescent_times)
}

build_coal_tree <- function(sampling_times, coalescent_times){
        
}