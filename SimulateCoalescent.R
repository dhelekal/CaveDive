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
                        rate <- choose(extant_lineages,2)/pop_size*extant_lineages
                        
                        if (idx+1 > length(sampling_times)) {
                                delta_t <- Inf
                        } else {
                                delta_t <- t-times_desc[idx+1]
                        }
                        
                        p_coal <- 1-exp(-rate*delta_t)
                        r <- runif(1,0,1)
                        
                        if  (r < p_coal) {
                                r_w <- runif(1,0,1)
                                w_t <- inv_t_conditional(r_w,rate, delta_t)
                                
                                lh <- lh * p_coal * cond_likelihood(rate, r_w, delta_t)
                                
                                t <- t + w_t
                                extant_lineages <- extant_lineages-1 
                                coalescent_times[coal_idx] <- t
                                coal_idx <- coal_idx+1
                        } else {
                                lh <- lh * (1-p_coal)
                                idx <- idx+1
                                extant_lineages <- extant_lineages+1
                                future_lineages <- future_lineages-1
                                t <- times_desc[idx]
                        }
                }
        }
        ### add return likelihood
        coal_sample<-list(coalescent_times=coalescent_times, likelihood=lh)

        return(coal_sample)
}

build_coal_tree <- function(sampling_times, coalescent_times){
        
}