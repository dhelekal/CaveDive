library("ape")

inv_t_conditional <- function(u,rate,delta_t){
        return((-1/rate)*log(1-u*(1-exp(-rate*delta_t))))
}

cond_likelihood <- function(rate, u, delta_t){
        return(rate*(1/(1-exp(-rate*delta_t)) - u))
}

exp_likelihood <- function(rate, t){
        return(rate*exp(-rate*t))
}

poi_0_likelihood <- function(rate, t){
        return(exp(-rate*t))
}

simulate <- function(sampling_times, pop_size){
        
        log_lh <- 0
                
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
                        rate <- choose(extant_lineages,2)/pop_size
                        
                        
                        if (idx+1 > length(sampling_times)) {
                                delta_t <- Inf
                        } else {
                                delta_t <- t-times_desc[idx+1]
                        }
                        
                        p_coal <- 1-exp(-rate*delta_t)
                        r <- runif(1,0,1)
                        
                        if  (r <= p_coal) {
                                r_w <- runif(1,0,1)
                                w_t <- inv_t_conditional(r_w,rate, delta_t)
                                
                                log_lh <- log_lh + log(p_coal) + log(cond_likelihood(rate, r_w, delta_t))
                                
                                t <- t - w_t
                                extant_lineages <- extant_lineages-1 
                                coalescent_times[coal_idx] <- t
                                coal_idx <- coal_idx+1
                        } else {
                                log_lh <- log_lh + log((1-p_coal))
                                idx <- idx+1
                                extant_lineages <- extant_lineages+1
                                future_lineages <- future_lineages-1
                                t <- times_desc[idx]
                        }
                }
        }
        
        ### add return likelihood
        coal_sample<-list(coalescent_times=coalescent_times,
                          log_likelihood=log_lh)

        return(coal_sample)
}

build_coal_tree <- function(sampling_times, coalescent_times){
        coal_times_desc <- coalescent_times[order(-coalescent_times)]
        times_desc <- sampling_times[order(-sampling_times)]
        tree_nodes <- seq(1,length(sampling_times))
        tree_nodes <- sapply(tree_nodes, function (x) return (paste0("S", x)))
        extant_entries <- c(1)
        extant_times <- c(times_desc[1])

        coal_idx <- 1 
        s_idx <- 1
        t <- times_desc[1]
        
        while (coal_idx <= length(coalescent_times)) {
                if (s_idx<length(times_desc) && (times_desc[s_idx+1] > coal_times_desc[coal_idx])){
                        s_idx <- s_idx+1  
                        t <- times_desc[s_idx]
                        extant_entries <- c(extant_entries, s_idx)
                        extant_times <- c(extant_times, t)
                } else {
                        t <- coal_times_desc[coal_idx]
                        coal_node1_idx <- trunc(runif(1,1,length(extant_entries)+1))
                        
                        coal_node1 <- extant_entries[coal_node1_idx]
                        ct1 <- extant_times[coal_node1_idx]
                        
                        extant_times <- extant_times[-coal_node1_idx]
                        extant_entries <- extant_entries[-coal_node1_idx]
                        
                        br_len_1 <- ct1-t
                        entry1 <- tree_nodes[coal_node1]
                        
                        coal_node2_idx <- trunc(runif(1,1,length(extant_entries)+1))
                        coal_node2 <- extant_entries[coal_node2_idx]
                        ct2 <- extant_times[coal_node2_idx]
                        br_len_2 <- ct2-t
                        entry2<-tree_nodes[coal_node2]
                        
                        tree_nodes[coal_node2] <- paste0("(",
                                                         entry1,
                                                         ":",br_len_1,
                                                         ",",
                                                         entry2,
                                                         ":",br_len_2,
                                                         ")")
                        extant_times[coal_node2_idx] <- t
                        coal_idx <- coal_idx+1
                }
        }
        tree_str<-paste0(tree_nodes[extant_entries[1]], ";")
        return(tree_str)
}

tree_likelihood <- function(sampling_times, coalescent_times, Ne){
        coal_times_desc <- coalescent_times[order(-coalescent_times)]
        times_desc <- sampling_times[order(-sampling_times)]
        
        n_sample <- length(times_desc)
        n_coal <- length(coal_times_desc)
        
        assertthat::are_equal(n_sample-1, n_coal)
        
        coal_idx <- 1
        sample_idx <- 1
        
        t <- times_desc[sample_idx]
        
        log_lh <- 0
        j <- 1
        
        while (coal_idx <=  n_coal) {
             if (sample_idx < n_sample &&
                 coal_times_desc[coal_idx] < times_desc[sample_idx+1]){
                     sample_idx <- sample_idx+1
                     delta_t <- t - times_desc[sample_idx]
                     rate <- choose(j,2)/Ne
                     log_lh <- log_lh+log(poi_0_likelihood(rate, delta_t))
                     
                     j <- j+1
                     t <- times_desc[sample_idx]
             } else {
                     delta_t <- t - coal_times_desc[coal_idx]
                     rate <- choose(j,2)/Ne
                     log_lh <- log_lh+log(exp_likelihood(rate, delta_t))
                     
                     j <- j-1
                     t <- coal_times_desc[coal_idx]
                     
                     coal_idx <-coal_idx+1
             }
        }
        return(log_lh)
}