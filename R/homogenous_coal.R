#'Simulate homogenous coalescent,
#' 
#' @param sampling_times times of leaves.
#' @param pop_size Effective population size \code{Neg}.
#' @return A list consisting of the simulate coalescent times \code{coalescent_times} and the log-likelihood of the simulation \code{log_likelihood}.
#' @export
homogenous_coal.simulate <- function(sampling_times, pop_size) {
        log_lh <- 0
        
        times_desc <- sampling_times[order(-sampling_times)]
        future_lineages <- length(sampling_times) - 1
        extant_lineages <- 1
        coalescent_times <- rep(0, future_lineages)
        t <- 0
        t0 <- times_desc[1]
        idx <- 1
        coal_idx <- 1
        
        while (future_lineages > 0 || extant_lineages > 1) {
                if (extant_lineages < 2) {
                        #If one lineage continue to next sampling event
                        t <- t0 - times_desc[idx + 1]
                        idx <- idx + 1
                        extant_lineages <- extant_lineages + 1
                        future_lineages <- future_lineages - 1
                } else {
                        #Pick waiting time
                        rate <-
                                choose(extant_lineages, 2) / pop_size
                        
                        if (idx + 1 > length(sampling_times)) {
                                delta_t <- Inf
                        } else {
                                delta_t <- t0 - t - times_desc[idx + 1]
                        }
                        
                        p_coal <- 1 - exp(-rate * delta_t)
                        r <- runif(1, 0, 1)
                        
                        if (r <= p_coal) {
                                r_w <- runif(1, 0, 1)
                                w_t <-
                                        inv_t_conditional_exp(r_w, rate, delta_t)
                                
                                log_lh <-
                                        log_lh + exp.loglh(rate, w_t) - log(choose(extant_lineages, 2))
                                
                                t <- t + w_t
                                extant_lineages <-
                                        extant_lineages - 1
                                coalescent_times[coal_idx] <- t0 - t
                                coal_idx <- coal_idx + 1
                        } else {
                                log_lh <- log_lh + log((1 - p_coal))
                                t <- t0 - times_desc[idx + 1]
                                extant_lineages <-
                                        extant_lineages + 1
                                future_lineages <-
                                        future_lineages - 1
                                idx <- idx + 1
                        }
                }
        }
        
        ### add return likelihood
        coal_sample <- list(coalescent_times = coalescent_times,
                            log_likelihood = log_lh)
        
        return(coal_sample)
}

#'Compute likelihood of a realisation of a homogeneous coalescent.
#' 
#' @param sampling_times times of leaves.
#' @param pop_size Effective population size \code{Neg}.
#' @return the log likelihood for the given parameters.
#' @export
homogenous_coal.log_lh <- function(sampling_times,
                                   coalescent_times,
                                   pop_size) {
        coal_times_desc <- coalescent_times[order(-coalescent_times)]
        times_desc <- sampling_times[order(-sampling_times)]
        
        n_sample <- length(times_desc)
        n_coal <- length(coal_times_desc)
        
        assertthat::are_equal(n_sample - 1, n_coal)
        
        coal_idx <- 1
        sample_idx <- 1
        
        t <- 0
        t0 <- times_desc[1]
        
        log_lh <- 0
        j <- 1
        
        while (coal_idx <=  n_coal) {
                if (sample_idx < n_sample &&
                    coal_times_desc[coal_idx] < times_desc[sample_idx + 1]) {
                        delta_t <- t0 -
                                t -
                                times_desc[sample_idx + 1]
                        sample_idx <- sample_idx + 1
                        rate <- choose(j, 2) / pop_size
                        log_lh <-
                                log_lh + poi_0.loglh(rate, delta_t)
                        
                        j <- j + 1
                        t <- t + delta_t
                } else {
                        delta_t <- t0 - t - coal_times_desc[coal_idx]
                        rate <- choose(j, 2) / pop_size
                        log_lh <- log_lh + exp.loglh(rate, delta_t) - log(choose(j, 2))
                        
                        j <- j - 1
                        t <- t + delta_t
                        
                        coal_idx <- coal_idx + 1
                }
        }
        return(log_lh)
}