#'Simulate inhomogenous coalescent,
#' 
#' @param sampling_times times of leaves.
#' @param Neg.rate function equal to 1/Neg(t).
#' @param Neg.rate.int integral of 1/Neg(t) from \code{t} to \code{t+s}.
#' @param Neg.rate.int_inv function that inverts the integral of 1/Neg(t) from \code{t} to \code{t+s} for \code{s}.
#' @return A list consisting of the simulate coalescent times \code{coalescent_times} and the log-likelihood of the simulation \code{log_likelihood}.
#' @export
inhomogenous_coal.simulate <- function(sampling_times,
                                       Neg.rate,
                                       Neg.rate.int,
                                       Neg.rate.int_inv=NULL) {
  log_lh <- 0
  
  times_desc <- sampling_times[order(-sampling_times)]
  future_lineages <- length(sampling_times) - 1
  extant_lineages <- 1
  coalescent_times <- rep(0, future_lineages)
  t <- 0
  t0 <- times_desc[1]
  idx <- 2
  coal_idx <- 1

  if(is.null(Neg.rate.int_inv)){
    Neg.rate.int_inv <- function(t,s) inv_rates(function (s) Neg.rate.int(t,s), s, dfun=Neg.rate)
  }
  
  while (future_lineages > 0 || extant_lineages > 1) {
    if (extant_lineages < 2) {
      #If one lineage continue to next sampling event
      t <- t0 - times_desc[idx]
      extant_lineages <- extant_lineages + 1
      future_lineages <- future_lineages - 1
      idx <- idx + 1
    } else {
      #Pick waiting time
      c <- choose(extant_lineages, 2)
      if (idx> length(sampling_times)) {
        s <- Inf
      } else {
        s <- t0 - t - times_desc[idx]
      }
      
      p_coal <-
        inhomogenous_exp.prob(function(t, s)
          c * Neg.rate.int(t, s), t, s)
      
      r <- runif(1, 0, 1)
      
      if (r <= p_coal) {
        w_t <- inv_t_inhomogenous_exp_conditional(function(t, s)
          c * Neg.rate.int(t, s),
          Neg.rate.int_inv,
          c,
          t,
          s)
        
        log_lh <- log_lh + inhomogenous_exp.loglh(function(s)
            c * Neg.rate(s),
            function(t, s)
              c * Neg.rate.int(t, s), t, w_t) - log(c)
        
        
        t <- t + w_t
        extant_lineages <- extant_lineages - 1
        coalescent_times[coal_idx] <- t0 - t
        coal_idx <- coal_idx + 1
      } else {
        log_lh <- log_lh + log(1 - p_coal)
        t <- t0 - times_desc[idx]
        extant_lineages <- extant_lineages + 1
        future_lineages <- future_lineages - 1
        idx <- idx + 1
      }
    }
  }
  
  ### add return likelihood
  coal_sample <- list(coalescent_times = coalescent_times,
                      log_likelihood = log_lh)
  
  return(coal_sample)
  
}


#'Compute likelihood of a realisation of an inhomogeneous coalescent.
#' 
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events.
#' @param Neg.rate function equal to 1/Neg(t).
#' @param Neg.rate.int integral of 1/Neg(t) from \code{t} to \code{t+s}.
#' @return the log likelihood for the given parameters.
#' @export
inhomogenous_coal.log_lh <- function(sampling_times,
                                     coalescent_times,
                                     Neg.rate,
                                     Neg.rate.int) {
  coal_times_desc <- coalescent_times[order(-coalescent_times)]
  times_desc <- sampling_times[order(-sampling_times)]
  
  n_sample <- length(times_desc)
  n_coal <- length(coal_times_desc)
    
  coal_idx <- 1
  sample_idx <- 1
  
  t <- 0
  t0 <- times_desc[1]
  
  log_lh <- 0
  j <- 1
  
  while (coal_idx <= n_coal) {
    rate.int <- function (t, s)
      Neg.rate.int(t, s) * choose(j, 2)
    if (sample_idx < n_sample &&
        coal_times_desc[coal_idx] < times_desc[sample_idx + 1]) {
      s <- t0 -
        t -
        times_desc[sample_idx + 1]
      sample_idx <- sample_idx + 1
      log_lh <- log_lh + inhomogenous_poi_0.loglh(rate.int, t, s)
      
      j <- j + 1
      t <- t + s
    } else {
      if (j <= 1) {
        warning(
          "Number of extant lineages at time of
                coalescent event less than 2.\n
                Supplied values not valid coalescent times.\n
                Numerical precision may have been exhausted"
        )
      }
      
      s <- t0 - t - coal_times_desc[coal_idx]
      rate <- function(s)
        Neg.rate(s) * choose(j, 2)
      log_lh <-
        log_lh + inhomogenous_exp.loglh(rate, rate.int, t, s) - log(choose(j, 2))

      
      j <- j - 1
      t <- t + s
      
      coal_idx <- coal_idx + 1
    }
  }
  return(log_lh)
}
