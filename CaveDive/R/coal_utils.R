inv_t_conditional_exp <- function(u, rate, delta_t) {
  return((-1 / rate) * log(1 - u * exp.prob(rate, delta_t)))
}

exp.loglh <- function(rate, t) {
  return(log(rate) -
           rate * t)
}

exp.prob <- function(rate, t) {
  return(1 - exp(-rate * t))
}

poi_0.loglh <- function(rate, t) {
  return(-rate * t)
}

inhomogenous_exp.loglh <- function(rate, rate.int, t, s) {
  return(log(rate(t + s)) - rate.int(t, s))
}

inhomogenous_exp.prob <- function(rate.int, t, s) {
  return(1 - exp(-rate.int(t, s)))
}

inv_t_inhomogenous_exp_conditional <- function(rate.int, rate.inv_int , exp.rate, t, s) {
    u <- runif(1, 0, 1)
    Q <- inhomogenous_exp.prob(rate.int, t, s)
    wt <- (-1 / exp.rate) * log(1 - u * Q)
    return(rate.inv_int(t , wt))
  }

inhomogenous_poi_0.loglh <- function(rate.int, t, s) {
  return(-rate.int(t, s))
}

#'Rescale a homogenous process to one with exponentially growing population.
#' 
#' @param c_t coalescent times.
#' @param beta growth rate.
#' @param t0 time of most recent leave.
#' @return rescaled waiting times.
#' @export
rescale_to_exponential <- function(c_t, beta, t0) {
  w_t.exp <- rep(1, length(c_t))
  
  t <- t0
  t.r <- t0
  
  v <- 0
  
  for (j in c(1:length(c_t))) {
    w <- t - c_t[j]
    w_t <- 1 / beta * log(1 + beta * w * exp(-beta * v))
    w_t.exp[j] <- t.r - w_t
    t.r <- w_t.exp[j]
    v <- v + w_t
    t <- c_t[j]
  }
  return(w_t.exp)
}

#'Precompute waiting times and lineage counts for intervals between events in a homogenous process.
#' 
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events / internal nodes.
#' @return a list of two entries \code{intervals} waiting times, \code{lineages} lineage counts.
#' @export
transform_to_intervals <-function(sampling_times, coalescent_times) {
  coal_times_desc <- coalescent_times[order(-coalescent_times)]
  times_desc <- sampling_times[order(-sampling_times)]
        
  n_sample <- length(times_desc)
  n_coal <- length(coal_times_desc)
        
        
  coal_idx <- 1
  sample_idx <- 1
  idx <- 1
        
  t <- 0
  t0 <- times_desc[1]

  j <- 1
        
  out.intervals <- rep(0, 2*n_sample-1)
  out.lineages <- rep(0, 2*n_sample-1)

  while (coal_idx <= n_coal) {
                if (sample_idx < n_sample &&
                    coal_times_desc[coal_idx] < times_desc[sample_idx + 1]) {
                        delta_t <- t0 -
                                t -
                                times_desc[sample_idx + 1]

                        out.intervals[idx] <- delta_t
                        out.lineages[idx] <- j

                        j <- j + 1
                        t <- t + delta_t

                        idx <- idx + 1
                        sample_idx <- sample_idx + 1

                } else {
                        delta_t <- t0 - t - coal_times_desc[coal_idx]

                        out.intervals[idx] <- delta_t
                        out.lineages[idx] <- j
                        
                        j <- j - 1
                        t <- t + delta_t
                        
                        idx <- idx + 1
                        coal_idx <- coal_idx + 1
                }
        }
        return(list(intervals=out.intervals, lineages=out.lineages))
}