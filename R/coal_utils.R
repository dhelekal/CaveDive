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