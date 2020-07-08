library(CaveDive)

gt.N = 1000
gt.lambda = 40

sam <- seq(-2,0, by=0.1)
sam <- sam[order(-sam)]

co <- plot_exp_growth(sam, gt.lambda, gt.N)
times <- co$coalescent_times
times <- times[order(-times)]

log_lh <- function(n, lambda){
  return(exponential_coalescent_loglh(sam, times, n, lambda))
}