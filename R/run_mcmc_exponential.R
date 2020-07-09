library(CaveDive)
library(ape)

gt.N = 100
gt.lambda = 1000

sam <- rep(0,50)
sam <- sam[order(-sam)]

co <- plot_exp_growth(sam, gt.lambda, gt.N)
times <- co$coalescent_times
times <- times[order(-times)]

log_lh <- function(x){
  lambda <- x[1]
  n <- x[2]
  if (lambda > 0 && n > 0) {
    lh <- exponential_coalescent_loglh(sam, times, lambda, n)
    } else {
      lh <- -Inf
    }
  return(lh)
}

prop.sampler <-function (x_prev){
  return(c(rnorm(1, mean=x_prev[1], 1), rnorm(1, mean=x_prev[2], 1)))
}

proposal.cond_lh <-function(x_prev, x_cand){
  return(log(dnorm(x_cand[1],x_prev[1],1)) + log(dnorm(x_cand[2],x_prev[2],1)))
}

x0 <- c(1,1)
n_it <- 5000000
burn_in <-500000

o <- run_mcmc(log_lh, proposal.cond_lh, prop.sampler, x0, n_it)

o.df <- as.data.frame(t(as.data.frame(o)))

colnames(o.df) <- c("lambda", "n")
o.df <- o.df[burn_in:n_it, ]

pdf(file="lambdahist.pdf")
hist(o.df$lambda, col="darkgreen")
dev.off()

pdf(file="nhist.pdf")
hist(o.df$n, col="darkgreen")
dev.off()

