library(CaveDive)
library(rmutil)

set.seed(1)

sam <- runif(100, 0, 10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

gt.N = runif(1, 1, 100)
gt.lambda = runif(1, 0.1, 10)

x0 <- c(runif(1,-30, 30), runif(1,0.1,1000))

co <- plot_exp_growth(sam, gt.lambda, gt.N)
times <- co$coalescent_times
times <- times[order(-times)]

log_lh <- function(x){
  lambda <- x[1]
  n <- x[2]

  if (n > 0){
    prior_lambda <- dlaplace(lambda, m=0, s=1, log=TRUE)
    prior_n <- dexp(n, rate = 1, log = TRUE)
    lh <- exponential_coalescent_loglh(sam, times, lambda, n) + prior_n + prior_lambda
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

n_it <- 3000000
burn_in <-100000

o <- run_mcmc(log_lh, proposal.cond_lh, prop.sampler, x0, n_it)

o.df <- as.data.frame(t(as.data.frame(o)))

colnames(o.df) <- c("lambda", "n")
o.df <- o.df[burn_in:n_it, ]

pdf(file="lambdahist.pdf")
hist(o.df$lambda, col="darkgreen")
abline(v=gt.lambda)
dev.off()

pdf(file="nhist.pdf")
hist(o.df$n, col="darkgreen")
abline(v=gt.N)
dev.off()

