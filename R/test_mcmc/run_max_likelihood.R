library(CaveDive)
library(ggplot2)
library(rmutil)

set.seed(1)

sam <- runif(100, 0, 10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

gt.N = runif(1, 1, 100)
gt.lambda = runif(1, 0.1, 10)

x0 <- c(runif(1,-30, 30), runif(1,0.1,200))

pdf(file="tree.pdf")
co <- plot_exp_growth(sam, gt.lambda, gt.N)
dev.off()

times <- co$coalescent_times
times <- times[order(-times)]

log_lh <- function(x){
  lambda <- x[1]
  n <- x[2]

  if (n > 0){
    prior_lambda <- dlaplace(lambda, m=0, s=1, log=TRUE)
    prior_n <- dexp(n, rate = 1, log = TRUE)
    lh <- -exponential_coalescent_loglh(sam, times, lambda, n) - prior_lambda - prior_n
  } else {
    lh <- Inf
  }
  return(lh)
}

out<-optim(x0, log_lh, control = list(maxit = 2000000)) #method = "L-BFGS-B", lower = c(-100, 1), upper = c(100, 10000),)
print(out)
