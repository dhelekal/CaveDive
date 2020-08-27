library(CaveDive)
library(ggplot2)
library(viridis)
library(rmutil)

set.seed(123)

sam <- runif(100, 0, 10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

gt.N = runif(1, 1, 100)
gt.lambda = runif(1, 0.1, 10)

colours <- trunc(runif(100, 1, 4))

N <- rexp(n, rate = 100)
A <- rexp(n, rate = 10)
K <- rexp(n, rate = 100)
div_times <- c(-1*runif(2,5,50), -Inf)

div_cols <- c(1, 2, 3)
rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))
co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)
print(co)

tree.div.str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)
tree.plt <- plot_structured_tree(tree.div, 3)

pdf(file="tree_strucutred.pdf")
plot(tree.plt)
dev.off()

tree.str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
tree <- read.tree(text = tree.str$full)

pre <- structured_coal.preprocess_phylo(tree)

# x[[1]] rates
# x[[2]] K
# x[[3]] neutral size
# x[[4]] div.times
# x[[5]] div.branch
log_lh <- function(x){
  rates <- x[[1]]
  K <- x[[2]]
  N <- x[[3]]
  div.times <- x[[4]]
  div.branch <- x[[5]]

  prior_rates <- sum(dexp(rates, rate = 1, log = TRUE))
  prior_K <- sum(dexp(K, rate = 100, log = TRUE))
  prior_N <- dexp(N, rate = 100, log = TRUE)


  if (n > 0){
    prior_lambda <- dlaplace(lambda, m=0, s=1, log=TRUE)
    prior_n <- dexp(n, rate = 1, log = TRUE)
    lh <- exponential_coalescent_loglh(sam, times, lambda, n) + prior_lambda + prior_n
  } else {
    lh <- -Inf
  }
  return(lh)
}

prop.sampler <-function (x_prev){
  return(c(rnorm(1, mean=x_prev[1], 1), rnorm(1, mean=x_prev[2], 1)))
}

proposal.cond_lh <-function(x_cand, x_prev){
  return(log(dnorm(x_cand[1],x_prev[1],1)) + log(dnorm(x_cand[2],x_prev[2],1)))
}


MLE <- optim(c(1,1), mle_log_lh , control = list(maxit = 2000000, fnscale=-1))
MAP <- optim(c(1,1), log_lh, control = list(maxit = 2000000, fnscale=-1))

n_it <- 1000000
burn_in <-100000

o <- run_mcmc(log_lh, proposal.cond_lh, prop.sampler, x0, n_it)
o.df <- as.data.frame(t(as.data.frame(o)))

colnames(o.df) <- c("lambda", "n")
o.df$iteration <- c(1:nrow(o.df))
