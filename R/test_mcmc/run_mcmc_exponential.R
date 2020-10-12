library(CaveDive)
library(ggplot2)
library(viridis)
library(rmutil)

set.seed(123)

sam <- runif(100,0,10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

t0 <- -50
r <- 0.1
K <- 100
           
Neg_t <- function (s) return(sat.rate(s, K, r ,t0))
Neg_t.int <- function (t, s) return (sat.rate.int(t, s, K, r, t0))
Neg_t.inv_int <- function(t, s) return(sat.rate.int_inv(t, s, K, r , t0))
            
co <-
     inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
     log_lh_tree <- co$log_likelihood
     times <- co$coalescent_times

#sam <- runif(100, 0, 10)
#sam <- sam - max(sam)
#sam <- sam[order(-sam)]

#gt.N = runif(1, 1, 100)
#gt.lambda = runif(1, 0.1, 10)

x0 <- c(runif(1,-30, 30), runif(1,0.1,400))

#pdf(file="tree.pdf")
#co <- plot_exp_growth(sam, gt.lambda, gt.N)
#dev.off()

times <- co$coalescent_times
times <- times[order(-times)]

log_lh <- function(x){
  r <- x[1]
  K <- x[2]

  if (n > 0){
    prior_lambda <- dlaplace(r, m=0, s=0.2, log=TRUE)
    prior_n <- dexp(K, rate = 100, log = TRUE)
    lh <- sat_coalescent_loglh(sam, times, t0, r, K, 0) + prior_lambda + prior_n
  } else {
    lh <- -Inf
  }
  return(lh)
}

mle_log_lh <- function(x){
  lambda <- x[1]
  n <- x[2]

  if (n > 0){
    lh <- exponential_coalescent_loglh(sam, times, lambda, n)
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

png(file="trace.png", width=600, height=600)
plt <- ggplot(o.df, aes(x=n, y=lambda, colour=iteration)) +
       geom_line(alpha = 0.3)+
       geom_point(alpha = 0.3)+
       geom_vline(xintercept = gt.N, colour="red", linetype = "longdash") +
       geom_hline(yintercept = gt.lambda, colour="red", linetype = "longdash") +
       geom_hline(yintercept = MLE$par[1], colour="orange", linetype = "longdash") +
       geom_hline(yintercept = MAP$par[1], colour="blue", linetype = "longdash") + 
       geom_vline(xintercept = MLE$par[2], colour="orange", linetype = "longdash") +
       geom_vline(xintercept = MAP$par[2], colour="blue", linetype = "longdash") +
       scale_color_viridis() + theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()

o.df <- o.df[burn_in:n_it, ]

pdf(file="lambdahist.pdf", width = 5, height = 5)
plt <- ggplot(o.df, aes(x=lambda)) +
       geom_histogram(colour="darkgreen", fill="white", binwidth = 0.01) +
       geom_vline(xintercept = MLE$par[1], colour="orange", linetype = "longdash") +
       geom_vline(xintercept = MAP$par[1], colour="blue", linetype = "longdash") + 
       theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file="nhist.pdf", width=5, height=5)
plt <- ggplot(o.df, aes(x=n)) + 
       geom_histogram(colour="darkgreen", fill="white", binwidth = 1) +
       geom_vline(xintercept = MLE$par[2], colour="orange", linetype = "longdash") +
       geom_vline(xintercept = MAP$par[2], colour="blue", linetype = "longdash") +
       theme(aspect.ratio=1)
plot(plt)
dev.off()

png(file="trace_zoom.png", width=600, height=600)
plt <- ggplot(o.df, aes(x=n, y=lambda, colour=iteration)) +
       geom_line(alpha = 0.3)+
       geom_point(alpha = 0.3)+
       geom_vline(xintercept = gt.N, colour="red", linetype = "longdash")+
       geom_hline(yintercept = gt.lambda, colour="red", linetype = "longdash") +
       geom_hline(yintercept = MLE$par[1], colour="orange", linetype = "longdash") +
       geom_hline(yintercept = MAP$par[1], colour="blue", linetype = "longdash") + 
       geom_vline(xintercept = MLE$par[2], colour="orange", linetype = "longdash") +
       geom_vline(xintercept = MAP$par[2], colour="blue", linetype = "longdash") +
       scale_color_viridis() + theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()
