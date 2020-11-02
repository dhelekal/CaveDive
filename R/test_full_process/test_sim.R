library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

set.seed(12345678)

n_tips <- 100
poi_rate <- 3
concentration <- 2

sam <- runif(n_tips, 0, 0.1)

r_mean <- 1
r_sd <- 1

K_mean <- 5
K_sd <- 0.5

time_shape <- 30
time_rate <- 5**(-1)

out <- simulate_outbreaks(poi_rate, concentration, sam, r_mean, r_sd, K_mean, K_sd, time_rate, time_shape)
co <- out$co

tree.div.str <- build_coal_tree.structured(sam, co$times, out$colours, co$colours, out$div_times, out$div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)
tree.plt <- plot_structured_tree(tree.div, out$n_exp)

pdf(file="tree_structured.pdf")
plot(tree.plt)
dev.off()

log_lh <- function(x){
  Neg <- x

  if (Neg > 0){
    prior_neg <- dlnorm(Neg, meanlog = K_mean, sdlog = K_sd, log=TRUE)
    lh <- coalescent_loglh(sam,
                        co$times,
                        Neg,
                        0) + prior_neg
  } else {
    lh <- -Inf
  }
  return(lh)
}

MAP <- optim(c(1), log_lh, control = list(maxit = 2000000, fnscale=-1))

print(paste0("MAP optimisation results for no-expansion model: ", MAP))
