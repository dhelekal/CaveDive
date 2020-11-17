library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

#set.seed(12345678)

n_tips <- 100
n_exp <- 3
concentration <- rep(2,n_exp)

K_mean <- 5
K_sd <- 0.5

para <- clonal_tree_process.simulate_params(n_exp, n_tips, concentration, K_mean, K_sd, sampling_scale=c(0,20), r_mean=0, r_sd=1, time_mean_sd_mult=8)
out <- clonal_tree_process.simulate_tree(para$n_exp, para$N, para$K, para$A, para$sampling_times, para$tip_colours, para$div_times, para$div_cols, para$exp_probs)

co <- out$co

tree.div.str <- build_coal_tree.structured(para$sampling_times, co$times, para$tip_colours, co$colours, para$div_times, para$div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)
tree.plt <- plot_structured_tree(tree.div, para$n_exp)

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
