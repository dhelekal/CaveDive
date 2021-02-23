#!/usr/bin/env Rscript
library(CaveDive)
library(rjson)

source("make_report.R")

args = commandArgs(trailingOnly=TRUE)

n_it <- as.integer(args[1])
thinning <- as.integer(args[2])
simulation_in <- as.character(args[3])
dir_out <- as.character(args[4])

sim_data <- fromJSON(file=simulation_in)
sim_data <- lapply(sim_data, unlist)
sim_data <- lapply(sim_data, as.numeric)

seed <- sim_data$seed
n_exp <- sim_data$n_exp+1
n_tips <- sim_data$n_tips
tip_times <- sim_data$tip_times
tip_colours <- sim_data$tip_colours
div_times <- sim_data$div_times
exp_probs <- sim_data$exp_probs

set.seed(seed)
setwd(file.path(".", dir_out))

N <- sim_data$N
K <- sim_data$K
A <- sim_data$A

div_cols <- c(1:n_exp)

sim_tree <- clonal_tree_process.simulate_tree(n_exp, N, K, A, tip_times, tip_colours, div_times, div_cols, exp_probs)
co <- sim_tree$co

sim <- sim_data
sim$co <- co

tree.div.str <- build_coal_tree.structured(tip_times, co$times, tip_colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)

pdf(file="tree.pdf", width=5, height=5)
tree.plt <- plot_structured_tree(tree.div, n_exp)
plot(tree.plt)
dev.off()

tree.str <- build_coal_tree.structured(tip_times, co$times, tip_colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE, aux_root = FALSE)
tree <- read.tree(text = tree.str$full)
pre <- structured_coal.preprocess_phylo(tree) 

r_mean <- 0 ## growth rate lognormal prior mean
K_mean <- 4 ## carrying capacity rate lognormal prior mean

r_sd <- 1 ## growth rate lognormal prior sd
K_sd <- 4 ## carrying capacity  rate lognormal prior sd

prior_i <- function(x) dpois(x, 1, log = TRUE) ### poisson 1 prior

prior_N <- function(x) dlnorm(x, meanlog = K_mean, sdlog = K_sd, log = TRUE)
prior_N.sample <- function() rlnorm(1, meanlog = K_mean, sdlog = K_sd) 

prior_r <- function(x) dlnorm(x, meanlog = r_mean, sdlog = r_sd, log = TRUE) 
prior_r.sample <- function() rlnorm(1, meanlog = r_mean, sdlog = r_sd) 

prior_K <- function(x) dlnorm(x, meanlog = K_mean, sdlog = K_sd, log = TRUE)
prior_K.sample <- function() rlnorm(1, meanlog = K_mean, sdlog = K_sd) 

prior_t <- function(x) {
       if (all(x < max(pre$nodes.df$times)) && all(x > min(pre$nodes.df$times))) {
              out <- length(x)*log(1/abs(max(pre$nodes.df$times)-min(pre$nodes.df$times)))
       } else {
              out <- -Inf 
       }
       return(out) ### Uniform time prior
} 

prior_t.sample <- function() runif(1, min(pre$nodes.df$times), max(pre$nodes.df$times)) ### Uniform time prior

o <- outbreaks_infer(tree, prior_i,  prior_N,  prior_N.sample, 
                     prior_r, prior_r.sample,  prior_K,  prior_K.sample,  prior_t,
                     prior_t.sample, 2, n_it=n_it, thinning=thinning, debug=F)

make_report(o, pre, sim)