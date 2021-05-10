library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(egg)

set.seed(1)

run_mcmc <- F

n_it <- 1e7
thinning <- n_it/1e4
burn_in <- 0.2

base_dir <- "./zika"
if(run_mcmc) {
    dir.create(file.path(base_dir))
}

priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=3, 
                    N_sd_log=3, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/2)

if (run_mcmc) {
    tree <- read.tree(file = paste0(base_dir,"/","zikaRed.nwk"))
    tree <- makeNodeLabel(tree)
    expansions <- run_expansion_inference(tree, priors, 1, n_it=n_it, thinning=thinning)
    saveRDS(expansions, file = paste0(base_dir, "/expansions.rds"))
} else {
    expansions <- readRDS(file = paste0(base_dir, "/expansions.rds"))
}

expansions <- discard_burn_in(expansions, proportion=burn_in)

png("fig_zika_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence")
dev.off()

unique_br <- unique(expansions$expansion_data$br)
mode_br <-unique_br[which.max(sapply(unique_br, function(br) length(which(expansions$expansion_data$br==br))))]
br_marginal <- expansions$expansion_data[which(expansions$expansion_data$br==mode_br),]

png("fig_zika_pop.png",width=1600, height=1600)
plot(ggplot(br_marginal, aes(K)) +
         geom_histogram(aes(y = stat(count / sum(count))), bins=100) +
         theme_bw() + 
         theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
               text = element_text(size=20)))
dev.off()