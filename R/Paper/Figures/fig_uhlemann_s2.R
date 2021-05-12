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
library(coda)

set.seed(2)

run_mcmc <- T

n_it <- 1e8
thinning <- n_it/1e4
burn_in <- 0.3

base_dir <- "./uhlemann2014"
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
    tree <- read.tree(file = paste0(base_dir,"/","dated.nwk"))
    tree <- makeNodeLabel(tree)
    tree <- ladderize(tree, F)
    expansions <- run_expansion_inference(tree, priors, 1, n_it=n_it, thinning=thinning)
    saveRDS(expansions, file = paste0(base_dir, "/expansions_s2.rds"))
} else {
    expansions <- readRDS(file = paste0(base_dir, "/expansions_s2.rds"))
}

expansions <- discard_burn_in(expansions, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/amr-profile.csv"),row.names=1)
corr_df <- corr_df[, c("AMI","TOB","KAN","MET","PEN","ERY","CIP")]
corr_df[] <- lapply(corr_df, factor)
png("fig_uhlemann_corr_s2.png", width=1600, height=1600)
plot(expansions, mode="persistence",k_modes=3, correlates=corr_df, 
                                     corr_axis_title="Antibiotic",
                                     corr_legend_title="Resistance")
dev.off()

png("fig_uhlemann_summary_s2.png", width=1600, height=1600)
plot(expansions, mode="summary",k_modes=3)
dev.off()# 
