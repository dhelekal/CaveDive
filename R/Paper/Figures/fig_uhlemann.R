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

set.seed(1)

run_mcmc <- F

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
    tree <- read.tree(file = paste0(base_dir,"/","uhlemann2014.nwk"))
    tree <- makeNodeLabel(tree)
    expansions <- run_expansion_inference(tree, priors, 1, n_it=n_it, thinning=thinning)
    saveRDS(expansions, file = paste0(base_dir, "/expansions.rds"))
} else {
    expansions <- readRDS(file = paste0(base_dir, "/expansions.rds"))
    expansions2 <- readRDS(file = paste0(base_dir, "/expansions_s2.rds"))
}

expansions <- discard_burn_in(expansions, proportion=burn_in)
expansions2 <- discard_burn_in(expansions2, proportion=burn_in)

png("fig_uhlemann_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence",k_modes=3)#, correlates=corr_df, 
                                     #corr_axis_title="Antibiotic",
                                     #corr_legend_title="Resistance")
dev.off()

png("fig_uhlemann_summary.png", width=1600, height=1600)
plot(expansions, mode="summary",k_modes=3)#, correlates=corr_df, 
                                     #corr_axis_title="Antibiotic",
                                     #corr_legend_title="Resistance")
dev.off()

png("fig_uhlemann_summary_s2.png", width=1600, height=1600)
plot(expansions2, mode="summary",k_modes=3)#, correlates=corr_df, 
                                     #corr_axis_title="Antibiotic",
                                     #corr_legend_title="Resistance")
dev.off()

m1 <- mcmc(expansions$model_data[,2:3],thin=expansions$metadata$thinning)
m2 <- mcmc(expansions2$model_data[,2:3],thin=expansions2$metadata$thinning)

grs <- gelman.diag(mcmc.list(m1,m2))
sink("GelmanRubinFigUhlemann.txt")
print(grs$psrf)
print(grs$mpsrf)
sink()