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

base_dir <- "./gpsc23"
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
    saveRDS(expansions, file = paste0(base_dir, "/expansions.rds"))
} else {
    expansions <- readRDS(file = paste0(base_dir, "/expansions.rds"))
}

expansions <- discard_burn_in(expansions, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/microreact-project-gpsGPSC23-data.csv"),row.names=1)
corr_df <- as.data.frame(corr_df)

corr_df <- corr_df[,c("Continent"),drop=F]
#corr_df$erm <- sapply(corr_df$erm, function(x) if (x=="neg") "absent" else x)
#colnames(corr_df) <- "erm gene"

#for(c in colnames(corr_df)) corr_df[,c] <- as.logical(corr_df[,c])

#corr_df <- corr_df[expansions$phylo_preprocessed$phy$tip.label,]
png("fig_gpsc23_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence", k_modes=3, correlates=corr_df, 
                                     corr_axis_title="",
                                     corr_legend_title="")
dev.off()

png("fig_gpsc23_param.png", width=1600, height=1600)
plot(expansions, mode="modes", k_modes=3, 
                                     corr_axis_title="",
                                     corr_legend_title="")
dev.off()

