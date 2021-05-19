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

base_dir <- "./ebola"
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
    tree <- read.tree(file = paste0(base_dir,"/","ebola.nwk"))
    tree <- makeNodeLabel(tree)
    tree <- ladderize(tree, F)
    expansions <- run_expansion_inference(tree, priors, 1, n_it=n_it, thinning=thinning)
    saveRDS(expansions, file = paste0(base_dir, "/expansions.rds"))
} else {
    expansions <- readRDS(file = paste0(base_dir, "/expansions.rds"))
}

expansions <- discard_burn_in(expansions, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/microreact-ebola.csv"),row.names=2)
corr_df[] <- lapply(corr_df, factor)
corr_df <- corr_df[c("location")]

df2 <- do.call(rbind,lapply(corr_df$location, function (x) levels(corr_df$location)==x))
rownames(df2) <- rownames(corr_df)
colnames(df2) <- levels(corr_df$location)
df2 <- as.data.frame(df2)
png("fig_ebola_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence", correlates=df2, 
                                     corr_axis_title="Location",
                                     corr_legend_title="Country")
dev.off()
