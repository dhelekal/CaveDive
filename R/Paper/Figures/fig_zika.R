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
burn_in <- 0.2

base_dir <- "../Trees/zika/mcmc_out"

expansions <- readRDS(file = paste0(base_dir, "/expansions.rds"))
expansions <- discard_burn_in(expansions, proportion=burn_in)

png("fig_zika_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence",k_modes=1)
dev.off()

png("fig_zika_pop.png",width=1600, height=900)
plot(expansions, mode="modes", k_modes=1)
dev.off()