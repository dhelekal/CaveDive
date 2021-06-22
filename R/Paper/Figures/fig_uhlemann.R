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

burn_in <- 0.3

base_dir <- "../Trees/Uhlemann2014"

expansions <- readRDS(file = paste0(base_dir, "/mcmc_out/expansions.rds"))
expansions2 <- readRDS(file = paste0(base_dir, "/mcmc_out/expansions_s2.rds"))

expansions <- discard_burn_in(expansions, proportion=burn_in)
expansions2 <- discard_burn_in(expansions2, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/metadata.csv"),row.names=1)
corr_df <- data.frame(corr_df)
corr_df <- corr_df[, c("MRSA", "ACME")]

for(c in colnames(corr_df)) corr_df[,c] <- sapply(corr_df[,c], function (x) if (x==0) "absent" else "present")

corr_df <- corr_df[expansions$phylo_preprocessed$phy$tip.label,]
png("fig_uhlemann_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence",k_modes=3, correlates=list(corr_df))
dev.off()

png("fig_uhlemann_summary.png", width=1600, height=1600)
plot(expansions, mode="summary",k_modes=3)
dev.off()# 

png("fig_uhlemann_modes.png", width=1600, height=1600)
plot(expansions, mode="modes",k_modes=3)
dev.off()# 

png("fig_uhlemann_trace.png", width=1600, height=1600)
plot(expansions, mode="traces")
dev.off()

png("fig_uhlemann_mtrace.png", width=1600, height=1600)
plot(expansions, mode="mtraces",k_modes=3)
dev.off()

m1 <- mcmc(expansions$model_data[,2:3],thin=expansions$metadata$thinning)
m2 <- mcmc(expansions2$model_data[,2:3],thin=expansions2$metadata$thinning)

grs <- gelman.diag(mcmc.list(m1,m2))
sink("GelmanRubinUhlemann.txt")
print(grs$psrf)
print(grs$mpsrf)
print(effectiveSize(m1))
sink()



