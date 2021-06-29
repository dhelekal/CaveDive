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

base_dir <- "../Trees/GPSC9"
expansions <- readRDS(file = paste0(base_dir, "/mcmc_out/expansions.rds"))
expansions2 <- readRDS(file = paste0(base_dir, "/mcmc_out/expansions_s2.rds"))

expansions <- discard_burn_in(expansions, proportion=burn_in)
expansions2 <- discard_burn_in(expansions2, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/microreact-project-gpsGPSC9-data.csv"),row.names=1)
corr_df <- as.data.frame(corr_df)

corr_df1 <- corr_df[,"erm",drop=F]
corr_df1$erm <- sapply(corr_df1$erm, function(x) if (x=="neg") "No erm gene" else x)
corr_df1$erm <- factor(corr_df1$erm, levels=c("No erm gene", "ermB1"), ordered=T)
colnames(corr_df1) <- "erm gene"
corr_df2 <- corr_df[,"Vaccine_Status",drop=F]
corr_df2$Vaccine_Status <- sapply(corr_df2$Vaccine_Status, function(x) if (x=="NVT") "Non-vaccine type" else "Vaccine type")
corr_df2$Vaccine_Status <- factor(corr_df2$Vaccine_Status, levels=c("Non-vaccine type", "Vaccine type"), ordered=T)
colnames(corr_df2) <- "type"

corr_df3 <- corr_df[,"Continent",drop=F]
colnames(corr_df3) <- "Continent"


#for(c in colnames(corr_df)) corr_df[,c] <- as.logical(corr_df[,c])

#corr_df <- corr_df[expansions$phylo_preprocessed$phy$tip.label,]
png("fig_gpsc9_corr.png", width=2000, height=2000)
plot(expansions, mode="persistence", k_modes=3, correlates=list(corr_df3, corr_df1, corr_df2), 
                                     corr_axis_title=list(),
                                     corr_legend_title=list(), no_y_text=T)
dev.off()

png("fig_gpsc9_param.png", width=1600, height=1600)
plot(expansions, mode="modes", k_modes=3)
dev.off()

png("fig_gpsc9_summary.png", width=1600, height=1600)
plot(expansions, mode="summary", k_modes=3)
dev.off()

png("fig_gpsc9_trace.png", width=1600, height=1600)
plot(expansions, mode="traces")
dev.off()

png("fig_gpsc9_mtrace.png", width=1600, height=1600)
plot(expansions, mode="mtraces",k_modes=3)
dev.off()

png("fig_gpsc9_popfn.png", width=1600, height=800)
plot(expansions, mode="popfunc",k_modes=3,t_max=c(50,50,50))+xlab("Time (Years)")

dev.off()

m1 <- mcmc(expansions$model_data[,2:3],thin=expansions$metadata$thinning)
m2 <- mcmc(expansions2$model_data[,2:3],thin=expansions2$metadata$thinning)

grs <- gelman.diag(mcmc.list(m1,m2))
sink("GelmanRubinGPSC9.txt")
print(grs$psrf)
print(grs$mpsrf)
print(effectiveSize(m1))
sink()


