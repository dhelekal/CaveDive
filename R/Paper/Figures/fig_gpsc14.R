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

base_dir <- "../GPSC_Trees/GPSC14"
expansions <- readRDS(file = paste0(base_dir, "/mcmc_out/expansions.rds"))
expansions <- discard_burn_in(expansions, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/microreact-project-gpsGPSC14-data.csv"),row.names=1)
corr_df <- as.data.frame(corr_df)

corr_df1 <- corr_df[,c("In_Silico_ST__autocolour"),drop=F]
colnames(corr_df1) <- c("ST")
corr_df1$ST <- sapply(corr_df1$ST, function(x) if(x==6279 || x==2059) x else "Other")
corr_df1$ST <- factor(corr_df1$ST) 
#corr_df <- corr_df[,c("Continent"),drop=F]

#corr_df$erm <- sapply(corr_df$erm, function(x) if (x=="neg") "absent" else x)
#colnames(corr_df) <- "erm gene"

#for(c in colnames(corr_df)) corr_df[,c] <- as.logical(corr_df[,c])

#corr_df <- corr_df[expansions$phylo_preprocessed$phy$tip.label,]
png("fig_gpsc14_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence", k_modes=1, correlates=list(corr_df1))
dev.off()

png("fig_gpsc14_param.png", width=1600, height=1600)
plot(expansions, mode="modes", k_modes=1)
dev.off()

