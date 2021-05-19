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
base_dir <- "./grad2016"

run_mcmc <- F

n_it <- 1e8
thinning <- n_it/1e4
burn_in <- 0.3

readResist=function() {
  metadata=read.table(paste0(base_dir,"/grad2016.tab"),sep='\t',comment.char='',header=T,as.is=T)
  abx=3:9  #PEN,TET, SPC,CFX,  CRO,  CIP, AZI
  cutoff1=c(1   ,1  ,64 ,0.125,0.125,0.06,0.5 )#EUCAST higher breakpoints - resistant only if >  cutoff
  cutoff2=c(0.06,0.5,64 ,0.125,0.125,0.03,0.25)#EUCAST lower  breakpoints - sensitive only if <= cutoff
  cutoff1=c(1   ,1  ,64 ,0.125,0.06 ,0.5 ,1   )#used   higher breakpoints - resistant only if >  cutoff
  abxnames=colnames(metadata)[abx]
  resist=matrix(1,nrow(metadata),length(abx))
  for (i in 1:length(abx)) {
    resist[which(is.na(metadata[,abx[i]])),i]=NA
    resist[which(metadata[,abx[i]]> cutoff1[i]),i]=2#Resistant
    resist[which(metadata[,abx[i]]<=cutoff2[i]),i]=0#Sensitive
  }
  colnames(resist)=abxnames
  rownames(resist)=metadata[,1]
  return(resist)
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

corr_df <- data.frame(readResist())
png("fig_grad_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence", correlates=corr_df, 
                                     corr_axis_title="Antibiotic",
                                     corr_legend_title="Resistance")
dev.off()

