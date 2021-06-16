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

compute_ci <- function(x, conf=0.95) {
  ci <- c()
  x_ord <- order(x)
  if(length(x)%%2==0) {
    l<-length(x)/2
    p1 <- x[x_ord][(l+1):length(x)]
    p2 <- x[x_ord][1:l]
  } else {
    l <-floor(length(x)/2)
    p1 <- x[x_ord][(l+2):length(x)]
    p2 <- x[x_ord][1:l]
  }
  ci[1] <- p2[l*(1-conf)]
  ci[2] <- p1[l*conf]
  return(ci)
}

base_dir <- "../GPSC_Trees/GPSC18"
expansions <- readRDS(file = paste0(base_dir, "/mcmc_out/expansions.rds"))
expansions <- discard_burn_in(expansions, proportion=burn_in)

corr_df <- read.csv(paste0(base_dir, "/microreact-project-gpsGPSC18-data.csv"),row.names=1)
corr_df <- as.data.frame(corr_df)

corr_df1 <- corr_df[,c("In_Silico_Serotype"),drop=F]
colnames(corr_df1) <- "Serotype"
#corr_df <- corr_df[,c("Continent"),drop=F]

#corr_df$erm <- sapply(corr_df$erm, function(x) if (x=="neg") "absent" else x)
#colnames(corr_df) <- "erm gene"

#for(c in colnames(corr_df)) corr_df[,c] <- as.logical(corr_df[,c])

#corr_df <- corr_df[expansions$phylo_preprocessed$phy$tip.label,]

which_br <- 385
t_max<-300
eval_pts <- 100

event.df <- expansions$expansion_data
mcmc.df <- expansions$model_data

mode_br_df <- event.df[which(event.df$br==which_br),]
mode_br_mcmc_df <- mcmc.df[mcmc.df$it %in% mode_br_df$it, ]
min_x <- min(mode_br_df$time)
if(is.null(t_max)) {
   max_x <- 0.3*abs(min_x)
} else {
   max_x <- t_max
}
X <- seq(from=min_x, to=max_x, length.out=eval_pts)
Y_med <- rep(0, eval_pts)
Y_min <- rep(0, eval_pts)
Y_max <- rep(0, eval_pts)

funcs <- lapply(c(1:nrow(mode_br_df)), 
   function (i) function (s) 1/sat.rate(s, mode_br_df$K[i], (1/mode_br_df$t_mid[i])**2, mode_br_df$time[i]))

for (i in c(1:eval_pts)) {
   f_vals <- sapply(c(1:length(funcs)), function(j) funcs[[j]](-X[i]))
   Y_med[i] <- median(f_vals)
   ci <- compute_ci(f_vals)
   Y_min[i] <- ci[1]
   Y_max[i] <- ci[2]
}

Y_med_neut <- rep(0, eval_pts)
Y_min_neut <- rep(0, eval_pts)
Y_max_neut <- rep(0, eval_pts)

neut_f  <- lapply(c(1:nrow(mode_br_mcmc_df)), 
   function (i) function (s) 1/constant.rate(s, mode_br_mcmc_df$N[i]))

for (i in c(1:eval_pts)) {
   f_vals <- sapply(c(1:length(neut_f)), function(j) neut_f[[j]](-X[i]))
   Y_med_neut[i] <- median(f_vals)
   ci <- compute_ci(f_vals)
   Y_min_neut[i] <- ci[1]
   Y_max_neut[i] <- ci[2]
}


df <- data.frame(t=c(X,X), 
        y_med=c(Y_med, Y_med_neut),
        y_min=c(Y_min, Y_min_neut), 
        y_max=c(Y_max, Y_max_neut), 
        Population=c(rep("Expansion", eval_pts), rep("Background", eval_pts)))
gg <- ggplot(df, aes(group=Population)) +
geom_ribbon(aes(x=t,ymin=y_min, ymax=y_max, fill=Population), alpha=0.3) +
geom_line(data=subset(df, t <= 0), aes(group=Population, x=t, y=y_med, color=Population), linetype="solid",lwd=2) +
geom_line(data=subset(df, t > 0), aes(group=Population, x=t, y=y_med, color=Population), linetype="longdash",lwd=2) + 
scale_color_brewer(palette="Dark2") +
scale_fill_brewer(palette="Dark2") +
theme_bw() +
xlab("Time") +
ylab("Neg") +
theme(text = element_text(size=20), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio=1,
      legend.position = c(0.2, 0.8))

png("fig_gpsc18_corr.png", width=1600, height=1600)
plot(expansions, mode="persistence", k_modes=1, correlates=list(corr_df1))
dev.off()

png("fig_gpsc18_param.png", width=1600, height=1600)
plot(expansions, mode="modes", k_modes=1)
dev.off()

png("fig_gpsc18_trace.png", width=1600, height=1600)
plot(expansions, mode="traces")
dev.off()

png("fig_gpsc18_summary.png", width=1600, height=1600)
plot(expansions, mode="summary",k_modes=1)
dev.off()


png("fig_gpsc18_mtrace.png", width=1600, height=1600)
plot(expansions, mode="mtraces",k_modes=1)
dev.off()

png("fig_gpsc18_plotfn.png", width=800, height=800)
plot(gg)
dev.off()




