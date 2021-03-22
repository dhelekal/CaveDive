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

set.seed(3)

run_mcmc <- FALSE

n_it <- 2e7
thinning <- n_it/1e4
burn_in <- 1e3 #10%

data_dir <- "./grad2016"
if(run_mcmc) {
    dir.create(file.path(data_dir))
}

priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=3, 
                    N_sd_log=3, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/2)

tree <- read.tree(file = paste0(data_dir,"/grad2016.nwk"))
tree <- makeNodeLabel(tree)
#tree <- ladderize(tree, right = T)
if (run_mcmc) {
    expansions <- run_expansion_inference(tree, priors, 1, n_it=n_it, thinning=thinning)
    pre <- expansions$phylo_preprocessed ### keep a hold of this, it's useful for plotting
    dfs <- mcmc2data.frame(expansions$mcmc_out)
    mcmc.df <- dfs$mcmc.df
    event.df <- dfs$event.df
    write.csv(mcmc.df, paste0(data_dir,"/mcmc_df.csv"))
    write.csv(event.df, paste0(data_dir,"/event_df.csv"))
} else {
    mcmc.df <- read.csv(paste0(data_dir,"/mcmc_df.csv"))
    event.df <- read.csv(paste0(data_dir,"/event_df.csv"))
    pre <- structured_coal.preprocess_phylo(tree, order_edges_by_node_label=F)
}
### discard burn in
mcmc.df <- mcmc.df[burn_in:nrow(mcmc.df),]
event.df <- event.df[which(event.df$it >= burn_in), ]

### marginalise mode model
unique_dims <- unique(mcmc.df$dim)
dim_counts <- sapply(unique_dims, function(i) length(which(mcmc.df$dim==i)))
mode_dim <- unique_dims[which.max(dim_counts)]

mode_dim_marginal <- mcmc.df[which(mcmc.df$dim==mode_dim),]
mode_dim_marginal_it <- mode_dim_marginal$it
event_mode_dim_marginal <- event.df[event.df$it %in% mode_dim_marginal_it,]

unique.br <- unique(event_mode_dim_marginal$br)
branch_counts <- sapply(unique.br, function(x) length(which(event_mode_dim_marginal$br == x)))
br_count_ord <- order(-branch_counts)
k_mode_br <- unique.br[br_count_ord[c(1:mode_dim)]]

mode_br_df <- event_mode_dim_marginal[event_mode_dim_marginal$br %in% k_mode_br, ]
mode_br_mcmc_df <- mode_dim_marginal[mode_dim_marginal$it %in% mode_br_df$it, ]

compute_persistence <- function(pre, df) {
    p_mat <- matrix(data = 0, nrow = pre$n_tips, ncol = pre$n_tips, byrow = FALSE,
       dimnames = list(pre$phy$tip.label,pre$phy$tip.label))
    for(i in unique(df$it)) {
        subs_it <- df[which(df$it == i), ]
        partitions <- extract_lineage_times(pre, pre$phy$node.label[(c(pre$edges.df$node.child[subs_it$br],pre$root_idx)-pre$n_tips)], c(subs_it$time, -Inf), return_partitions=TRUE)$partitions
        for(p in partitions) {
            for(t1 in p){
                p_mat[t1,p] <- p_mat[t1,p] + 1
            }
        }
    }
    p_mat <- p_mat/length(unique(df$it))
    return(p_mat)
}


p_mat <- compute_persistence(pre, event.df)
#ord <- hclust(as.dist(1-p_mat, diag = T, upper = T), method = "ward.D")$order

readResist=function() {
  metadata=read.table(paste0(data_dir,"/grad2016.tab"),sep='\t',comment.char='',header=T,as.is=T)
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

plot_tree<-function(pre,event.df){
   tree <- pre$phy
   freq <- table(event.df$br)

   labs <- c(tree$node.label, tree$tip.label)
   tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
   ids <- nodeid(tree, labs)
   id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) 0.0 else if (is.na(freq[paste0(pre$incoming[[i]])])) 0.0 else freq[paste0(pre$incoming[[i]])])

   ldf <- data.frame(node = ids, frequency = id_freq, tip=tip)
   ldf$edge_id <- sapply(ldf$node, function(i) pre$incoming[[i]])
   #ldf$lab <- sapply(ldf$node, function (x) {
   #                 a <- MRCA_lab[which(MRCA_lab==x)]
   #                 if(length(a) > 0) return(a[1]) else return(NA)
   #            })
   tree.full <- full_join(tree, ldf, by = 'node')

   x_max <- -min(pre$nodes.df$times)

   p1 <- ggtree(tree.full, aes(color=frequency), size=1, ladderize=F) +
   geom_point() +
   #geom_text2(aes(label=edge_id, 
   #              subset=!is.na(lab), 
   #              x=branch), color="red", size=12, vjust=-1) +
   scale_size_manual(values=c(1)) +
   scale_x_continuous(limits=c(0, x_max)) +
   scale_color_viridis(option="plasma") +
   theme_tree2() + 
   theme_minimal() +
   theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          legend.position = "none")

}

p_mat[upper.tri(p_mat)]<-NA
p_df <- melt(p_mat)
names(p_df) <- c("sample_1", "sample_2", "value")
p_df$sample_1 <- factor(x = p_df$sample_1,
                               levels = pre$phy$tip.label,#[ord], 
                               ordered = TRUE)
p_df$sample_2 <- factor(x = p_df$sample_2,
                               levels = pre$phy$tip.label,#[ord], 
                               ordered = TRUE)

heatmap <- ggplot(data = p_df, aes(x = sample_1, y = sample_2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis_c(option= "plasma", na.value = "white") +
  labs(fill = "Pairwise Probability")+
  guides(fill=guide_legend(title.position = "right", vjust=0.5)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=30),
        legend.position = c(0.1,0.8),
        legend.title = element_text(angle = -90))

r_df <- as.data.frame(readResist())
r_df$x <- rownames(r_df) 
r_df <- melt(r_df)
r_df$x <- factor(x = r_df$x,
                     levels = pre$phy$tip.label,#[ord], 
                     ordered = TRUE)

resistmap <- ggplot(data = r_df, aes(x = x, y = variable)) +
                    geom_tile(aes(fill = factor(value, levels=c(0,1,2), labels = c("Susceptible", "Intermediate", "Resistant"), ordered=TRUE))) +
                    scale_fill_viridis(option= "viridis", na.value="gray50" , discrete=T) +
                    theme_minimal() +
                    ylab("Antimicrobial")+
                    labs(fill = "Resistance Level")+
                    guides(fill=guide_legend(title.position = "top"))+
                    theme(axis.title.x = element_blank(), 
                          axis.text.x = element_blank(), 
                          axis.ticks.x = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          text = element_text(size=30),
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          legend.position="bottom")

hist_dim <- ggplot(mcmc.df, aes(dim)) +  
   geom_histogram(aes(y = stat(count / sum(count))), binwidth=1) + 
   theme_minimal() +
   xlab("Number of Expansions")+
   ylim(0,1)+
   scale_fill_brewer(palette="Dark2")  + 
   theme(axis.title.y = element_blank(), 
         axis.text.y = element_blank(), 
         axis.ticks.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.margin = margin(1, 0, 0, 0, "cm"),
         text = element_text(size=30))

treemap <- plot_tree(pre, event.df) + scale_x_reverse()
summary_panel <- ggarrange(
        heatmap, treemap, resistmap,hist_dim,
        widths=c(16,12),heights=c(16,4))

png("grad2016_heatmap.png", width=32, height=24, "in", res=300, bg="white")
summary_panel
dev.off() 
