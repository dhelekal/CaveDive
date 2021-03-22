library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(gridExtra)
library(RColorBrewer)

set.seed(3)

run_mcmc <- F

data_dir <- "./tree_sim"
if(run_mcmc) {
    dir.create(file.path(data_dir))
}

tips <- 200
priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=4, 
                    N_sd_log=1/2, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/4)
given <- list(n_exp=3)

sam <- runif(200,-10, 0)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

sim <- expansions_simulate(priors, sam, 2, given=given)

co <- sim$co
params <- sim$params

set.seed(3)
phy.txt.div <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
phy.div <- read.tree(text=phy.txt.div$full)
phy <- collapse.singles(phy.div)

###Lets use less informative priors for inference
inference_priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=4, 
                    N_sd_log=4, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/2)

n_it <- 1e7
thinning <- n_it/1e4
burn_in <- 1e3 #10%
if(run_mcmc){
    start <- proc.time()
    expansions <- run_expansion_inference(phy, inference_priors, 1, n_it=n_it, thinning=thinning)
    elapsed <-proc.time() - start

    print(paste0(n_it," RjMCMC iterations completed. Time elapsed:"))
    print(elapsed)

    pre <- expansions$phylo_preprocessed ### keep a hold of this, it's useful for plotting
    dfs <- mcmc2data.frame(expansions$mcmc_out)
    mcmc.df <- dfs$mcmc.df
    event.df <- dfs$event.df
    write.csv(mcmc.df, paste0(data_dir,"/mcmc_df.csv"))
    write.csv(event.df, paste0(data_dir,"/event_df.csv"))
} else {
    mcmc.df <- read.csv(paste0(data_dir,"/mcmc_df.csv"))
    event.df <- read.csv(paste0(data_dir,"/event_df.csv"))
    pre <- structured_coal.preprocess_phylo(phy)
}

### Get expansion root edges exploiting the fact that the simulator names expansion nodes "N_[A-Z]".
root_set <- rep(NA, given$n_exp)
if(given$n_exp > 0) {
    for (i in c(1:given$n_exp)) {
        L <- LETTERS[i]
        N_set <- nodeid(pre$phy, pre$phy$node.label[grep(paste0("N_",L), pre$phy$node.label)]) 
        set_root <- N_set[which.min(pre$nodes.df$times[N_set])]
        set_root.edge <- pre$incoming[[set_root]]
        root_set[i] <- set_root.edge
    }
}

### discard burn in
mcmc.df <- mcmc.df[burn_in:nrow(mcmc.df),]
event.df <- event.df[which(event.df$it >= burn_in), ]

event.df$clade <- sapply(event.df$br, function (x) {
    a <- root_set[which(root_set==x)]
    if(length(a) > 0) return(a[1]) else return(NA)
})
event.df$is.gt <- !sapply(event.df$clade, is.na)

### marginalise mode model
unique_dims <- unique(mcmc.df$dim)
dim_counts <- sapply(unique_dims, function(i) length(which(mcmc.df$dim==i)))
mode_dim <- unique_dims[which.max(dim_counts)]

mode_dim_marginal <- mcmc.df[which(mcmc.df$dim==mode_dim),]
mode_dim_marginal_it <- mode_dim_marginal$it

event_mode_dim_marginal <- event.df[unlist(sapply(mode_dim_marginal_it, function (x) which(event.df$it==x))),]
unique.br <- unique(event_mode_dim_marginal$br)
branch_counts <- sapply(unique.br, function(x) length(which(event_mode_dim_marginal$br == x)))
br_count_ord <- order(-branch_counts)
k_mode_br <- unique.br[br_count_ord[c(1:mode_dim)]]

mode_br_df <- event_mode_dim_marginal[event_mode_dim_marginal$br %in% k_mode_br, ]
mode_br_mcmc_df <- mode_dim_marginal[mode_dim_marginal$it %in% mode_br_df$it, ]

hist_dim <- ggplot(mcmc.df, aes(dim)) +  
   geom_histogram(aes(y = stat(count / sum(count))), binwidth=1) + 
   theme_bw() +
   scale_fill_brewer(palette="Dark2")  + 
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
         text = element_text(size=20))
hist_N <- ggplot(mcmc.df, aes(N)) +
    geom_histogram(aes(y = stat(count / sum(count))), bins=100) +
    theme_bw() + 
    scale_fill_brewer(palette="Dark2")  + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          text = element_text(size=20))
hist_br <- ggplot(event.df, aes(x=factor(br), fill=is.gt)) +  
   geom_bar(aes(y = stat(count / sum(count)))) + 
   geom_text(stat="count", aes(label = clade, y= ((..count..)/sum(..count..))), vjust = -.25, hjust=-0.1, size=11, color="red")+
   theme_bw() + 
   scale_fill_brewer(palette="Dark2") + 
   labs(x="Branch Number",fill="Expansion Root") +
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = c(0.8, 0.2),
         text = element_text(size=20))

lab_layer <- geom_label(aes(x=branch, label=clade)) 

tree_freq <- plot_tree_freq(mcmc.df, event.df, pre, prior_t_given_N=NULL, highlight_node=NULL, MRCA_lab=pre$edges.df$node.child[root_set])

grid_layout <- rbind(c(1, 2), c(3,4))
grid_width <- c(2,2)
grid_heigth <- c(2,2)
summary_panel <- arrangeGrob(
        grobs=list(hist_N, hist_dim, hist_br, tree_freq),
        layout_matrix = grid_layout,
        widths = grid_width,
        heights = grid_heigth)

dummy_gt <- data.frame(br=unique(mode_br_df$br))
dummy_gt$gt.K <- sapply(dummy_gt$br, function (x) params$K[which(root_set==x)])
dummy_gt$gt.t_mid <- sapply(dummy_gt$br, function (x) params$t_mid[which(root_set==x)])

br.labs <- sapply(unique(mode_br_df$br), function (x) paste0("Branch: ",x))
names(br.labs) <- unique(mode_br_df$br)

prior_mixture <- function(prior, cond_values) {
     f_mixture <- function (X) sapply(X, function (x) (1/length(cond_values))*sum(sapply(cond_values, function(y) prior(x, y))))
     return(stat_function(fun=f_mixture, colour="purple", size=2))
}

K_facet <- ggplot(mode_br_df, aes(x=K)) + 
           geom_histogram(aes(y = stat(count / sum(count))), bins=100) +
           prior_mixture(function(x,N) exp(inference_priors$prior_K_given_N(x,N)),mode_br_mcmc_df$N) +
           facet_wrap(~br, labeller=labeller(br = br.labs)) +
           geom_vline(data = dummy_gt, aes(xintercept = gt.K), color="red",lwd=3) +
           labs(x="Carrying Capacity") +
           theme_bw() +
           theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),text = element_text(size=20))

t_mid_facet <- ggplot(mode_br_df, aes(x=t_mid)) + 
           geom_histogram(aes(y = stat(count / sum(count))), bins=100) +
           prior_mixture(function(x,N) exp(inference_priors$prior_t_mid_given_N(x,N)),mode_br_mcmc_df$N) +
           facet_wrap(~br, labeller=labeller(br = br.labs)) +
           geom_vline(data = dummy_gt, aes(xintercept = gt.t_mid), color="red",lwd=3) +
           labs(x="Time to Midpoint") +
           theme_bw() +
           theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),text = element_text(size=20))

event_panel <- arrangeGrob(
        grobs=list(K_facet, t_mid_facet),
        nrow=2,
        heights = c(1,1))
png("fig2a.png", width=1800, height=1600)
plot(event_panel)
dev.off()

png("fig2b.png", width=1600, height=1600)
plot(summary_panel)
dev.off()

