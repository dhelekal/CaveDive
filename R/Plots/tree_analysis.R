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

run_mcmc <- FALSE

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
phy.txt.div <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from, include_div_nodes=TRUE)
phy.div <- read.tree(text=phy.txt.div$full)
phy <- collapse.singles(phy.div)

tree_plt <- plot_structured_tree(phy.div, 4) + scale_color_brewer(palette="Dark2") +
            theme(legend.position="right",
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank())
png("tree2.png",1600,1600)
plot(tree_plt)
dev.off()

tree_plt <- plot_structured_tree(phy, 4) + scale_color_brewer(palette="Dark2") +
            theme(legend.position="right",
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank())
png("tree3.png",1600,1600)
plot(tree_plt)
dev.off()


###Lets use less informative priors for inference
inference_priors <- standard_priors(expansion_rate=1, 
                    N_mean_log=4, 
                    N_sd_log=4, 
                    t_mid_rate=5, 
                    K_sd_log=1, 
                    exp_time_nu=1/2, 
                    exp_time_kappa=1/2)

n_it <- 1e7
thinning <- 1e7/1e4
if(run_mcmc){
    phy <- makeNodeLabel(phy)
    start <- proc.time()
    expansions <- run_expansion_inference(phy, inference_priors, 1, n_it=n_it, thinning=thinning)
    elapsed <-proc.time() - start

    print(paste0(n_it," RjMCMC iterations completed. Time elapsed:"))
    print(elapsed)

    pre <- expansions$phylo_preprocessed ### keep a hold of this, it's useful for plotting
    dfs <- mcmc2data.frame(expansions$mcmc_out)
    mcmc.df <- dfs$mcmc.df
    event.df <- dfs$event.df
    write.csv(mcmc.df, paste0("./mcmc_df.csv"))
    write.csv(event.df, paste0("./event_df.csv"))
} else {
    mcmc.df <- read.csv("./mcmc_df.csv")
    event.df <- read.csv("./event_df.csv")
    pre <- structured_coal.preprocess_phylo(phy)
}

png("dim_panel.png",1600,1600)
plot(plot_dim_panel(mcmc.df, prior_N=NULL))
dev.off()

png("tree_freq.png",1600,1600)
plot(plot_tree_freq(mcmc.df, event.df, pre, prior_t_given_N=NULL, highlight_node=NULL))
dev.off()
