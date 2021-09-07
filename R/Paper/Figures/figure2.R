library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(coda)

set.seed(3)

run_mcmc <- F

data_dir <- "./tree_sim"

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

if(run_mcmc) {
    dir.create(file.path(data_dir))
    set.seed(3)
    phy.txt.div <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
    phy.div <- read.tree(text=phy.txt.div$full)
    phy <- collapse.singles(phy.div)

    n_it <- 1e7
    thinning <- n_it/1e4

    set.seed(1)
    start <- proc.time()
    expansions <- run_expansion_inference(phy, priors, 2, n_it=n_it, thinning=thinning)
    elapsed <-proc.time() - start

    print(paste0(n_it," RjMCMC iterations completed. Time elapsed:"))
    print(elapsed)

    sink("TimeFig2.txt")
    print(paste0(n_it," RjMCMC iterations completed. Time elapsed:"))
    print(elapsed)
    sink()

    expansions <- discard_burn_in(expansions, proportion=0.1)
    
    saveRDS(expansions, file = paste0(data_dir, "/expansions.rds"))

    set.seed(2)
    start <- proc.time()
    expansions2 <- run_expansion_inference(phy, priors, 2, n_it=n_it, thinning=thinning)
    elapsed <-proc.time() - start

    print(paste0(n_it," RjMCMC iterations completed. Time elapsed:"))
    print(elapsed)

    expansions2 <- discard_burn_in(expansions2, proportion=0.1)
    
    saveRDS(expansions2, file = paste0(data_dir, "/expansions2.rds"))
} else {
    expansions <- readRDS(file = paste0(data_dir, "/expansions.rds"))
    expansions2 <- readRDS(file = paste0(data_dir, "/expansions2.rds"))
}

png("fig2_pop_fn.png", width=1600, height=800)
plot(expansions, mode="popfunc", k_modes=3, gt.K=params$K[c(1,3,2)], gt.t_mid=params$t_mid[c(1,3,2)], gt.time=params$div_times[c(1,3,2)], t_max=c(20,20,20), tree_scale="Years")

dev.off()

png("fig2a.png", width=1600, height=1600)
plot(expansions, mode="summary", k_modes=3)
dev.off()

png("fig2b.png", width=1800, height=1600)
plot(expansions, mode="modes", k_modes=3, gt.K=params$K[c(1,3,2)], gt.t_mid=params$t_mid[c(1,3,2)])
dev.off()

png("fig2_trace.png", width=1600, height=1600)
plot(expansions, mode="traces")
dev.off()

png("fig2_mode_trace.png", width=1600, height=1600)
plot(expansions, mode="mtraces", k_modes=3)
dev.off()

m1 <- mcmc(expansions$model_data[,2:3],thin=expansions$metadata$thinning)
m2 <- mcmc(expansions2$model_data[,2:3],thin=expansions2$metadata$thinning)

grs <- gelman.diag(mcmc.list(m1,m2))
sink("GelmanRubinFig2.txt")
print(grs$psrf)
print(grs$mpsrf)
print(effectiveSize(m1))
sink()