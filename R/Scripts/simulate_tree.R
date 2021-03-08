#!/usr/bin/env Rscript
library(DirichletReg)
library(rjson)
library(optparse)
library(CaveDive)
library(ape)
library(treeio)

option_list <- list(
    make_option(c("-e", "--nexp"), type="integer", default=NULL, 
      help="Number of expansions", metavar="integer"),
    make_option(c("-t", "--tips"), type="integer", default=NULL, 
      help="Number of tips", metavar="integer"),
    make_option(c("-s", "--seed"), type="integer", default=c(1), 
      help="(Optional) RNG seed [default= %default]", metavar="integer"),
    make_option(c("-o", "--out"), type="character", default="simulation", 
      help="(Optional) Output folder name [default= %default]", metavar="character"),
    make_option(c("-N", "--popscale"), type="double", default=NULL, 
      help="(Optional) Background population size [default=randomise]", metavar="double"),
    make_option(c("-K", "--capacities"), type="character", default=NULL, 
      help="(Optional) Comma separated list of carrying capacities [default=randomise]",  metavar="character"),
    make_option(c("-T", "--t_mid"), type="character", default=NULL, 
      help="(Optional) Comma separated list of non-dimensional times to midpoint, A = 1/sqrt(r) with r being the non-dimensional
      rate [default=randomise]",  metavar="character"),
    make_option(c("-S", "--samscale"), type="double", default=-1, 
      help="(Optional) Sampling scale lower bound (upper bound always 0) [default=%default]", metavar="double"),
    make_option(c("-l", "--lambdar"), type="double", default=10, 
      help="(Optional) Lambda_r value for non-dimensional time to midpoint distribution [default=%default]", metavar="double"),
    make_option(c("--nu"), type="double", default=1/2, 
       help="(Optional) Expansion time prior nu [default=%default",  metavar="double"),
    make_option(c("--kappa"), type="double", default=1/8, 
       help="(Optional) Expansion time prior kappa [default=%default",  metavar="double"),
    make_option(c("--sdk"), type="double", default=1, 
       help="(Optional) Expansion size prior sd [default=%default]", metavar="double"),
    make_option(c("--meanscale"), type="double", default=3, 
       help="(Optional) Background population size prior mean [default=%default]", metavar="double"),
    make_option(c("--sdscale"), type="double", default=3, 
       help="(Optional) Background population size prior sd [default=%default",  metavar="double"),
    make_option(c("--metadata"), type="character", default=NULL, 
      help="(Optional) Metadata to be included in the simulation file in a comma separated list
      with entries in the format of [name]:[value] pairs [default=%default]", metavar="double")
    ) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
given <- list()

### Required arguments

if (is.null(opt$nexp)){
  print_help(opt_parser)
  stop("Number of expansions must be supplied.", call.=FALSE)
} else {
    given$n_exp <- opt$nexp
}
if (is.null(opt$tips)){
  print_help(opt_parser)
  stop("Number of tips must be supplied.", call.=FALSE)
} else {
  given$n_tips <- opt$tips
}
if (is.null(opt$metadata)){
    meta <- NA
} else {
    meta_entries <- unlist(strsplit(opt$metadata, ","))
    meta_names <- sapply(c(1:length(meta_entries)), function (x) unlist(strsplit(meta_entries[x], ":"))[1])
    meta_vals <- sapply(c(1:length(meta_entries)), function (x) unlist(strsplit(meta_entries[x], ":"))[2])
    meta <- as.list(meta_vals)
    names(meta) <- meta_names 
}

### Optional Arguments
dir_out <- opt$out
seed <- opt$seed
sampling_scale <- opt$samscale

set.seed(seed)
dir.create(file.path(".", dir_out))
setwd(file.path(".", dir_out))

given$N <- opt$popscale  

sam <- runif(given$n_tips, sampling_scale, 0)
sam <- sam[order(-sam)]
sam <- sam - max(sam)

concentration <- 2

lambda_r <- opt$lambdar
sigma_k <- 1

if (given$n_exp > 0) {
    if (is.null(opt$t_mid)){
        given$t_mid <- NULL
    } else {
        t_mid <- as.double(unlist(strsplit(opt$t_mid, ",")))
        given$t_mid <- t_mid*given$pop_scale
    }

    if (is.null(opt$capacities)){
        given$K <- NULL
    } else {
        given$K <- as.double(unlist(strsplit(opt$capacities, ",")))
    }
}

priors <- standard_priors(expansion_rate=1, 
    N_mean_log=opt$meanscale, 
    N_sd_log=opt$sdscale, 
    t_mid_rate=lambda_r, 
    K_sd_log=opt$sdk, 
    exp_time_nu=opt$nu,   
    exp_time_kappa=opt$kappa)
out <- expansions_simulate(priors, sam, concentration, given=given)
params <- out$params
co <- out$co
print(co)
print(params)

phy <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from, include_div_nodes=FALSE)
phy.div_nodes <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from, include_div_nodes=TRUE)

tree  <- read.tree(text = phy$full)
pre <- structured_coal.preprocess_phylo(tree)

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

sim_data <- params
sim_data$seed <- seed
sim_data$tip_times <- sam
sim_data$root_set <- root_set
sim_data$meta <- meta
sim_data_txt <- toJSON(sim_data)

tree.div  <- read.tree(text = phy.div_nodes$full)
pdf(file="tree.pdf", width=5, height=5)
tree.plt <- plot_structured_tree(tree.div, given$n_exp+1)
plot(tree.plt)
dev.off()

write(sim_data_txt, paste0("./","tree_params.json"))
write(phy$full, "./tree.nwk")
write(phy.div_nodes$full, "./tree_div.nwk")