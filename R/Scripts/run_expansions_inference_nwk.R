#!/usr/bin/env Rscript
library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(rjson)
library(optparse)

option_list <- list(
   make_option(c("-n", "--niterations"), type="integer", default=NULL, 
       help="Number of iterations", metavar="integer"),
   make_option(c("-t", "--thinning"), type="integer", default=NULL, 
       help="Thinning", metavar="integer"),
   make_option(c("-f", "--filein"), type="character", default=NULL, 
       help="Path to .nwk file", metavar="character"),
   make_option(c("-s", "--seed"), type="integer", default=1, 
       help="(Optional) Random seed", metavar="integer"),
   make_option(c("-o", "--out"), type="character", default="mcmc_out", 
       help="(Optional) Output directory name [default= %default]", metavar="character"),
   make_option(c("--meanscale"), type="double", default=3, 
       help="(Optional) Background population size prior mean [default=%default]", metavar="double"),
   make_option(c("--sdscale"), type="double", default=3, 
       help="(Optional) Background population size prior sd [default=%default",  metavar="double"),
   make_option(c("--lambdar"), type="double", default=5, 
       help="(Optional) Growth rate / time to mid point prior lambda [default=%default",  metavar="double"),
   make_option(c("--nu"), type="double", default=1/2, 
       help="(Optional) Expansion time prior nu [default=%default",  metavar="double"),
   make_option(c("--kappa"), type="double", default=1/2, 
       help="(Optional) Expansion time prior kappa [default=%default",  metavar="double"),
   make_option(c("--sdk"), type="double", default=1, 
       help="(Optional) Expansion size prior sd [default=%default]", metavar="double")
   ) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$niterations)){
       n_it <- opt$niterations
} else {
       print_help(opt_parser)
       stop("Number of iterations must be supplied.", call.=FALSE)
}

if (!is.null(opt$thinning)) {
       thinning <- opt$thinning
} else {
       print_help(opt_parser)
       stop("Thinning must be supplied.", call.=FALSE)
}


if (!is.null(opt$filein)) {
       tree_in <- opt$filein
} else {
       print_help(opt_parser)
       stop("Simulation file must be supplied.", call.=FALSE)
}


if (!is.null(opt$out)) {
       dir_out <- opt$out
} else {
       print_help(opt_parser)
       stop("Output directory must be supplied.", call.=FALSE)
}

tree <- read.tree(file = tree_in)
tree <- makeNodeLabel(tree)

set.seed(opt$seed)
setwd(file.path(".", dir_out))

N_mean <- opt$meanscale ## carrying capacity rate lognormal prior mean

N_sd <- opt$sdscale ## carrying capacity  rate lognormal prior sd
K_sd <- opt$sdk

kappa <- opt$kappa
nu <- opt$nu

lambda_r <- opt$lambdar

priors <- standard_priors(expansion_rate=1, 
    N_mean_log=N_mean, 
    N_sd_log=N_sd, 
    t_mid_rate=lambda_r, 
    K_sd_log=K_sd, 
    exp_time_nu=nu,   
    exp_time_kappa=kappa)

expansions <- run_expansion_inference(tree, priors, 1, n_it=n_it, thinning=thinning)
saveRDS(expansions, "./expansions.rds")
warnings()