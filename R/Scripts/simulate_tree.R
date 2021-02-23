#!/usr/bin/env Rscript
library(DirichletReg)
library(rjson)
library(optparse)

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
    make_option(c("-A", "--rates"), type="character", default=NULL, 
              help="(Optional) Comma separated list of non-dimensional times to midpoint, A = 1/sqrt(r) with r being the non-dimensional
                    rate [default=randomise]",  metavar="character"),
    make_option(c("-S", "--samscale"), type="double", default=-1, 
              help="(Optional) Sampling scale lower bound (upper bound always 0) [default=%default]", metavar="double"),
    make_option(c("-l", "--lambdar"), type="double", default=10, 
              help="(Optional) Lambda_r value for non-dimensional time to midpoint distribution [default=%default]", metavar="double"),
    make_option(c("--metadata"), type="double", default=NULL, 
              help="(Optional) Metadata to be included in the simulation file in a comma separated list
                with entries in the format of [name]:[value] pairs [default=%default]", metavar="double")
) 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

### Required arguments

if (is.null(opt$nexp)){
  print_help(opt_parser)
  stop("Number of expansions must be supplied.", call.=FALSE)
} else {
    n_exp <- opt$nexp
}
if (is.null(opt$tips)){
  print_help(opt_parser)
  stop("Number of tips must be supplied.", call.=FALSE)
} else {
    n_tips <- opt$tips
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

if (is.null(opt$popscale)){
    pop_scale <- rlnorm(1,3,3)
} else {
    pop_scale <- opt$popscale
}


colouring <- c()
exp_probs <- c()
clade_sizes <- c()

sam <- runif(n_tips, sampling_scale, 0)
sam <- sam[order(-sam)]
sam <- sam - max(sam)

concentration <- c(rep(2,(n_exp)),2)

lambda_r <- opt$lambdar
kappa <- 1/10
nu <- 1/2
sigma_k <- 1

prior_t_given_N.sample <- function(k, N) (-rgamma(k, shape=(nu^2)/kappa, scale = kappa * N / nu))

max_it <- 50
it <- 0

while((length(unique(colouring)) < (n_exp+1)) || (!all(clade_sizes>1) && n_exp > 0)) {
    if (it > max_it) {
        stop("Maximum sampling iterations exceeded.")
    }
    exp_probs <- rdirichlet(1, concentration)
    colouring <- rmultinom(n_tips, 1, exp_probs)
    colouring <- sapply(c(1:n_tips), function (i) which(colouring[,i]>0))
    clade_sizes <- sapply(c(1:(n_exp+1)), function(i) length(which(colouring==i)))
    it <- it + 1
}


if (n_exp > 0) {

    if (is.null(opt$rates)){
        t_mid <- rexp(n_exp, lambda_r/pop_scale)
    } else {
        t_mid <- as.double(unlist(strsplit(opt$rates, ",")))
        t_mid <- t_mid*pop_scale
    }

    if (is.null(opt$capacities)){
        K <- rlnorm(n_exp, log(pop_scale), sigma_k)
    } else {
        K <- as.double(unlist(strsplit(opt$capacities, ",")))
    }

    div.times <- rep(Inf, n_exp+1)

    max_it <- 50
    it <- 0

    while(!all(sapply(c(1:(n_exp+1)), function(i) all(div.times[i] < sam[which(colouring==i)])))){
        if (it > max_it) {
            stop("Maximum sampling iterations exceeded.")
        }
        ## div time now
        div.times <- prior_t_given_N.sample(n_exp, pop_scale)

        div.times <- div.times
        div.times <- c(div.times,-Inf)

        ## Re-order everything so that divergence event numbering corresponds to their order of occurence
        div.times.ord <- order(-div.times) 
        div.times <- div.times[div.times.ord]

        colouring <- sapply(colouring, function(i) which(div.times.ord==i))
        it <- it+1
    }    
} else {
    t_mid <- c()
    K <- c()
    div.times <- c(-Inf)
}


out <- list(seed=seed, n_tips=n_tips, n_exp=n_exp, exp_probs=exp_probs, tip_times=sam, tip_colours=colouring, N=pop_scale, K=K, t_mid=t_mid, div_times=div.times, meta=meta)
out_txt <- toJSON(out)

write(out_txt, paste0("./","tree_params.json"))