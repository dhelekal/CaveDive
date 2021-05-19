priorList <- function(prior_i,
                    prior_i.sample,
                    prior_N,
                    prior_N.sample,
                    prior_t_mid_given_N,
                    prior_t_mid_given_N.sample,
                    prior_K_given_N,
                    prior_K_given_N.sample,
                    prior_t_given_N,
                    prior_t_given_N.sample) {
     return(structure(list(prior_i=prior_i, 
                            prior_i.sample=prior_i.sample,
                            prior_N=prior_N,
                            prior_N.sample=prior_N.sample,
                            prior_t_mid_given_N=prior_t_mid_given_N,
                            prior_t_mid_given_N.sample=prior_t_mid_given_N.sample,
                            prior_K_given_N=prior_K_given_N,
                            prior_K_given_N.sample=prior_K_given_N.sample,
                            prior_t_given_N=prior_t_given_N,
                            prior_t_given_N.sample=prior_t_given_N.sample), class="priorList"))
}

#' Constructs object of class expansionsMCMC
#' 
#' @param phylo_preprocessed an object of type preprocessedPhy
#' @param priors an object of type priorList
#' @param model_data dataframe containing model data as returned by mcmc2data.frame
#' @param expansion_data dataframe containing expansion data as returned by mcmc2data.frame
#' @param metadata mcmc metadata
#' @return an object of type expansionsMCMC
#' @export
expansionsMCMC <- function(phylo_preprocessed, priors, model_data, expansion_data, metadata) {
   stopifnot("pre must be of type preprocessed phylogeny"= class(phylo_preprocessed) == "preprocessedPhy")
   stopifnot("priors must be of type priorList"=class(priors)== "priorList")
   stopifnot("invalid model data"=all(colnames(model_data)==c("it", "dim", "N", "pr", "lh", "prior")))
   stopifnot("invalid expansion_data data"=all(colnames(expansion_data)==c("it", "t_mid", "K", "time", "br", "pr")))
   out <- list(phylo_preprocessed=phylo_preprocessed, priors=priors, model_data=model_data, expansion_data=expansion_data, metadata=metadata)
   attr(out, "class") <- "expansionsMCMC"
   return(out)
}

#' Converts MCMC output into two data frames, one conaining global model parameters and one containing expansion data
#' 
#' @param o MCMC output
#' @return A list with names mcmc.df and event.df. MCMC dataframe contains iteration numbers, 
#'         model indicator at given iteration, background population size and likelihood and prior values
#'         event.df contains rows describing different expansions. Columns consist of which iteration does an expansion belong to, 
#'         and the associated t_mid/K/time/branch/probability values.
#' @export
mcmc2data.frame <- function(o) {

   dim <- sapply(o$dims, function(x) x)
   N <- sapply(o$para, function(x) x[[1]])
   lh <- sapply(o$log_lh, function(x) x)
   prior <- sapply(o$log_prior, function(x) x)
   prob.neutral <- sapply(o$para, function(x) x[[2]][length(x[[2]])])

   it <- c(1:length(dim))
   mcmc.df <- data.frame(it=it, dim=dim, N=N, pr=prob.neutral, lh=lh, prior=prior)

   events <- lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(c(1:length(o$para[[x]][c(-1, -2)])), function(i){
     y <- (o$para[[x]][c(-1, -2)])[[i]]
     probs <- o$para[[x]][[2]][-length(o$para[[x]][[2]])]
     return(list(t_mid=y[[1]], K=y[[2]], time=y[[3]], br=y[[4]], pr=probs[i], it=x))}) else NA)
   events.t_mid <- unlist(lapply(events, function(x) if(all(is.na(x))) NA else lapply(x, function(y) y$t_mid)))
   events.K <- unlist(lapply(events, function(x) if(all(is.na(x))) NA else lapply(x, function(y) y$K)))
   events.time <- unlist(lapply(events, function(x) if(all(is.na(x))) NA else lapply(x, function(y) y$time)))
   events.br <- unlist(lapply(events, function(x) if(all(is.na(x))) NA else lapply(x, function(y) y$br)))
   events.pr <- unlist(lapply(events, function(x) if(all(is.na(x))) NA else lapply(x, function(y) y$pr)))
   events.it <-  unlist(lapply(events, function(x) if(all(is.na(x))) NA else lapply(x, function(y) y$it)))
   event.df <- data.frame(it=events.it, t_mid=events.t_mid, K=events.K, time=events.time, br=events.br, pr=events.pr)

   return(list(mcmc.df=mcmc.df, event.df=event.df))
}

#' @export
discard_burn_in <- function(x, ...) UseMethod("discard_burn_in")

#' Discard burn in iterations
#' @param x expansionsMCMC object
#' @param proportion (Optional) The proportion of iterations to be discarded. Either 'proportion' or 'k_it' must be provided 
#' @param k_it (Optional) The number of iterations to be discarded. Either 'proportion' or 'k_it' must be provided 
#' @return expansionsMCMC object with iterations designated as burn in discarded
#' @export
discard_burn_in.expansionsMCMC <- function(x, proportion=NULL, k_it=NULL, ...){
     stopifnot("x must be of type expansionsMCMC"=class(x)=="expansionsMCMC")
     stopifnot("Only one of arguments proportion, k_it must be provided"=xor(is.null(proportion), is.null(k_it)))

     metadata <- x$metadata
     if (!is.null(proportion)){
               stopifnot("proportion must be within interval [0,1)"=(proportion>=0 && proportion < 1))
               burnin <- metadata$effective_it*proportion
     } else {
          stopifnot("k_it must be within interval [0, effective_it)"=(k_it>=0 && k_it < effective_it))
          burnin <- k_it
     }

     model_data_burn_in <- x$model_data[(burnin+1):metadata$effective_it, ] 
     expansion_data_burn_in <- x$expansion_data[x$expansion_data$it %in% model_data_burn_in$it, ]

     metadata$effective_it <- metadata$effective_it-burnin 
     metadata$burn_in <- burnin 

     return(expansionsMCMC(x$phylo_preprocessed, x$priors, model_data_burn_in, expansion_data_burn_in, metadata))
}

#' @export
print.expansionsMCMC <- function(x, ...) {
     cat(paste("\nClonal Expansions MCMC result\n\n"))
     cat(paste("\nFields:", names(x)))
     cat(paste("\nNumber of mcmc iterations: ", x$metadata$n_it))
     cat(paste("\nThinning applied: ", x$metadata$thinning,"\n"))
}

#' @export
plot.expansionsMCMC <- function(x, ..., mode=c("summary", "modes", "persistence", "traces", "mtraces"), k_modes=NULL, correlates=NULL, corr_axis_title="Variable",
                                  corr_legend_title="Value", gt.K=NULL, gt.t_mid=NULL) {
     mode <- match.arg(mode)

     expansion_data <- x$expansion_data
     model_data <- x$model_data
     if (nrow(expansion_data > 0)) {
          expansion_data$mode_clade <- NA
          expansion_data$is.mode <- NA
     }

     modes <- NULL

     if(!is.null(k_modes)) {
          unique_br <- unique(expansion_data$br)

          stopifnot("Cannot select modes - no expansions were detected"=(length(unique_br) > 0))

          br_counts <- sapply(unique_br, function (x) length(which(expansion_data$br==x)))

          k <- min(length(unique_br), k_modes)
          mode_ord <- order(-br_counts)

          modes <- unique_br[mode_ord][c(1:k)]

          expansion_data$mode_clade <- sapply(expansion_data$br, function (x) {
              a <- modes[which(modes==x)]
              if(length(a) > 0) return(a[1]) else return(NA)
          })
          expansion_data$is.mode <- !sapply(expansion_data$mode_clade, is.na)
     }

     if (mode=="summary") {
          stopifnot("No expansions were detected"=nrow(expansion_data)>0)
          if (!is.null(correlates)) warning("Unused argument: correlates")
          plot_summary(model_data, expansion_data, x$phylo_preprocessed, x$priors, modes)

     } else if (mode == "modes") {
          stopifnot("No expansions were detected"=nrow(expansion_data)>0)
          stopifnot("Number of modes must be supplied"=!is.null(k_modes))
          if (!is.null(correlates)) warning("Unused argument: correlates")
          plot_mode_summary(model_data, expansion_data, x$priors, k_modes, gt.K, gt.t_mid)

     } else if (mode=="persistence") {
          stopifnot("No expansions were detected"=nrow(expansion_data)>0)
          plot_persistence(model_data,
                                  expansion_data, 
                                  x$phylo_preprocessed, 
                                  corr_axis_title,
                                  corr_legend_title,
                                  correlates=correlates,
                                  modes=modes)

     } else if (mode=="traces") {
          if (!is.null(k_modes)) warning("Unused argument: k_modes")
          if (!is.null(correlates)) warning("Unused argument: correlates")
          plot_traces(model_data, expansion_data)
     } else if (mode=="mtraces") {
          stopifnot("No expansions were detected"=nrow(expansion_data)>0)
          stopifnot("Number of modes must be supplied"=!is.null(k_modes))
          plot_mode_traces(model_data, expansion_data, k_modes)
     } else {
          stop("Invalid plotting options")
     }
}