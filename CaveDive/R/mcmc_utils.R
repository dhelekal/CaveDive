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
   prob.neutral <- sapply(o$para, function(x) x[[2]][1])

   it <- c(1:length(dim))
   mcmc.df <- data.frame(it=it, dim=dim, N=N, pr=prob.neutral, lh=lh, prior=prior)

   events <- lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(c(1:length(o$para[[x]][c(-1, -2)])), function(i){
     y <- (o$para[[x]][c(-1, -2)])[[i]]
     probs <- o$para[[x]][[2]][-1]
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

#' Plots mcmc output as several panels
#' 
#' @param mcmc.df mcmc.df returned by mcmc2data.frame
#' @param event.df event.df returned by mcmc2data.frame
#' @param which_br Which branch to generate marginals for
#' @param pre Preprocessed phylogeny
#' @param prior_N (Optional) Background population size prior, mutually exclusive with passing list of priors. If either is supplied priors will be overlayed in plotting.
#' @param prior_t_mid_given_N (Optional) Time to midpoint  prior, mutually exclusive with passing list of priors. If either is supplied priors will be overlayed in plotting.
#' @param prior_K_given_N (Optional) Carrying capacity prior, mutually exclusive with passing list of priors. If either is supplied priors will be overlayed in plotting.
#' @param prior_K_given_N (Optional) Expansion time prior, mutually exclusive with passing list of priors. If either is supplied priors will be overlayed in plotting.
#' @return a list of 3 plot panels
#' @export
plot_event_summary <- function(mcmc.df, event.df, which_br, pre, 
   prior_N=NULL, 
   prior_t_mid_given_N=NULL, 
   prior_K_given_N=NULL, 
   prior_t_given_N=NULL) 
{

   correct_dim <- mcmc.df[which(mcmc.df$dim > 0),]
   correct_dim_it <- correct_dim$it
   event_dim_marginal <- event.df[unlist(sapply(correct_dim_it, function (x) which(event.df$it==x))),]
   br_subs <- which(event_dim_marginal$br == which_br)
   event_br_marginal <- event_dim_marginal[br_subs,]
   event_br_marginal$it <- correct_dim_it[br_subs]
   event_br_marginal$idx <- if(length(event_br_marginal$it) > 0) c(1:length(event_br_marginal$it)) else c()

   dim_panel <- plot_dim_panel(mcmc.df, prior_N)

   trace_K <- ggplot(event_br_marginal, aes(x=idx, y=K)) +
   geom_line(alpha = 0.3)+
   theme_bw() +
   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
   trace_t_mid <- ggplot(event_br_marginal, aes(x=idx, y=t_mid)) +
   geom_line(alpha = 0.3)+
   theme_bw() +
   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
   hist_K <- ggplot(event_br_marginal, aes(K)) +  
   geom_histogram(aes(y = stat(density)), colour="orange", fill="orange", bins=100) +
   theme_bw() +
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
   hist_t_mid <- ggplot(event_br_marginal, aes(t_mid)) +  
   geom_histogram(aes(y = stat(density)), colour="orange", fill="orange", bins=100) + 
   theme_bw() +
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

   if (!is.null(prior_K_given_N)) {
     hist_K <- hist_K + prior_mixture(prior_K_given_N, mcmc.df$N)
   }

   if (!is.null(prior_t_mid_given_N)) {
      hist_t_mid <- hist_t_mid + prior_mixture(prior_t_mid_given_N, mcmc.df$N)
   }

   event_panel <- arrangeGrob(
        grobs=list(hist_K, hist_t_mid, trace_K, trace_t_mid),
        layout_matrix = grid_layout,
        widths = grid_width,
        heights = grid_heigth)

   trace_br <- ggplot(event.df, aes(x=it, y=br)) + 
   geom_point(alpha=0.1,size=0.1) + 
   theme_bw() + 
   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
   tree_freq <- plot_tree_freq(mcmc.df, event.df, pre, prior_t_given_N=prior_t_given_N)

   hist_br <- ggplot(event.df, aes(br)) +  
   geom_bar(aes(y = stat(count / sum(count))), colour="orange", fill="orange") + 
   theme_bw() + 
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

   grid_layout <- rbind(c(1,2), c(3,3))
   grid_width <- c(3,3)
   grid_heigth <- c(3,2)

   tree_panel <- arrangeGrob(
        grobs=list(tree_freq, hist_br, trace_br),
        layout_matrix = grid_layout,
        widths = grid_width,
        heights = grid_heigth)
   return(list(dim_panel=dim_panel, event_panel=event_panel, tree_panel=tree_panel))
}


#' Plots a panel summarising the inferred number of expansions and base population size
#' 
#' @param mcmc.df mcmc.df returned by mcmc2data.frame
#' @param prior_N (Optional) Background population size prior, mutually exclusive with passing list of priors. If either is supplied priors will be overlayed in plotting.
#' @return a panel containing histograms and traces for the inferred number of expansions and the base population size
#' @export
plot_dim_panel <- function(mcmc.df, prior_N=NULL) {
   trace_N <- ggplot(mcmc.df, aes(x=it, y=N)) +
   geom_line(alpha = 0.3) +
   theme_bw() + 
   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
   trace_dim <- ggplot(mcmc.df, aes(x=it, y=dim)) +
   geom_line(alpha = 0.3) +
   theme_bw() +
   theme(axis.title.x = element_blank(), axis.text.x = element_blank())
   hist_N <- ggplot(mcmc.df, aes(N)) +  
   geom_histogram(aes(y = stat(density)), colour="orange", fill="orange", bins=100) +
   theme_bw() + 
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
   hist_dim <- ggplot(mcmc.df, aes(dim)) +  
   geom_histogram(aes(y = stat(density)), colour="orange", fill="orange", binwidth=1) + 
   theme_bw() +
   theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
   
   if (!is.null(prior_N)) {
     hist_N <- hist_N + stat_function(fun=prior_N, colour="purple", size=2)
   }

   grid_layout <- rbind(c(1, 2), c(3,3), c(4,4))
   grid_width <- c(2,2)
   grid_heigth <- c(2,1,1)

   dim_panel <- arrangeGrob(
        grobs=list(hist_N, hist_dim, trace_N, trace_dim),
        layout_matrix = grid_layout,
        widths = grid_width,
        heights = grid_heigth)
   return(dim_panel)
}


#' Plots tree branch frequency histogram 
#' 
#' @param mcmc.df mcmc.df returned by mcmc2data.frame
#' @param event.df event.df returned by mcmc2data.frame
#' @param pre Preprocessed phylogeny
#' @param prior_t_given_N (Optional) Expansion time prior. If supplied prior will be overlayed.
#' @return a plot object
#' @export
plot_tree_freq <- function(mcmc.df, event.df, pre, prior_t_given_N=NULL) {
   tree <- pre$phy
   freq <- table(event.df$br)

   labs <- c(tree$node.label, tree$tip.label)
   tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
   ids <- nodeid(tree, labs)
   id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) NA else if (is.na(freq[paste0(pre$incoming[[i]])])) 0.0 else freq[paste0(pre$incoming[[i]])])

   ldf <- data.frame(node = ids, frequency = id_freq, tip=tip)
   tree.full <- full_join(tree, ldf, by = 'node')

   x_max <- -min(pre$nodes.df$times)

   p1 <- ggtree(tree.full, aes(color=frequency, linetype=tip), ladderize=TRUE) +
   geom_point() +
   scale_linetype(c("solid","dashed","dotted"), na.value = "blank") +
   scale_size_manual(values=c(1)) +
   scale_x_continuous(limits=c(0, x_max)) +
   scale_color_viridis() +
   theme_tree2() + 
   theme_bw() +
   theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


   temp <- ggplotGrob(p1)
   leg_index <- which(sapply(temp$grobs, function(x) x$name) == "guide-box")
   legend <- temp$grobs[[leg_index]]
   p1 <- p1 + theme(legend.position="none")

   event.df_invtime <- event.df
   event.df_invtime$time <- x_max+event.df_invtime$time

   p2 <- ggplot(event.df_invtime, aes(time)) +  
   geom_histogram(aes(y = stat(density)), colour="orange", fill="orange", breaks=seq(0, x_max,length.out=100)) +
   scale_x_continuous(limits=c(0, x_max))
   if (!is.null(prior_t_given_N)) {
     p2 <- p2 + prior_mixture(function(t,n) prior_t_given_N(t-x_max, n), mcmc.df$N)
   }
   p2 <- p2 + theme_bw() + theme(axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank())

   grid_layout <- rbind(c(1,3), c(2,NA))
   grid_width <- c(5,1)
   grid_heigth <- c(5,1)

   p <-arrangeGrob(
        grobs=list(p1, p2,legend),
        layout_matrix = grid_layout,
        widths = grid_width,
        heights = grid_heigth)
   return(p)
}

prior_mixture <- function(prior, cond_values) {
     f_mixture <- function (X) sapply(X, function (x) (1/length(cond_values))*sum(sapply(cond_values, function(y) prior(x, y))))
     return(stat_function(fun=f_mixture, colour="purple", size=2))
}