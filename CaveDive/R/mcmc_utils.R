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

#' Converts MCMC output into two data frames, one conaining global model parameters and one containing expansion data
#' 
#' @param o MCMC output
#' @return a plot object
#' @export
plot_event_summary <- function(mcmc.df, event.df, which_br, prior_N=NULL, 
                                                            prior_tmid_given_N=NULL, 
                                                            prior_K_given_N=NULL, 
                                                            prior_t_given_N=NULL,
                                                            pre=NULL) {

    correct_dim <- mcmc_data[which(mcmc_data$dim > 0),]
    correct_dim_it <- correct_dim$it
    event_dim_marginal <- chain_data[unlist(sapply(correct_dim_it, function (x) which(chain_data$it==x))),]
    mode_subs <- which(event_dim_marginal$br == root_set[[1]])
    event_br_marginal <- event_dim_marginal[mode_subs,]
    event_br_marginal$it <- correct_dim_it[mode_subs]
    event_br_marginal$idx <- c(1:length(event_br_marginal$it))

        png(file="trace_K.png", width=800, height=800)
        plt <- ggplot(event_br_marginal, aes(x=idx, y=K)) +
               geom_line(alpha = 0.3)+
               theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
        plot(plt)
    dev.off()

    png(file="trace_A.png", width=800, height=800)
        plt <- ggplot(event_br_marginal, aes(x=idx, y=t_mid)) +
               geom_line(alpha = 0.3)+
               theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
        plot(plt)
    dev.off()
    

}

#' Plots tree branch frequency histogram 
#' 
#' @param mcmc.df mcmc.df returned by mcmc2data.frame
#' @param event.df event.df returned by mcmc2data.frame
#' @param pre preprocessed phylogeny
#' @return a plot object
#' @export
plot_tree_freq <- function(mcmc.df, event.df, pre) {
    tree <- pre$phy
    tt.br <- table(event.df$br)
    freq <- sapply(c(1:length(tt.br)), function (i) tt.br[i])

    labs <- c(tree$node.label, tree$tip.label)
    tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
    ids <- nodeid(tree, labs)
    id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) NA else if (is.na(freq[paste0(pre$incoming[[i]])])) 0 else freq[paste0(pre$incoming[[i]])])
    
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
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())

    temp <- ggplotGrob(p1)
    leg_index <- which(sapply(temp$grobs, function(x) x$name) == "guide-box")
    legend <- temp$grobs[[leg_index]]
    p1 <- p1 + theme(legend.position="none")

    event.df_invtime <- event.df
    event.df_invtime$time <- x_max+event.df_invtime$time

    p2 <- ggplot(event.df_invtime, aes(time)) + 
    geom_histogram(aes(y=..ncount..), colour="blue", fill="blue", breaks=seq(0, x_max,length.out=100))
    p_2 <- p2 + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    grid_layout <- rbind(c(1,3), c(2,NA))
    grid_width <- c(5,1)
    grid_heigth <- c(5,1)

    p <-grid.arrange(
      grobs=list(p1, p2,legend),
      layout_matrix = grid_layout,
      widths = grid_width,
      heights = grid_heigth)
    return(p)
}