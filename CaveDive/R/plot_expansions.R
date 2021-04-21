plot_persistence <- function(mcmc.df, event.df, pre, prior_t_given_N=NULL, correlates=NULL) {
     p_mat <- compute_persistence(pre, event.df)

     p_mat[lower.tri(p_mat)]<-NA
     p_df <- melt(p_mat)
     names(p_df) <- c("sample_1", "sample_2", "value")
     p_df$sample_1 <- factor(x = p_df$sample_1,
                                    levels = pre$phy$tip.label,#[ord], 
                                    ordered = TRUE)
     p_df$sample_2 <- factor(x = p_df$sample_2,
                                    levels = pre$phy$tip.label,#[ord], 
                                    ordered = TRUE)

     heat_map <- ggplot(data = p_df, aes(x = sample_1, y = sample_2)) +
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
             legend.position = c(0.8,0.1),
             legend.title = element_text(angle = -90),
             aspect.ratio=1)
     corr_map <- NULL
     if(!is.null(correlates)) {
          stopifnot("Number correlate rows must match number of tips"=nrow(correlates)==pre$n_tips)
          stopifnot("Correlates must be a data.frame with rownames equal to tip labels"=rownames(correlates)[order(rownames(correlates))]==pre$phy$tip.label[order(pre$phy$tip.label)])
          r_df <- correlates
          r_df$tip_id <- rownames(r_df) 
          r_df <- melt(r_df)
          r_df$tip_id <- factor(x = r_df$tip_id,
                          levels = pre$phy$tip.label, 
                          ordered = TRUE)

          corr_map <- ggplot(data = r_df, aes(x = tip_id, y = variable)) +
                         geom_tile(aes(fill = value)) +
                         scale_fill_viridis(option= "viridis", na.value="gray50" , discrete=!is.double(r_df$value)) +
                         theme_minimal() +
                         guides(fill=guide_legend(title.position = "right"))+
                         coord_flip()+
                         theme(axis.title.y = element_blank(), 
                               axis.text.y = element_blank(), 
                               axis.ticks.y = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               text = element_text(size=30),
                               legend.position="right",
                               legend.title = element_text(angle = -90))
     }
     tree_map <- plot_tree(pre, event.df)
     tree_map_t <- tree_map + coord_flip() + scale_x_reverse()
     blank <- ggplot() + theme_minimal()
     ggarrange(blank, tree_map_t, blank, tree_map, heat_map, corr_map, widths=c(1,3,0.75), heights=c(1,3))
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
   tree.full <- full_join(tree, ldf, by = 'node')

   x_max <- -min(pre$nodes.df$times)

   p1 <- ggtree(tree.full, aes(color=frequency), size=0.75, ladderize=F) +
   geom_point() +
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
   return(p1)
}

plot_summary <- function (model_data, expansion_data, phylo_preprocessed, priors, modes=NULL) {
     hist_dim <- ggplot(model_data, aes(dim)) +  
        geom_histogram(aes(y = stat(count / sum(count))), binwidth=1) + 
        theme_bw() +
        scale_fill_brewer(palette="Dark2")  + 
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              text = element_text(size=20))
     hist_N <- ggplot(model_data, aes(N)) +
         geom_histogram(aes(y = stat(count / sum(count))), bins=100) +
         theme_bw() + 
         scale_fill_brewer(palette="Dark2")  + 
         theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
               text = element_text(size=20))


     hist_br <- ggplot(expansion_data, aes(x=factor(br)))
        if(is.null(modes)) {
          hist_br <- hist_br + geom_bar(aes(y = stat(count / sum(count)))) 
        } else {
          hist_br <- hist_br + geom_bar(aes(y = stat(count / sum(count)), fill=is.mode)) + 
          scale_fill_brewer(palette="Dark2") + 
                     geom_text(stat="count", aes(label = mode_clade, y= ((..count..)/sum(..count..))), vjust = -.25, hjust=-0.1, size=11, color="red")
        }

        hist_br <- hist_br + theme_bw() + 
        labs(x="Branch Number",fill="Expansion Root") +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = c(0.8, 0.2),
              text = element_text(size=20))

     if (is.null(modes)) MRCA_lab=NULL else MRCA_lab=phylo_preprocessed$edges.df$node.child[modes]
     tree_freq <- plot_tree_freq(model_data, 
                              expansion_data, 
                              phylo_preprocessed, 
                              prior_t_given_N=function (x, n) exp(priors$prior_t_given_N(x,n)), 
                              highlight_node=NULL, 
                              MRCA_lab=MRCA_lab)
     
     grid_layout <- rbind(c(1, 2), c(3,4))
     grid_width <- c(2,2)
     grid_heigth <- c(2,2)
     summary_panel <- grid.arrange(
        grobs=list(hist_N, hist_dim, hist_br, tree_freq),
        layout_matrix = grid_layout,
        widths = grid_width,
        heights = grid_heigth)

     return(summary_panel)
}

plot_traces <- function(model_data, expansion_data) {
     max_it <- max(model_data$it)
     min_it <- min(model_data$it)
     trace_lh <- ggplot(model_data, aes(x=it, y=lh)) +
          geom_line(alpha = 0.3) +
          theme_bw() + 
          ylab("log-likelihood") +
          xlim(c(min_it,max_it)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
     trace_prior <- ggplot(model_data, aes(x=it, y=prior)) +
          geom_line(alpha = 0.3) +
          theme_bw() + 
          ylab("log-prior") +
          xlim(c(min_it,max_it)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
     trace_N <- ggplot(model_data, aes(x=it, y=N)) +
          geom_line(alpha = 0.3) +
          theme_bw() + 
          ylab("N") +
          xlim(c(min_it,max_it)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
     trace_dim <- ggplot(model_data, aes(x=it, y=dim)) +
          geom_line(alpha = 0.3) +
          theme_bw() +
          ylab("Number of Expansions") +
          xlim(c(min_it,max_it)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
     trace_br <- ggplot(expansion_data, aes(x=it, y=br)) + 
          geom_point(alpha=0.1,size=0.1) + 
          theme_bw() + 
          ylab("Branch") +
          xlim(c(min_it,max_it)) +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())

     g <- gridExtra::gtable_rbind(ggplotGrob(trace_lh), 
                           ggplotGrob(trace_prior), 
                           ggplotGrob(trace_N),
                           ggplotGrob(trace_dim),
                           ggplotGrob(trace_br))

     panels <- g$layout$t[grep("panel", g$layout$name)]
     g$heights[panels] <- unit(c(1,1,1,1,3), "null")
     
     grid.newpage()
     grid.draw(g)
} 

prior_mixture <- function(prior, cond_values) {
     f_mixture <- function (X) sapply(X, function (x) (1/length(cond_values))*sum(sapply(cond_values, function(y) prior(x, y))))
     return(stat_function(fun=f_mixture, colour="purple", size=2))
}

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

#' Plots tree branch frequency histogram 
#' 
#' @param mcmc.df mcmc.df returned by mcmc2data.frame
#' @param event.df event.df returned by mcmc2data.frame
#' @param pre Preprocessed phylogeny
#' @param prior_t_given_N (Optional) Expansion time prior. If supplied prior will be overlayed.
#' @param highlight_node (Optional) Node to highlight
#' @return a plot object
#' @export
plot_tree_freq <- function(mcmc.df, event.df, pre, prior_t_given_N=NULL, highlight_node=NULL, MRCA_lab=c()) {
   tree <- pre$phy
   freq <- table(event.df$br)

   labs <- c(tree$node.label, tree$tip.label)
   tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
   ids <- nodeid(tree, labs)
   id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) 0.0 else if (is.na(freq[paste0(pre$incoming[[i]])])) 0.0 else freq[paste0(pre$incoming[[i]])])

   ldf <- data.frame(node = ids, frequency = id_freq, tip=tip)
   ldf$edge_id <- sapply(ldf$node, function(i) pre$incoming[[i]])
   ldf$lab <- sapply(ldf$node, function (x) {
                    a <- MRCA_lab[which(MRCA_lab==x)]
                    if(length(a) > 0) return(a[1]) else return(NA)
               })
   tree.full <- full_join(tree, ldf, by = 'node')

   x_max <- -min(pre$nodes.df$times)

   p1 <- ggtree(tree.full, aes(color=frequency), size=0.75, ladderize=F) +
   geom_point() +
   geom_text2(aes(label=edge_id, 
                 subset=!is.na(lab), 
                 x=branch), color="red", size=12, vjust=-1) +
   scale_size_manual(values=c(1)) +
   scale_x_continuous(limits=c(0, x_max)) +
   scale_color_viridis(option="plasma") +
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

plot_mode_summary <- function(mcmc.df, event.df, priors) {
  mode_br_df <- event.df[which(event.df$is.mode),]
  mode_br_mcmc_df <- mode_dim_marginal[mode_dim_marginal$it %in% mode_br_df$it, ]

  dummy_gt <- data.frame(br=unique(mode_br_df$br))
  dummy_gt$median.K <- sapply(dummy_gt$br, function (x) median(mode_br_df$K[which(mode_br_df$br==x)]))
  dummy_gt$median.t_mid <- sapply(dummy_gt$br, function (x) median(mode_br_df$t_mid[which(mode_br_df$br==x)]))

  K_ci <- lapply(dummy_gt$br, function (x) compute_ci(mode_br_df$K[which(mode_br_df$br==x)]))
  t_mid_ci <- lapply(dummy_gt$br, function (x) compute_ci(mode_br_df$t_mid[which(mode_br_df$br==x)]))

  dummy_gt$ci_lo.K <- sapply(K_ci, function (x) x[1])
  dummy_gt$ci_hi.K <- sapply(K_ci, function (x) x[2])

  dummy_gt$ci_lo.t_mid <- sapply(t_mid_ci, function (x) x[1])
  dummy_gt$ci_hi.t_mid <- sapply(t_mid_ci, function (x) x[2])

  br.labs <- sapply(unique(mode_br_df$br), function (x) paste0("Branch: ",x))
  names(br.labs) <- unique(mode_br_df$br)

  K_facet <- ggplot(mode_br_df) + 
             geom_histogram(aes(x=K, y = stat(count / sum(count))), bins=100) +
             prior_mixture(function(x,N) exp(priors$prior_K_given_N(x,N)),mode_br_mcmc_df$N) +
             geom_rect(data = dummy_gt, aes(xmin = ci_lo.K, xmax = ci_hi.K), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.3) +
             facet_wrap(~br, labeller=labeller(br = br.labs)) +
             labs(x="Carrying Capacity") +
             theme_bw() +
             theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),text = element_text(size=20))

  t_mid_facet <- ggplot(mode_br_df) + 
             geom_histogram(aes(x=t_mid, y = stat(count / sum(count))), bins=100) +
             prior_mixture(function(x,N) exp(priors$prior_t_mid_given_N(x,N)),mode_br_mcmc_df$N) +
             geom_rect(data = dummy_gt, aes(xmin = ci_lo.t_mid, xmax = ci_hi.t_mid), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.3) +
             facet_wrap(~br, labeller=labeller(br = br.labs)) +
             labs(x="Time to Midpoint") +
             theme_bw() +
             theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),text = element_text(size=20))

  grid.arrange(
          grobs=list(K_facet, t_mid_facet),
          nrow=2,
          heights = c(1,1))
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

   grid_layout <- rbind(c(1, 2), c(3,3), c(4,4))
   grid_width <- c(2,2)
   grid_heigth <- c(2,1,1)

   correct_dim <- mcmc.df[which(mcmc.df$dim > 0),]
   correct_dim_it <- correct_dim$it
   event_dim_marginal <- event.df[unlist(sapply(correct_dim_it, function (x) which(event.df$it==x))),]
   br_subs <- which(event_dim_marginal$br == which_br)
   event_br_marginal <- event_dim_marginal[br_subs,]
   event_br_marginal$it <- correct_dim_it[br_subs]
   event_br_marginal$idx <- if(length(event_br_marginal$it) > 0) c(1:length(event_br_marginal$it)) else c()

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

   mrca <- pre$edges.df$node.child[which_br]
   tree_freq <- plot_tree_freq(mcmc.df, event.df, pre, prior_t_given_N=prior_t_given_N, highlight_node=mrca)

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
   return(list(event_panel=event_panel, tree_highlight_panel=tree_panel))
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