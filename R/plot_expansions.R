plot_persistence <- function(mcmc.df, event.df, pre, axis_titles=list(), legend_titles=list(), correlates=list(), modes=NULL) {
     phy <- pre$phy
     if (is.null(modes)) MRCA_lab=NULL else MRCA_lab=pre$edges.df$node.child[modes]
     tree_map <- plot_tree(pre, event.df, MRCA_lab)

     dat <- tree_map[["data"]]
     dat <- dat[dat$isTip,]
     tip.ord <- order(dat$y)

     p_mat <- compute_persistence(pre, event.df)
     p_mat <- p_mat[dat$label[tip.ord],dat$label[tip.ord]]
     p_mat[lower.tri(p_mat)]<-NA
     p_df <- melt(p_mat)
     names(p_df) <- c("sample_1", "sample_2", "value")
     p_df$sample_1 <- factor(p_df$sample_1,
                              levels = dat$label[tip.ord], 
                              ordered = TRUE)
     p_df$sample_2 <- factor(p_df$sample_2,
                              levels = dat$label[tip.ord], 
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
             legend.position = c(0.8,0.2),
             legend.title = element_text(angle = -90, hjust=0.5),
             aspect.ratio=1)

     blank <- ggplot() + theme_minimal()
     corr_maps <- list(blank)
     n_cor <- length(correlates)
     if(n_cor > 0) {
          corr_maps <- sapply(c(1:n_cor), function(i) list(build_correlate_map(correlates[[i]], pre, dat, tip.ord, unlist(axis_titles[i]), unlist(legend_titles[i]))))
     }
     tree_map_t <- tree_map + coord_flip() + scale_x_reverse()

     do.call("ggarrange",c(c(list(blank), list(tree_map_t), rep(list(blank), max(1, n_cor)), list(tree_map), list(heat_map), corr_maps), 
                              list(widths=c(1,3,rep(0.75, max(1,n_cor))), heights=c(1,3))))
}

build_correlate_map <- function(correlate, pre, dat, tip.ord, axis_title, leg_title) {
     stopifnot("Number correlate rows must match number of tips"=nrow(correlate)==pre$n_tips)
     stopifnot("Correlates must be a data.frame with rownames equal to tip labels"=rownames(correlate)[order(rownames(correlate))]==pre$phy$tip.label[order(pre$phy$tip.label)])
     r_df <- correlate
     r_df$tip_id <- rownames(r_df) 
     r_df <- melt(r_df, id.vars="tip_id")
     r_df$tip_id <- factor(x = r_df$tip_id,
                     levels = dat$label[tip.ord], 
                     ordered = TRUE)

     col_scheme <- NULL
     categ=!is.numeric(r_df$value)
     if (categ) {
          if (length(unique(r_df$value))==2){
               rdbu <- brewer.pal(n = 3, name = "RdBu")
               col_scheme <- scale_fill_manual(values=c(rdbu[1], rdbu[3]), na.value=rdbu[2])
          } else if (length(unique(r_df$value))<=12){
               col_scheme <- scale_fill_brewer(palette = "Paired", na.value="white")
          } else {
               col_scheme <- scale_fill_viridis(option= "viridis", na.value="white" , discrete=T)
          }
     } else {
          col_scheme <- scale_fill_viridis(option= "viridis", na.value="white" , discrete=!is.double(r_df$value))
     }
     axis_title_str <- " "
     leg_title_str <- " "

     if(!is.null(axis_title)) axis_title_str <- axis_title
     if(!is.null(leg_title)) leg_title_str <- leg_title

     corr_map <- ggplot(data = r_df, aes(x = tip_id, y = variable)) +
                    geom_tile(aes(fill = value)) +
                    col_scheme +
                    theme_minimal() +
                    guides(fill=guide_legend(title.position = "left"))+
                    coord_flip() +
                    labs(y=axis_title_str, fill=leg_title_str) +
                    theme(axis.title.y = element_blank(), 
                          axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          text = element_text(size=30),
                          axis.text.x = element_text(size=18, angle=45, hjust=1),
                          legend.position="bottom", legend.direction="vertical",
                          legend.title = element_text(angle = -90))
     return(corr_map)
}

plot_tree<-function(pre,event.df, MRCA_lab=NULL){
   tree <- pre$phy
   freq <- table(event.df$br) 

   labs <- c(tree$node.label, tree$tip.label)
   tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
   ids <- nodeid(tree, labs)
   id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) 0.0 else if (is.na(freq[paste0(pre$incoming[[i]])])) 0.0 else freq[paste0(pre$incoming[[i]])])

   ldf <- data.frame(node = ids, frequency = id_freq, tip=tip)
   ldf <- ldf[order(ldf$node),]
   ldf$edge_id <- sapply(ldf$node, function(i) pre$incoming[[i]])
   ldf$lab <- sapply(ldf$node, function (x) {
                    a <- MRCA_lab[which(MRCA_lab==x)]
                    if(length(a) > 0) return(a[1]) else return(NA)
               })
   tree.full <- full_join(tree, ldf, by = 'node')

   x_max <- -min(pre$nodes.df$times)

   p1 <- ggtree(tree.full, aes(color=frequency), size=0.75, ladderize=T) +
   geom_point() +
   geom_text2(aes(label=edge_id, 
                 subset=!is.na(lab), 
                 x=branch), color="red", size=12, vjust=-1, hjust=1) +
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
     
     hist_dim <- ggplot(model_data, aes(x=dim)) +
                 geom_bar(aes(y = ..prop..), stat="count") + 
                 geom_text(aes( label = scales::percent(..prop..), y= ..prop.. ), stat= "count", vjust = -.5, size=12) +
                 theme_bw() + 
                 xlab("Number of Expansions") + 
                 scale_y_continuous(labels=percent, limits=c(0,1)) +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                                  axis.title.y = element_blank(),
                                  text = element_text(size=20))
     hist_N <- ggplot(model_data, aes(N)) +
         geom_histogram(aes(y = stat(count / sum(count))), bins=100) +
         theme_bw() + 
         xlab("N") +
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
              text = element_text(size=14),
              axis.text.x = element_text(size=12, angle=45))

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
    p_mat <- matrix(data = 0.0, nrow = pre$n_tips, ncol = pre$n_tips, byrow = FALSE,
       dimnames = list(pre$phy$tip.label, pre$phy$tip.label))
    for(i in unique(df$it)) {
        subs_it <- df[which(df$it == i), ]
        partitions <- extract_lineage_times(pre, pre$phy$node.label[(c(pre$edges.df$node.child[subs_it$br],pre$root_idx)-pre$n_tips)], c(subs_it$time, -Inf), return_partitions=TRUE)$partitions
        for(p in partitions) {
          p_mat[p,p] <- p_mat[p,p] + 1
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
   ldf <- ldf[order(ldf$node),]
   ldf$edge_id <- sapply(ldf$node, function(i) pre$incoming[[i]])
   ldf$lab <- sapply(ldf$node, function (x) {
                    a <- MRCA_lab[which(MRCA_lab==x)]
                    if(length(a) > 0) return(a[1]) else return(NA)
               })
   tree.full <- full_join(tree, ldf, by = 'node')

   x_max <- -min(pre$nodes.df$times)

   p1 <- ggtree(tree.full, aes(color=frequency), size=0.75, ladderize=T) +
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

plot_mode_summary <- function(mcmc.df, event.df, priors, k_modes, gt.K=NULL, gt.t_mid=NULL) {
  mode_br_df <- event.df[which(event.df$is.mode),]
  mode_br_mcmc_df <- mcmc.df[mcmc.df$it %in% mode_br_df$it, ]

  dummy_gt <- data.frame(br=unique(mode_br_df$br))
  dummy_gt$median.K <- sapply(dummy_gt$br, function (x) median(mode_br_df$K[which(mode_br_df$br==x)]))
  dummy_gt$median.t_mid <- sapply(dummy_gt$br, function (x) median(mode_br_df$t_mid[which(mode_br_df$br==x)]))

  K_ci <- lapply(dummy_gt$br, function (x) compute_ci(mode_br_df$K[which(mode_br_df$br==x)]))
  t_mid_ci <- lapply(dummy_gt$br, function (x) compute_ci(mode_br_df$t_mid[which(mode_br_df$br==x)]))

  dummy_gt$ci_lo.K <- sapply(K_ci, function (x) x[1])
  dummy_gt$ci_hi.K <- sapply(K_ci, function (x) x[2])

  dummy_gt$ci_lo.t_mid <- sapply(t_mid_ci, function (x) x[1])
  dummy_gt$ci_hi.t_mid <- sapply(t_mid_ci, function (x) x[2])

  if(!is.null(gt.K)) {
     dummy_gt$gt.K <- gt.K
  }

  if(!is.null(gt.t_mid)) {
     dummy_gt$gt.t_mid <- gt.t_mid
  }

  br.labs <- sapply(unique(mode_br_df$br), function (x) paste0("Branch: ",x))
  names(br.labs) <- unique(mode_br_df$br)

  K_facet <- ggplot(mode_br_df) + 
             geom_histogram(aes(x=K, y = ..density..), bins=50) +
             prior_mixture(function(x,N) exp(priors$prior_K_given_N(x,N)),mode_br_mcmc_df$N) +
             geom_rect(data = dummy_gt, aes(xmin = ci_lo.K, xmax = ci_hi.K), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.3)+
             geom_vline(data=dummy_gt, aes(xintercept=median.K), colour="orange", linetype = "longdash", lwd=2)

  if(!is.null(gt.K)) K_facet <- K_facet + geom_vline(data=dummy_gt, aes(xintercept=gt.K), color="red", lwd=2)
 
  K_facet <- K_facet + 
             facet_wrap(~br, labeller=labeller(br = br.labs), scales="free") +
             labs(x="Carrying Capacity") +
             theme_bw() +
             theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),text = element_text(size=20))

  t_mid_facet <- ggplot(mode_br_df) + 
             geom_histogram(aes(x=t_mid, y = ..density..), bins=50) +
             prior_mixture(function(x,N) exp(priors$prior_t_mid_given_N(x,N)),mode_br_mcmc_df$N) +
             geom_rect(data = dummy_gt, aes(xmin = ci_lo.t_mid, xmax = ci_hi.t_mid), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.3) +
             geom_vline(data=dummy_gt, aes(xintercept=median.t_mid), colour="orange", linetype = "longdash",lwd=2)

  if(!is.null(gt.t_mid)) t_mid_facet <- t_mid_facet + geom_vline(data=dummy_gt, aes(xintercept=gt.t_mid),color="red", lwd=2)

  t_mid_facet <- t_mid_facet +
             facet_wrap(~br, labeller=labeller(br = br.labs), scales="free") +
             labs(x="Time to Midpoint") +
             theme_bw() +
             theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),text = element_text(size=20))

  grid.arrange(
          grobs=list(K_facet, t_mid_facet),
          nrow=2,
          heights = c(1,1))
}

plot_mode_traces <- function(mcmc.df, event.df, k_modes) {
  mode_br_df <- event.df[which(event.df$is.mode),]
  mode_br_mcmc_df <- mcmc.df[mcmc.df$it %in% mode_br_df$it, ]

  br.labs <- sapply(unique(mode_br_df$br), function (x) paste0("Branch: ",x))
  names(br.labs) <- unique(mode_br_df$br)
  
  max_it <- max(mode_br_mcmc_df$it)
  min_it <- min(mode_br_mcmc_df$it)

  K_facet <- ggplot(mode_br_df, aes(x=it, y=K)) +
               geom_line(alpha = 0.3) +
               theme_bw() + 
               xlim(c(min_it,max_it)) +
               theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     text = element_text(size=20)) 
  K_facet <- K_facet + 
             facet_wrap(~br, labeller=labeller(br = br.labs), scales="free") +
             labs(y="Carrying Capacity")


  t_mid_facet <- ggplot(mode_br_df, aes(x=it, y=t_mid)) +
               geom_line(alpha = 0.3) +
               theme_bw() + 
               xlim(c(min_it,max_it)) +
               theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     text = element_text(size=20)) 
  
  t_mid_facet <- t_mid_facet +
             facet_wrap(~br, labeller=labeller(br = br.labs), scales="free") +
             labs(y="Time to Midpoint")
  grid.arrange(
          grobs=list(K_facet, t_mid_facet),
          nrow=2,
          heights = c(1,1))
}

plot_pop_fn <- function(mcmc.df, event.df, which_br, t_max=NULL, eval_pts=100) {
  mode_br_df <- event.df[which(event.df$br==which_br),]
  mode_br_mcmc_df <- mcmc.df[mcmc.df$it %in% mode_br_df$it, ]
  min_x <- min(mode_br_df$time)
  if(is.null(t_max)) {
     max_x <- 0.3*abs(min_x)
  } else {
     max_x <- t_max
  }
  X <- seq(from=min_x, to=max_x, length.out=eval_pts)
  Y_med <- rep(0, eval_pts)
  Y_min <- rep(0, eval_pts)
  Y_max <- rep(0, eval_pts)

  funcs <- lapply(c(1:nrow(mode_br_df)), 
     function (i) function (s) 1/sat.rate(s, mode_br_df$K[i], (1/mode_br_df$t_mid[i])**2, mode_br_df$time[i]))

  for (i in c(1:eval_pts)) {
     f_vals <- sapply(c(1:length(funcs)), function(j) funcs[[j]](-X[i]))
     Y_med[i] <- median(f_vals)
     ci <- compute_ci(f_vals)
     Y_min[i] <- ci[1]
     Y_max[i] <- ci[2]
  }

  pal <- brewer.pal(n = 3, name = "Dark2")
  df <- data.frame(t=X, y_med=Y_med, y_min=Y_min, y_max=Y_max)
  gg <- ggplot(df) +
  geom_ribbon(aes(x=t,ymin=y_min, ymax=y_max),fill="grey50", alpha=0.3) +
  geom_line(data=subset(df, t <= 0), aes(x=t, y=y_med), linetype="solid",lwd=2, color=pal[2]) +
  geom_line(data=subset(df, t > 0), aes(x=t, y=y_med), linetype="longdash",lwd=2, color=pal[3]) + 
  theme_bw() +
  xlab("Time") +
  ylab("Neg") +
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  plot(gg)
}

plot_pop_fn_facet <- function(mcmc.df, event.df, k_modes, eval_pts=100, t_max=NULL, gt.K=NULL, gt.t_mid=NULL, gt.time=NULL) {
  mode_br_df <- event.df[which(event.df$is.mode),]
  mode_br_mcmc_df <- mcmc.df[mcmc.df$it %in% mode_br_df$it, ]

  brs <- unique(mode_br_df$br)
  min_x <- sapply(brs, function(x) min(mode_br_df$time[mode_br_df$br==x]))
  
  if(is.null(t_max)) {
     max_x <- 0.3*abs(min_x)
  } else {
     max_x <- t_max
  }

  Xseq <- lapply(c(1:k_modes), function(i) seq(from=min_x[i], to=max_x[i], length.out=eval_pts))

  Y_med <- c()
  Y_min <- c()
  Y_max <- c()
  X <- c()
  br_v <- c()

  for (k in c(1:length(brs))) {
     Y_med_tmp <- rep(0,eval_pts)
     Y_min_tmp <- rep(0,eval_pts)
     Y_max_tmp <- rep(0,eval_pts)

     br_subs <- mode_br_df[which(mode_br_df$br==brs[k]),]

     funcs <- lapply(c(1:nrow(br_subs)), 
     function (i) function (s) 1/sat.rate(s, br_subs$K[i], (1/br_subs$t_mid[i])**2, br_subs$time[i]))

     for (i in c(1:eval_pts)) {
          f_vals <- sapply(c(1:length(funcs)), function(j) funcs[[j]](-Xseq[[k]][i]))
          Y_med_tmp[i] <- median(f_vals)
          ci <- compute_ci(f_vals)
          Y_min_tmp[i] <- ci[1]
          Y_max_tmp[i] <- ci[2]
     }
     Y_med <- c(Y_med, Y_med_tmp)
     Y_min <- c(Y_min, Y_min_tmp)
     Y_max <- c(Y_max, Y_max_tmp)
     X <- c(X, Xseq[[k]])
     br_v <- c(br_v, rep(brs[k], eval_pts))
  }

  Y_val <- c()
  
  if(!is.null(gt.K) && !is.null(gt.t_mid) && !is.null(gt.time)) {
     for (k in c(1:length(brs))) {
        func <- function (s) 1/sat.rate(s, gt.K[k], (1/gt.t_mid[k])**2, gt.time[k])
        Y_val_tmp <- sapply(-Xseq[[k]], func) 

        Y_val <- c(Y_val, Y_val_tmp)
        X <- c(X, Xseq[[k]])
        br_v <- c(br_v, rep(brs[k], eval_pts))
     }
     dummy_df <- data.frame(t=X, y=Y_val, br=br_v)
  }


  br.labs <- sapply(brs, function (x) paste0("Branch: ",x))
  names(br.labs) <- unique(brs)

  pal <- brewer.pal(n = 3, name = "Dark2")
  df <- data.frame(t=X, y_med=Y_med, y_min=Y_min, y_max=Y_max, br=br_v)
  gg <- ggplot(df) +
  geom_ribbon(aes(x=t,ymin=y_min, ymax=y_max),fill="grey50", alpha=0.3) +
  geom_line(data=subset(df, t <= 0), aes(x=t, y=y_med), linetype="solid",lwd=2, color=pal[2]) +
  geom_line(data=subset(df, t > 0), aes(x=t, y=y_med), linetype="longdash",lwd=2, color=pal[3]) + 
  theme_bw() +
  xlab("Time") +
  ylab("Neg") +
  ylim(c(0, max(Y_max)))
     
  if(!is.null(gt.K) && !is.null(gt.t_mid) && !is.null(gt.time)) gg <- gg + 
       geom_line(data=dummy_df, aes(x=t, y=y), linetype="dotted",lwd=3, color=pal[1])

  gg <- gg + facet_wrap(~br, labeller=labeller(br = br.labs), scales="free") +
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  plot(gg)
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

compute_ci <- function(x, conf=0.95) {
  ci <- c()
  x_ord <- order(x)
  if(length(x)%%2==0) {
    l<-length(x)/2
    p1 <- x[x_ord][(l+1):length(x)]
    p2 <- x[x_ord][1:l]
  } else {
    l <-floor(length(x)/2)
    p1 <- x[x_ord][(l+2):length(x)]
    p2 <- x[x_ord][1:l]
  }
  ci[1] <- p2[l*(1-conf)]
  ci[2] <- p1[l*conf]
  return(ci)
}
