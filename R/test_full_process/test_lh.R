library(CaveDive)
library(ape)
library(testthat)
set.seed(1123456)
    
    n_tips <- 100
    poi_rate <- 2
    concentration <- 2

    sam <- runif(n_tips, 0, 0.1)
    sam <- sam - max(sam)
    sam <- sam[order(-sam)]

    r_mean <- 0
    r_sd <- 1

    K_mean <- 6
    K_sd <- 0.2

    time_shape <- floor(n_tips/3)
    time_rate <- 5**(-1)

    out <- outbreaks_simulate(poi_rate, concentration, sam, r_mean, r_sd, K_mean, K_sd, time_rate, time_shape)
    co <- out$co

    tr.nodiv <- build_coal_tree.structured(sam, co$times, out$tip_colours, co$colours, out$div_times, out$div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)
    
    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[(n_tips+1):length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:out$n_exp), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs <- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- structured_coal.preprocess_phylo(tree.nodiv)
    comp_log_lh <- outbreaks_likelihood(pre, MRCAs, out$div_times, out$A, out$K, out$N, out$exp_probs)
    sim_log_lh <- out$full_lh

    comp <- structured_coal.likelihood(pre, MRCAs, out$div_times, out$A, out$K, out$N)

    for (i in out$div_cols){
      sam.gt <- sam[which(out$tip_colours==i)]
      coal.gt  <- co$times[which(co$colours==i)]

      for (j in c(1:length(co$div_from))) {
        if (co$div_from[j]==i) {
          sam.gt <- c(sam.gt, out$div_times[j])
        }
      }

      sam.gt <- sam.gt[order(-sam.gt)]
      coal.gt <- coal.gt[order(-coal.gt)]

    expect_equal(sam.gt, comp$sam.times[[i]])
    expect_equal(coal.gt, comp$coal.times[[i]])
  }
    expect_equal(comp_log_lh, sim_log_lh)