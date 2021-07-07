context("Utils")
test_that("Homogenous and Inhomogenous exponentials agree for constant rate",
          {
            Neg <- 5
            rate <- 1 / Neg
            s <- 2
            
            Neg_t <- function (s)
              return (1 / Neg)
            Neg_t.int <- function (t, s)
              return(s / Neg)
            Neg_t.inv_int <- function(t, s)
              return(s * Neg)
            
            expect_equal(exp.loglh(rate, s),
                         inhomogenous_exp.loglh(Neg_t, Neg_t.int, 0, s))
            expect_equal(poi_0.loglh(rate, s),
                         inhomogenous_poi_0.loglh(Neg_t.int, 0, s))
            expect_equal(exp.prob(rate, s), inhomogenous_exp.prob(Neg_t.int, 0, s))
            
          })

context("Homogenous Likelihood")
test_that("Computed likelihood matches ground truth", {
  expect_equal(homogenous_coal.log_lh(c(1, 2, 3, 4, 5), c(0.1, 1.1, 2.1, 3.1), 1.2), -3.729286, tolerance =
                 1e-6)
})

test_that("Computed likelihood matches simulation likelihood", {
  sam <- c(1, 2, 3, 4, 5, 7, 8, 9) * 100
  Neg <- 2000
  co <- homogenous_coal.simulate(sam, Neg)
  log_lh_tree <- co$log_likelihood
  log_lh <- homogenous_coal.log_lh(sam, co$coalescent_times, Neg)

  expect_equal(log_lh_tree, log_lh)
})

context("Inhomogenous Likelihood")
test_that("Computed likelihood matches simulation likelihood with constant Neg",
          {
            sam <- seq(1, 100, 2) * 100
            Neg <- 2000
            
            Neg_t <- function (s)
              return (1 / Neg)
            Neg_t.int <- function (t, s)
              return(s / Neg)
            Neg_t.inv_int <- function(t, s)
              return(s * Neg)
            
            co <-
              inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
            
            log_lh <- co$log_likelihood
            times <- co$coalescent_times
            comp_lh <- inhomogenous_coal.log_lh(sam,
                                                times,
                                                Neg_t,
                                                Neg_t.int)
            
            expect_equal(comp_lh, log_lh, tolerance = 1e-6)
          })

test_that("Computed likelihood matches simulation likelihood with exponential Neg",
          {
            sam <- seq(1, 100, 2) * 100
            lambda <- 1 / 400
            N <- 1e6
            
            Neg_t <- function (s)
              return (1 / N * exp(lambda * s))
            Neg_t.int <- function (t, s) {
              out <- Inf
              if (!(exp(lambda * (t + s)) == Inf)) {
                out <- (1 / (lambda * N)) * (exp(lambda * (t + s)) - exp(lambda * t))
              }
              return (out)
            }
            Neg_t.inv_int <- function(t, s)
              return((1 / lambda) * log(lambda * N * s * exp(-lambda * t) + 1))
            
            co <-
              inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
            log_lh <- co$log_likelihood
            times <- co$coalescent_times
            comp_lh <- inhomogenous_coal.log_lh(sam,
                                                times,
                                                Neg_t,
                                                Neg_t.int)
            
            expect_equal(comp_lh, log_lh)
          })

test_that("Homogenous process matches inhomogenous process for constant Neg",
          {
            sam <- seq(1, 100, 2) * 100
            Neg <- 2000
            set.seed(1)
            co <- homogenous_coal.simulate(sam, Neg)
            gt.log_lh <- co$log_likelihood
            gt.times <- co$coalescent_times
            
            Neg_t <- function (s)
              return (1 / Neg)
            Neg_t.int <- function (t, s)
              return(s / Neg)
            Neg_t.inv_int <- function(t, s)
              return(s * Neg)
            set.seed(1)
            in.co <-
              inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
            
            log_lh <- in.co$log_likelihood
            times <- in.co$coalescent_times
            
            comp_lh <- inhomogenous_coal.log_lh(sam,
                                                times,
                                                Neg_t,
                                                Neg_t.int)
            
            comp_lh.gt <- homogenous_coal.log_lh(sam, gt.times, Neg)
            
            expect_equal(gt.log_lh, log_lh)
            expect_equal(gt.times, times)
            expect_equal(comp_lh.gt, comp_lh)
          })


test_that("Linear growth likelihood matches precomputed ground truth",
          {
            sam <- c(3,2)
            times <- c(1)

            Neg_t <- function (t) 1 / (3*2.2 - 2.2*t)
            Neg_t.int <- function (t,s) -1/2.2*(log(2.2*(3-t-s))-log(2.2*(3-t)))

            comp_lh <- inhomogenous_coal.log_lh(sam,
                                                times,
                                                Neg_t,
                                                Neg_t.int)

            expect_equal(comp_lh, -1.103524, tolerance = 1e-6)
})

context("Native Likelihood")
test_that("Native homogenous likelihood matches simulation likelihood", { 
  sam <- runif(100,0,10)
  sam <- sam - max(sam)
  sam <- sam[order(-sam)]
  Neg <- 100
  co <- homogenous_coal.simulate(sam, Neg)
  log_lh_tree <- co$log_likelihood

  log_lh.native <- coalescent_loglh(sam, co$coalescent_times, Neg, 0)

  expect_equal(log_lh_tree, log_lh.native)
})

test_that("Native exponential likelihood matches simulation likelihood",
          {
            sam <- seq(0,0.001, by=0.0001)
            lambda <- 1000
            N <- 10
            
            Neg_t <- function (s)
              return (1 / N * exp(lambda * s))
            Neg_t.int <- function (t, s) {
              out <- Inf
              if (!(exp(lambda * (t + s)) == Inf)) {
                out <- (1 / (lambda * N)) * (exp(lambda * (t + s)) - exp(lambda * t))
              }
              return (out)
            }
            Neg_t.inv_int <- function(t, s)
              return((1 / lambda) * log(lambda * N * s * exp(-lambda * t) + 1))
            
            co <-
              inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
            log_lh_tree <- co$log_likelihood
            times <- co$coalescent_times

            log_lh.r <- inhomogenous_coal.log_lh(sam,
                                                times,
                                                Neg_t,
                                                Neg_t.int)

            log_lh.native <- exponential_coalescent_loglh(sam[order(-sam)], times[order(-times)], lambda, N)

            expect_equal(log_lh_tree, log_lh.r)  
            expect_equal(log_lh_tree, log_lh.native)  
          })

context("Trees")
test_that("Coalescent Tree matches precomputed tree", {
  set.seed(1)
  ord <-
    trunc(runif(2 * 4, 1, 2 + 1)) ### sample twice per coalescent event
  ord <- ord[seq(1, 8, by = 2)] ### discard even samples
  
  nodes <- seq(1, 5)
  nodes <- sapply(nodes, function (x)
    return (paste0("S", x)))
  
  sam <- c(0, 2, 4, 6, 7)
  co <- c(5, 3, 1,-1)
  
  set.seed(1)
  tr <- build_coal_tree(sam, co)
  gt_str <- nodes[1]
  
  for (i in c(2:length(sam))) {
    s <- c(1:2)
    s[ord[i - 1]] <- paste0(gt_str, ":2")
    s[3 - ord[i - 1]] <- paste0(nodes[i], ":1")
    gt_str <- paste0("(", s[1], ",", s[2], ")","N", i-1)
  }
  
  gt_tr <- paste0(gt_str, ";")
  expect_identical(tr, gt_tr)
})

context("Structured Coalescent")
test_that("preprocess_phylo works", {
    set.seed(1)
    
    sam <- runif(100, 0, 10)
    tmax <- max(sam)
    sam <- sam - tmax
    sam <- sam[order(-sam)]
    colours <- trunc(runif(100, 1, 4))
    
    N <- 100
    A <- c(5.1, 0.3)
    K <- c(100,100)

    div_times <- c(-25, -60, -Inf)
    div_cols <- c(1, 2, 3)
    rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)

    tr <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)

    pre <- preprocess_phylo(tree.nodiv)

    expect_equal(all(sapply(pre$edges.df$id,
                function(i) pre$nodes.df$times[pre$edges.df$node.child[i]]>pre$nodes.df$times[pre$edges.df$node.parent[i]])),
                TRUE)
  })

  test_that("preprocess_phylo edge indexing is invariant of node indexing if order_edges_by_node_label is true.", {
    set.seed(1)
    
    sam <- runif(100, 0, 10)
    tmax <- max(sam)
    sam <- sam - tmax
    sam <- sam[order(-sam)]
    colours <- trunc(runif(100, 1, 4))
    
    N <- 100
    A <- c(5.1, 0.3)
    K <- c(100,100)

    div_times <- c(-25, -60, -Inf)
    div_cols <- c(1, 2, 3)
    rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)

    tr <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)

    pre <- preprocess_phylo(ladderize(tree.nodiv, right=T))
    pre_l <- preprocess_phylo(ladderize(tree.nodiv, right=F))

    expect_equal(pre$phy$node.label[pre$edges.df$node.child], pre_l$phy$node.label[pre_l$edges.df$node.child])
    expect_equal(pre$phy$tip.label[pre$edges.df$node.child], pre_l$phy$tip.label[pre_l$edges.df$node.child])
  })

test_that("extract_lineage_times works", {
    set.seed(1)
    
    sam <- runif(100, 0, 10)
    tmax <- max(sam)
    sam <- sam - tmax
    sam <- sam[order(-sam)]
    colours <- trunc(runif(100, 1, 4))
    
    N <- 100
    A <- c(5.1, 0.3)
    K <- c(100,100)

    div_times <- c(-25, -60, -Inf)
    div_cols <- c(1, 2, 3)
    rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)

    tr <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)
    
    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[101:length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:3), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs <- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- preprocess_phylo(tree.nodiv)
    subtrees <- lapply(c(1:length(MRCAs)), function (x) pre$clades.list[[(nodeid(pre$phy, MRCAs[x])-100)]]) 
    times <- extract_lineage_times(pre, MRCAs, div_times, T)
    times3 <- extract_lineage_times_native(pre, MRCAs, div_times)

    for (i in div_cols){
      sam.gt <- sam[which(colours==i)]
      coal.gt  <- co$times[which(co$colours==i)]

      for (j in c(1:length(co$div_from))) {
        if (co$div_from[j]==i) sam.gt<-c(sam.gt, div_times[j])
      }
      sam.gt <- sam.gt[order(-sam.gt)]

      expect_equal(sam.gt, times$sam.times[[i]])
      expect_equal(coal.gt, times$coal.times[[i]])
      expect_equal(times$sam.times[[i]], times3$sam.times[[i]])
      expect_equal(times$coal.times[[i]], times3$coal.times[[i]])
      expect_equal(times$partition_counts[[i]], times3$partition_counts[[i]])
    }
})

test_that("Structured Coalescent likelihood with only one clade matches neutral likelihood", {
    set.seed(1)

    n_tips <- 100
    
    sam <- runif(n_tips, 0, 10)
    tmax <- max(sam)
    sam <- sam - tmax
    sam.ord <- order(-sam)
    sam <- sam[sam.ord]
    colours <- trunc(runif(100, 1, 4))
    
    N <- 100
    A <- c(5.1, 0.3)
    K <- c(100,100)

    div_times <- c(-25, -60, -Inf)
    div_cols <- c(1, 2, 3)
    rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)

    tr <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)

    pre <- preprocess_phylo(tree.nodiv)

    root_MRCA <- pre$phy$node.label[pre$nodes.df$id[which.min(pre$nodes.df$times)]-n_tips]
    root_div <- -Inf

    comp <- structured_coal.likelihood(pre, root_MRCA, root_div, list(), list(), 100, type="Sat")
    hom_log_lh <- coalescent_loglh(sam, co$times, 100, 0)

    expect_equal(comp$log_lh, hom_log_lh)
})

test_that("Simulation Likelihood matches product of colour specific likelihoods", {
    set.seed(1)
    
    sam <- runif(100, 0, 10)
    tmax <- max(sam)
    sam <- sam - tmax
    sam.ord <- order(-sam)
    sam <- sam[sam.ord]
    colours <- trunc(runif(100, 1, 4))
    
    N <- 100
    A <- c(5.1, 0.3)
    K <- c(100,100)

    div_times <- c(-25, -60, -Inf)
    div_cols <- c(1, 2, 3)
    rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)

    tr <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)
    
    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[101:length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:3), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs <- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- preprocess_phylo(tree.nodiv)
    comp <- structured_coal.likelihood(pre, MRCAs, div_times, A, K, N, type="Sat")
    log_lh <- 0
    for (i in div_cols){
      sam.gt <- sam[which(colours==i)]
      coal.gt  <- co$times[which(co$colours==i)]

      for (j in c(1:length(co$div_from))) {
        if (co$div_from[j]==i) {
          sam.gt <- c(sam.gt, div_times[j])
        }
      }


      sam.gt <- sam.gt[order(-sam.gt)]
      coal.gt <- coal.gt[order(-coal.gt)]

      expect_equal(sam.gt, comp$sam.times[[i]])
      expect_equal(coal.gt, comp$coal.times[[i]])
      if (i < length(div_cols)){
        log_lh <- log_lh + sat_coalescent_loglh(sam.gt, coal.gt, div_times[i], A[i], K[i], 0)
        expect_equal(sat_coalescent_loglh(sam.gt, coal.gt, div_times[i], A[i], K[i], 0),
                      sat_coalescent_loglh(comp$sam.times[[i]], comp$coal.times[[i]], div_times[i], A[i], K[i], 0))
      } else {
        log_lh <- log_lh + coalescent_loglh(sam.gt, coal.gt, N, 0)
        expect_equal(coalescent_loglh(sam.gt, coal.gt, N, 0),
                      coalescent_loglh(comp$sam.times[[i]], comp$coal.times[[i]], N, 0))
      }
    }

    log_lh <- log_lh - lgamma(length(div_times))

    expect_equal(log_lh, co$log_lh)
    expect_equal(comp$log_lh-lgamma(length(div_times)), co$log_lh)
})

test_that("Expansion simulation likelihood matches computed expansion likelihood", {
    set.seed(1123456)
    
    n_tips <- 100
    concentration <- 2

    sam <- runif(n_tips, 0, 0.1)
    sam <- sam - max(sam)
    sam <- sam[order(-sam)]

    priors <- standard_priors()

    out <- expansions_simulate(priors, sam, concentration=concentration, given=list(N=100))
    co <- out$co
    params <- out$params

    tr <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)
    
    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[(n_tips+1):length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:(params$n_exp+1)), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs <- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- preprocess_phylo(tree.nodiv)
    sim_log_lh <- out$coal_log_lh

    A <- sapply(params$t_mid, function(x) (1/x)**2)

    comp <- structured_coal.likelihood(pre, MRCAs, params$div_times, A, params$K, params$N)

    expect_equal(comp$log_lh, sim_log_lh)

    for (i in params$div_cols){
      sam.gt <- sam[which(params$tip_colours==i)]
      coal.gt  <- co$times[which(co$colours==i)]

      for (j in c(1:length(co$div_from))) {
        if (co$div_from[j]==i) {
          sam.gt <- c(sam.gt, params$div_times[j])
        }
      }

      sam.gt <- sam.gt[order(-sam.gt)]
      coal.gt <- coal.gt[order(-coal.gt)]

    expect_equal(sam.gt, comp$sam.times[[i]])
    expect_equal(coal.gt, comp$coal.times[[i]])
  }
})

context("RJMCMC")
test_that("Transdimensional moves are balanced", {
    set.seed(1)
    
    sam <- runif(100, 0, 10)
    tmax <- max(sam)
    sam <- sam - tmax
    sam.ord <- order(-sam)
    sam <- sam[sam.ord]

    s_priors <- standard_priors()

    sim <- expansions_simulate(s_priors, sam, 1, given=list(N=100))

    co <- sim$co
    params <- sim$params

    tr <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)
    pre <- preprocess_phylo(tree.nodiv)

    edges <- pre$edges.df
    nodes <- pre$nodes.df

    fn_log_J <- function(i_prev, x_prev, x_next) {
      return(0)
    }

    fn_log_J_inv <- function(i_prev, x_prev, x_next, which_mdl_rm) {
      return(0)
    }

    max_t <- max(nodes$times)
    min_t <- min(nodes$times)

    max_t_nodes <- max(nodes$times[which(!nodes$is_tip)])

    tree_height <- max_t - min_t

    div_time_prop.sample <- function(N) runif(1, min_t, max_t_nodes)
    div_time_prop.lh <- function(x, N) log(1/(max_t_nodes-min_t))

    p.init <- function(N) para.initialiser(N, s_priors$prior_t_mid_given_N.sample, s_priors$prior_K_given_N.sample, div_time_prop.sample, pre)
    p.log_lh <- function(x, N) para.log_lh(x, N, s_priors$prior_t_mid_given_N, s_priors$prior_K_given_N, div_time_prop.lh, pre)

    i_0 <- 0 
    x_0 <- list()
    x_0[[1]] <- 100 
    x_0[[2]] <- c(1)

    x_prev <- x_0
    i_prev <- i_0

    ### Increase dim up to 10 always check QR both when adding a dim and then remove an event and check QR
    for (i in c(1:20)){
        prop_up <- transdimensional.sampler(x_prev, i_prev, pre, p.init, p.log_lh, fn_log_J, fn_log_J_inv, fixed_move=1)
        x_next <- prop_up$x_next
        i_next <- prop_up$i_next
    
        qr_comp <- prop_up$qr
        qr_test <- -prop_lh(x_prev, i_prev, x_next, i_next, pre, p.log_lh) + prop_lh(x_next, i_next, x_prev, i_prev, pre, p.log_lh)

        x_prev <- prop_up$x_next
        i_prev <- prop_up$i_next

        expect_equal(abs(qr_comp) < Inf, TRUE)
        expect_equal(qr_comp, qr_test)

        prop_down <- transdimensional.sampler(x_prev, i_prev, pre, p.init, p.log_lh, fn_log_J, fn_log_J_inv, fixed_move=2)
        qr_comp2 <- prop_down$qr
        x_next2 <- prop_down$x_next
        i_next2 <- prop_down$i_next

        qr_test2 <- -prop_lh(x_prev, i_prev, x_next2, i_next2, pre, p.log_lh) + prop_lh(x_next2, i_next2, x_prev, i_prev, pre, p.log_lh)
        
        expect_equal(abs(qr_comp2) < Inf, TRUE)
        expect_equal(qr_comp2, qr_test2)
    }
    ### x_prev now 10+1-dimensional, test each move 5 times and check QR
    for (k in c(1,2,3,4)) {
      ### 5 moves
      for (l in c(1:50)) {
          attempts <- 0
          x_next <- NA
          i_next <- NA
          qr_comp <- -Inf
          while(!(qr_comp > -Inf) && attempts < 10) {
            prop_within <- within_model.sampler(x_prev, i_prev, pre, tree_height/2, fixed_move=k) 
            x_next <- prop_within$x_next
            i_next <- i_prev
            qr_comp <- prop_within$qr
            attempts <- attempts + 1
          }

          if (attempts==10) {
            print(x_next)
          }

          expect_equal((attempts < 10), TRUE)
          qr_test <- -prop_lh(x_prev, i_prev, x_next, i_next, pre, p.log_lh, tree_height/2) + prop_lh(x_next, i_next, x_prev, i_prev, pre, p.log_lh, tree_height/2)
          expect_equal(qr_comp, qr_test)
      }
    }
})

test_that("Log-posterior returns correct values", {
    set.seed(1)

    n_tips <- 100
 
    sam <- runif(n_tips, 0, 0.1)
    sam <- sam - max(sam)
    sam <- sam[order(-sam)]

    priors <- standard_priors() 
    
    out <- expansions_simulate(priors, sam, concentration=2, given=list(N=100, n_exp=2))
    params <- out$params
    co <- out$co 

    prior_probs <- function(probs) {

        out <- -Inf
        if (abs(sum(probs)-1) < 1e-8 && all(probs > 0)) {
            out <- ddirichlet(t(matrix(probs)), alpha=rep(2, length(probs)), log=TRUE)
        }
        return(out) 
    }
    
    tr <- build_coal_tree.structured(sam, co$times, params$tip_colours, co$colours, params$div_times, params$div_cols, co$div_from)
    tree <- read.tree(text = tr$full)
    tree.nodiv <- collapse.singles(tree)
    pre <- preprocess_phylo(tree.nodiv)

    root_set <- rep(NA, params$n_exp)
    if(params$n_exp > 0) {
        for (i in c(1:params$n_exp)) {
            L <- LETTERS[i]
            N_set <- nodeid(pre$phy, pre$phy$node.label[grep(paste0("N_",L), pre$phy$node.label)]) 
            set_root <- N_set[which.min(pre$nodes.df$times[N_set])]
            set_root.edge <- pre$incoming[[set_root]]
            root_set[i] <- set_root.edge
        }
    }

    x <- list()

    i <- params$n_exp
    x[[1]] <- params$N
    x[[2]] <- params$exp_probs

    for(j in c(1:params$n_exp)) {
      expansion_para <- list()
      expansion_para[[1]] <- params$t_mid[j]
      expansion_para[[2]] <- params$K[j]
      expansion_para[[3]] <- params$div_times[j]
      expansion_para[[4]] <- root_set[j]
      x[[2+j]] <- expansion_para
    }

    log_p <- log_posterior(x,
      i, 
      priors$prior_i, 
      priors$prior_N, 
      priors$prior_t_mid_given_N, 
      priors$prior_K_given_N, 
      priors$prior_t_given_N, 
      prior_probs, 
      pre)

    expect_equal(log_p$prior, out$param_log_lh)
    expect_equal(log_p$lh, out$coal_log_lh)
})