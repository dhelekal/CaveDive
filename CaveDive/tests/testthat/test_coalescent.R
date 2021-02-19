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

test_that("Inhomogenous process with exponential Neg matches rescaled Homogenous process",
          {
            sam <- rep(0, 10)
            lambda <- 1000
            N <- 1
            
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
            
            set.seed(1)
            co <-
              inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
            log_lh <- co$log_likelihood
            times <- co$coalescent_times
            
            set.seed(1)
            co.gt <- homogenous_coal.simulate(sam, 1)
            times.rescaled <-
              rescale_to_exponential(co.gt$coalescent_times, lambda, 0)
            log_lh.rescaled <- inhomogenous_coal.log_lh(sam,
                                                        times.rescaled,
                                                        Neg_t,
                                                        Neg_t.int)
            expect_equal(log_lh.rescaled, log_lh)
            expect_equal(times.rescaled, times)
            
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

test_that("Native exponential likelihood matches simulation likelihood from plot_exp_growth", 
  {
    set.seed(1)
    sam <- runif(100, 0, 10)
    sam <- sam - max(sam)
    sam <- sam[order(-sam)]

    N = runif(1, 1, 100)
    lambda = runif(1, 0.1, 10)

    co <- plot_exp_growth(sam, lambda, N)
    log_lh <-co$log_likelihood
    times <- co$coalescent_times
    times <- times[order(-times)]

    log_lh_native <- exponential_coalescent_loglh(sam, times, lambda, N)
    expect_equal(log_lh, log_lh_native)
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
test_that("structured_coal.preprocess_phylo works", {
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

    tr.nodiv <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)

    pre <- structured_coal.preprocess_phylo(tree.nodiv)

    expect_equal(all(sapply(pre$edges.df$id,
                function(i) pre$nodes.df$times[pre$edges.df$node.child[i]]>pre$nodes.df$times[pre$edges.df$node.parent[i]])),
                TRUE)
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

    tr.nodiv <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)
    
    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[101:length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:3), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs <- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- structured_coal.preprocess_phylo(tree.nodiv)
    subtrees <- lapply(c(1:length(MRCAs)), function (x) pre$clades.list[[(nodeid(pre$phy, MRCAs[x])-100)]]) 
    times <- extract_lineage_times(pre, MRCAs, div_times) 

    for (i in div_cols){
      sam.gt <- sam[which(colours==i)]
      coal.gt  <- co$times[which(co$colours==i)]

      for (j in c(1:length(co$div_from))) {
        if (co$div_from[j]==i) sam.gt<-c(sam.gt, div_times[j])
      }
      sam.gt <- sam.gt[order(-sam.gt)]

      expect_equal(sam.gt, times$sam.times[[i]])
      expect_equal(coal.gt, times$coal.times[[i]])
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

    tr.nodiv <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)

    pre <- structured_coal.preprocess_phylo(tree.nodiv)

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

    tr.nodiv <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)
    
    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[101:length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:3), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs <- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- structured_coal.preprocess_phylo(tree.nodiv)
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

test_that("Outbreak simulation likelihood matches computer expansion likelihood", {
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

    out <- expansions_simulate(poi_rate, concentration, sam, r_mean, r_sd, K_mean, K_sd, time_rate, time_shape)
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
    comp_log_lh <- expansions_likelihood(pre, MRCAs, out$div_times, out$A, out$K, out$N, out$exp_probs)
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
    expect_equal(comp_log_lh-lgamma(length(out$div_times)), sim_log_lh)
})

context("RJMCMC")
test_that("Transdimensional moves are balanced", {
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

    tr.nodiv <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)
    pre <- structured_coal.preprocess_phylo(tree.nodiv)

    r_mean <- 0 ## growth rate lognormal prior mean
    N_mean <- 4 ## carrying capacity rate lognormal prior mean

    r_sd <- 2 ## growth rate lognormal prior sd
    K_sd <- 1 ## carrying capacity  rate lognormal prior sd
    N_sd <- 3 ## parent population size lognormal prior sd

    kappa <- 1/10
    nu <- 1/2

    prior_i <- function(x) dpois(x, 1, log = TRUE) ### poisson 1 prior

    prior_N <- function(x) dlnorm(x, meanlog = N_mean, sdlog = N_sd, log = TRUE)
    prior_N.sample <- function() rlnorm(1, meanlog = N_mean, sdlog = N_sd) 

    prior_r <- function(x) dlnorm(x, meanlog = r_mean, sdlog = r_sd, log = TRUE) 
    prior_r.sample <- function() rlnorm(1, meanlog = r_mean, sdlog = r_sd) 

    prior_K_given_N <- function(x, N) dlnorm(x, meanlog = log(N), sdlog = K_sd, log = TRUE)
    prior_K_given_N.sample <- function(N) rlnorm(1, meanlog = log(N), sdlog = K_sd) 

    prior_t_given_N <- function(x, N) {
           if (all(x < 0)) {
                  out <- dgamma(-x, shape=(nu^2)/kappa, scale = kappa * N / nu, log = TRUE)
           } else {
                  out <- -Inf
           }
           return(out) 
    } 

    prior_t_given_N.sample <- function(N) (-rgamma(1, shape=(nu^2)/kappa, scale = kappa * N / nu))

    para.initialiser <- function(N, prior_r, prior_K_given_N, prior_t_given_N, pre){
      edges <- pre$edges.df
      nodes <- pre$nodes.df
      x_next <- vector(mode = "list", length = 4)

      rates <- prior_r()
      K <- prior_K_given_N(N)
      div.times <- prior_t_given_N(N)

      ### which branches exist at time of divergence
      br_extant_before <- edges$id[which(nodes$times[edges$node.child] > div.times)]
      br_extant_after <- br_extant_before[which(nodes$times[edges$node.parent[br_extant_before]] < div.times)]

      ### filter out terminal branches as those have 0 prior mass
      extant_inner <- br_extant_after[which(edges$node.child[br_extant_after]>pre$n_tips)] 
      ### choose one at random

      if(length(extant_inner) < 1) extant_inner <- c(NA)


      div.branch <- extant_inner[sample.int(length(extant_inner),1)] 

      x_next[[1]] <- rates
      x_next[[2]] <- K
      x_next[[3]] <- div.times
      x_next[[4]] <- div.branch
      return(x_next)
    }

    para.log_lh <- function(x, N, prior_r, prior_K_given_N, prior_t_given_N, pre) {
      edges <- pre$edges.df
      nodes <- pre$nodes.df

      rates <- x[[1]]
      K <- x[[2]]
      div.times <- x[[3]]
      div.branch <- x[[4]]

      ### which branches exist at time of divergence
      br_extant_before <- edges$id[which(nodes$times[edges$node.child] > div.times)]
      br_extant_after <- br_extant_before[which(nodes$times[edges$node.parent[br_extant_before]] < div.times)]

      ### filter out terminal branches as those have 0 prior mass
      extant_inner <- br_extant_after[which(edges$node.child[br_extant_after]>pre$n_tips)] 
      if (length(extant_inner) < 1) extant_inner <- c(NA)

      out <- prior_r(rates) + prior_K_given_N(K,N) + prior_t_given_N(div.times,N) + log(1/length(extant_inner))
      return(out)
    }

    fn_log_J <- function(i_prev, x_prev, x_next) {
      return(0)
    }

    fn_log_J_inv <- function(i_prev, x_prev, x_next, which_mdl_rm) {
      return(0)
    }

    p.init <- function(N) para.initialiser(N, prior_r.sample, prior_K_given_N.sample, prior_t_given_N.sample, pre)
    p.log_lh <- function(x, N) para.log_lh(x, N, prior_r, prior_K_given_N, prior_t_given_N, pre)

    i_0 <- 0 
    x_0 <- list()
    x_0[[1]] <- 100 
    x_0[[2]] <- c(1)

    x_prev <- x_0
    i_prev <- i_0

    ### Increase dim up to 10 always check QR both when adding a dim and then remove an event and check QR
    for (i in c(1:10)){
        prop_up <- transdimensional.sampler(x_prev, i_prev, pre, p.init, p.log_lh, fn_log_J, fn_log_J_inv, fixed_move=1)
        x_next <- prop_up$x
        i_next <- prop_up$i
    
        qr_comp <- prop_up$qr
        qr_test <- -prop_lh(x_prev, i_prev, x_next, i_next, pre, p.log_lh) + prop_lh(x_next, i_next, x_prev, i_prev, pre, p.log_lh)

        expect_equal(qr_comp, qr_test)

        x_prev <- x_next
        i_prev <- i_next

        prop_down <- transdimensional.sampler(x_prev, i_prev, pre, p.init, p.log_lh, fn_log_J, fn_log_J_inv, fixed_move=2)
        qr_comp2 <- prop_down$qr
        x_next2 <- prop_down$x
        i_next2 <- prop_down$i

        qr_test2 <- -prop_lh(x_prev, i_prev, x_next2, i_next2, pre, p.log_lh) + prop_lh(x_next2, i_next2, x_prev, i_prev, pre, p.log_lh)

        expect_equal(qr_comp2, qr_test2)
    }
    ### x_prev now 10+1-dimensional, test each move 5 times and check QR
    for (k in c(1:5)) {
      ### 5 moves
      for (l in c(1:10)) {
          prop_within <- within_model.sampler(x_prev, i_prev, pre, 10, fixed_move=k) 
          x_next <- prop_within$x
          i_next <- i_prev
          qr_comp <- prop_within$qr

          qr_test <- -prop_lh(x_prev, i_prev, x_next, i_next, pre, p.log_lh) + prop_lh(x_next, i_next, x_prev, i_prev, pre, p.log_lh)
          expect_equal(qr_comp, qr_test)
      }
    }
})