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
            
            expect_equal(exp.lh(rate, s), inhomogenous_exp.lh(Neg_t, Neg_t.int, 0, s))
            expect_equal(poi_0.lh(rate, s), inhomogenous_poi_0.lh(Neg_t.int, 0, s))
            expect_equal(exp.prob(rate, s), inhomogenous_exp.prob(Neg_t.int, 0, s))
            
          })

context("Homogenous Likelihood")
test_that("Computed likelihood matches ground truth", {
  expect_equal(homogenous_coal.log_lh(c(1, 2, 3, 4, 5), c(0.1, 1.1, 2.1, 3.1), 1.2),-3.729286, tolerance =
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
test_that("Computed likelihood matches simulation likelihood with constant Neg", {
  sam <- seq(1,100,2) * 100
  Neg <- 2000
  
  Neg_t <- function (s)
    return (1 / Neg)
  Neg_t.int <- function (t, s)
    return(s / Neg)
  Neg_t.inv_int <- function(t, s)
    return(s * Neg)
  
  co <- inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
  
  log_lh <- co$log_likelihood
  times <- co$coalescent_times
  comp_lh <- inhomogenous_coal.log_lh(sam,
                                      times,
                                      Neg_t,
                                      Neg_t.int)
  
  expect_equal(comp_lh, log_lh, tolerance = 1e-6)
})

test_that("Computed likelihood matches simulation likelihood with exponential Neg", {
  
  sam <- seq(1,100,2) * 100
  lambda <- 1/400
  N <- 1e6
  
  Neg_t <- function (s)
    return (1/N * exp(lambda*s))
  Neg_t.int <- function (t, s) {
      out <- Inf
      if (!(exp(lambda*(t+s)) == Inf)) {
        out <- (1/(lambda*N))*(exp(lambda*(t+s)) - exp(lambda*t))
      }
    return (out)
  }
  Neg_t.inv_int <- function(t, s)
    return((1/lambda)*log(lambda*N*s*exp(-lambda*t) + 1))
  
  co <- inhomogenous_coal.simulate(sam, Neg_t, Neg_t.int, Neg_t.inv_int)
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
            sam <- seq(1,100,2) * 100
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
            expect_equal(comp_lh, comp_lh.gt)
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
  co <- c(5, 3, 1, -1)
  
  set.seed(1)
  tr <- build_coal_tree(sam, co)
  gt_str <- nodes[1]
  
  for (i in c(2:length(sam))) {
    s <- c(1:2)
    s[ord[i - 1]] <- paste0(gt_str, ":2")
    s[3 - ord[i - 1]] <- paste0(nodes[i], ":1")
    gt_str <- paste0("(", s[1], ",", s[2], ")")
  }
  
  gt_tr <- paste0(gt_str, ";")
  expect_identical(tr, gt_tr)
})