context("Likelihood")

test_that("Computed likelihood matches ground truth", {
  expect_equal(tree_likelihood(c(1,2,3,4,5),c(0.1, 1.1, 2.1, 3.1),1.2), -3.729286, tolerance=1e-6)
})

test_that("Computed likelihood matches simulation likelihood", {
  sam <- c(1,2,3,4,5,7,8,9)*100
  Neg <- 20
  co <-simulate(sam, Neg)
  log_lh_tree <- co$log_likelihood
  log_lh <- tree_likelihood(sam,co$coalescent_times,Neg)
  expect_equal(log_lh_tree,log_lh)
})