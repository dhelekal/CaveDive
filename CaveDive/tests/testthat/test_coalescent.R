context("Homogenous Likelihood")

test_that("Computed likelihood matches ground truth", {
  expect_equal(homogenous_coal.log_lh(c(1,2,3,4,5),c(0.1, 1.1, 2.1, 3.1),1.2), -3.729286, tolerance=1e-6)
})

test_that("Computed likelihood matches simulation likelihood", {
  sam <- c(1,2,3,4,5,7,8,9)*100
  Neg <- 20
  co <-homogenous_coal.simulate(sam, Neg)
  log_lh_tree <- co$log_likelihood
  log_lh <- homogenous_coal.log_lh(sam,co$coalescent_times,Neg)
  expect_equal(log_lh_tree,log_lh)
})

context("Trees")
test_that("Coalescent Tree matches precomputed tree", {
  set.seed(1)
  ord <- trunc(runif(2*4,1,2+1)) ### sample twice per coalescent event
  ord <- ord[seq(1,8,by=2)] ### discard even samples
  
  nodes <- seq(1,5)
  nodes <- sapply(nodes, function (x) return (paste0("S", x)))
  
  sam <- c(0,2,4,6,7)
  co <-c(5,3,1,-1)
  
  set.seed(1)
  tr <- build_coal_tree(sam,co)
  gt_str <- nodes[1]
  
  for (i in c(2:length(sam))) {
    s<-c(1:2)
    s[ord[i-1]] <- paste0(gt_str,":2")
    s[3-ord[i-1]] <- paste0(nodes[i], ":1")
    gt_str<-paste0("(",s[1], ",",s[2],")")
  }
  
  gt_tr <-paste0(gt_str,";")
  expect_identical(tr, gt_tr)
})