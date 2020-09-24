
prop.sampler <-function (x_prev, pre){
  which_move <- sample.int(3, size=1)
  if (which_move==1) {
    x_next <- move_1(x_prev, pre)
  } else if(which_move==2) {
    x_next <- move_2(x_prev, pre)
  }
  else if(which_move==3) {
    x_next <- move_3(x_prev, pre)
  } 
  return(x_next)
}

prop.cond_log_lh <- function(x, given, pre) {

  rates_given <- given[[1]]
  K_given <- given[[2]]
  N_given <- given[[3]]
  div.times_given <- given[[4]]
  div.branch_given <- given[[5]]

  rates <- x[[1]]
  K <- x[[2]]
  N <- x[[3]]
  div.times <- x[[4]]
  div.branch <- x[[5]]

  outgoing <- pre$outgoing
  incoming <- pre$incoming
  edges <- pre$edges.df
  nodes <- pre$nodes.df
  root <- pre$root
  which_half <- pre$which_half

  param_log_lh <- sum(dnorm(rates, mean=rates_given, sd=1, log=TRUE)) + 
         sum(dnorm(K, mean=K_given, sd=1, log=TRUE)) +
         sum(dnorm(N, mean=N_given, sd=1, log=TRUE))

  time_log_lh <- 0
  branch_log_lh <- 0

  for (i in c(1:length(div.branch))) {
    delta_t <- 0

    br <- div.branch[i]
    br_given <- div.branch_given[i]

    t <- div.times[i]
    t_given <- div.times_given[i]
    
    if (!(br == br_given) || !(t == t_given)) {
      root_traversed <- 1 ### 1 for FALSE
      if (which_half[br] != which_half[br_given]) { ### root traversed, switch to other half and update delta_t/t_given/br_given
        delta_t <- delta_t + (nodes$times[root] - t_given)
        br_given <- outgoing[[root]][which(which_half[outgoing[[root]]] != which_half[br_given])]
        t_given <- nodes$times[root]
        root_traversed <- -1 ### -1 for TRUE
      }
      delta_t <- delta_t + (t-t_given)*root_traversed
      if ((t-t_given) > 0){
        while (br!=br_given) {
          br <- incoming[[edges$node.parent[br]]]
          branch_log_lh <- branch_log_lh + log(1/2)
        }
      }
    }
    time_log_lh <- time_log_lh + dnorm(delta_t, mean=0, sd=1, log=TRUE)
  } 
  log_lh <- param_log_lh+time_log_lh+branch_log_lh 
  return(log_lh)
}

move_1 <- function(x_prev, pre) { ### update rates
  
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  N <- x_prev[[3]]
  div.times <- x_prev[[4]]
  div.branch <- x_prev[[5]]

  which_idx <- sample.int(length(rates), size=1)

  rates_upd <- rates
  K_upd <- K

  rates_upd[which_idx] <- rnorm(1, mean=rates[which_idx], sd=1) 
  K_upd[which_idx]  <- rnorm(1, mean=K[which_idx], sd=1)

  x_next[[1]] <- rates_upd
  x_next[[2]] <- K_upd
  x_next[[3]] <- N
  x_next[[4]] <- div.times
  x_next[[5]] <- div.branch

  return(x_next)
}

move_2 <- function(x_prev, pre) { ### update time
  
  log_lh <- 0 

  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  N <- x_prev[[3]]
  div.times <- x_prev[[4]]
  div.branch <- x_prev[[5]]

  edges <- pre$edges.df
  nodes <- pre$nodes.df
  outgoing <- pre$outgoing
  incoming <- pre$incoming
  root <- pre$root_idx

  which_idx <- sample.int(length(div.times), size=1)
  div.times_upd <- div.times
  div.branch_upd <- div.branch
  
  delta_t <- rnorm(1, mean = 0, sd = 1)
  t0 <- div.times[which_idx]
  br <- div.branch[which_idx]

  old_td <- delta_t

  root_traversed <- FALSE

  if (delta_t < 0) {
    while (nodes$times[edges$node.parent[br]] > (t0+delta_t) && edges$node.parent[br] != root) {
      br <- incoming[[edges$node.parent[br]]]
    }
    if (nodes$times[edges$node.parent[br]] > (t0+delta_t)) { ###reached root, but havent moved delta_t distance yet
      root_traversed <- TRUE
      if(edges$node.parent[br] != root) warning("Branch does not match time")
      delta_t <- nodes$times[root] - (t0+delta_t) ###update delta_t
      if (delta_t < 0) warning("Invalid delta_t")
      t0 <- nodes$times[root] ###change t0 to root
      br <- outgoing[[root]][which(outgoing[[root]]!=br)] ### switch to other root branch
    }
  }

  if (delta_t > 0) {
    while (!is.na(br) && nodes$times[edges$node.child[br]] < (t0+delta_t)) {
      r <- sample.int(2,size=1)
      br <- outgoing[[edges$node.child[br]]][r]
    }
  }

  new_time <- t0+delta_t
  new_branch <- br

  div.times_upd[which_idx] <- new_time
  div.branch_upd[which_idx] <- new_branch

  MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)
  MRCAs <- c(MRCAs, root_MRCA)
  div.times <- c(div.times, root_div)

  subtrees <- lapply(c(1:length(div.MRCA.nodes)), function (x) phylo.preprocessed$clades.list[[MRCA.idx[x]]]) 
  k_div <- length(subtrees)
  log_lh <- 0
  times <- extract_lineage_times(phylo.preprocessed, subtrees, div.MRCA.nodes, div.times)
  N_upd  <- rnorm(1, mean=N, sd=1)

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- N_upd
  x_next[[4]] <- div.times_upd
  x_next[[5]] <- div.branch_upd

  return(x_next)
}

move_3 <- function(x_prev, pre) { ### update N
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  N <- x_prev[[3]]
  div.times <- x_prev[[4]]
  div.branch <- x_prev[[5]]

  N_upd  <- rnorm(1, mean=N, sd=1)

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- N_upd
  x_next[[4]] <- div.times
  x_next[[5]] <- div.branch

  return(x_next)
}