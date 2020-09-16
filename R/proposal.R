
prop.sampler <-function (x_prev, pre){
  which_move <- sample.int(4, size=1)
  if (which_move==1) {
    x_next <- move_1(x_prev, pre)
  } else if(which_move==2) {
    x_next <- move_2(x_prev, pre)
  } else if(which_move==3) {
    x_next <- move_3(x_prev, pre)
  } else if(which_move==4) {
    x_next <- move_4(x_prev, pre)
  }
  return(x_next)
}

prop.cond_log_lh <- function(x_cand, x_prev, pre) {

  rates_prev <- x_prev[[1]]
  K_prev <- x_prev[[2]]
  N_prev <- x_prev[[3]]
  div.times_prev <- x_prev[[4]]
  div.branch_prev <- x_prev[[5]]

  rates_cand <- x_cand[[1]]
  K_cand <- x_cand[[2]]
  N_cand <- x_cand[[3]]
  div.times_cand <- x_cand[[4]]
  div.branch_cand <- x_cand[[5]]

  outgoing <- pre$outgoing
  edges <- pre$edges.df

  param_log_lh <- sum(dnorm(rates_cand, mean=rates_prev, sd=1, log=TRUE)) + 
         sum(dnorm(K_cand, mean=K_prev, sd=1, log=TRUE)) +
         sum(dnorm(N_cand, mean=N_prev, sd=1, log=TRUE))

  time_log_lh <- 0
  branch_log_lh <- 0
  for (i in c(1:length(div.branch_cand))) {
    if (div.branch_cand[i] == div.branch_prev[i]) {
      time_log_lh <- time_log_lh + dnorm(div.times_cand[i],
                                          mean=div.times_prev[i],
                                          sd=1)#sd=max(edges$length[div.branch_cand[i]], 1))
    } else {
      time_log_lh <- time_log_lh - log(edges$length[div.branch_cand[i]])
      if (!all(outgoing[[edges$node.child[div.branch_prev[i]]]] != div.branch_cand[i])) {
        branch_log_lh <- branch_log_lh +
         log(edges$length[div.branch_cand[i]]/sum(edges$length[outgoing[[edges$node.child[div.branch_prev[i]]]]]))
      }
    }
  }
  return(param_log_lh+time_log_lh+branch_log_lh)
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
  
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  N <- x_prev[[3]]
  div.times <- x_prev[[4]]
  div.branch <- x_prev[[5]]

  which_idx <- sample.int(length(div.times), size=1)
  div.times_upd <- div.times
  div.times_upd <- rnorm(1, mean=div.times[which_idx], sd=1)#sd=max(pre$edges.df$length[div.branch[which_idx]]))

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- N
  x_next[[4]] <- div.times_upd
  x_next[[5]] <- div.branch

  return(x_next)
}

move_3 <- function(x_prev, pre) { ### update branch
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

  div.branch_upd <- div.branch
  div.times_upd <- div.times
  which_idx <- sample.int(length(div.branch), size=1)
  ## decide whether moving up or down
  r1 <- runif(1,1,3)

  if (r1 < 2) { ###up
    if (edges$node.parent[div.branch[which_idx]]!=root) {
        div.branch_upd[which_idx] <- incoming[[edges$node.parent[div.branch[which_idx]]]]
      } else {
        div.branch_upd[which_idx] <- outgoing[[root]][which(outgoing[[root]]!=div.branch[which_idx])]
      }
  } else { ###down
    r2 <- runif(1,0,sum(edges$length[outgoing[[edges$node.child[div.branch[which_idx]]]]]))

    if (r2 < edges$length[outgoing[[edges$node.child[div.branch[which_idx]]]]][1]) {
        div.branch_upd[which_idx] <- outgoing[[edges$node.child[div.branch[which_idx]]]][1]
    } else {
        div.branch_upd[which_idx] <- outgoing[[edges$node.child[div.branch[which_idx]]]][2]  
    }
  }

  div.times_upd[which_idx] <- runif(1, nodes$times[edges$node.parent[div.branch_upd[which_idx]]], nodes$times[edges$node.child[div.branch_upd[which_idx]]])

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- N
  x_next[[4]] <- div.times_upd
  x_next[[5]] <- div.branch_upd

  return(x_next)
}

move_4 <- function(x_prev, pre) { ### update N
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