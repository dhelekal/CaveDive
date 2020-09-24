
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
  edges <- pre$edges.df

  param_log_lh <- sum(dnorm(rates, mean=rates_given, sd=1, log=TRUE)) + 
         sum(dnorm(K, mean=K_given, sd=1, log=TRUE)) +
         sum(dnorm(N, mean=N_given, sd=1, log=TRUE))

  time_log_lh <- 0
  branch_log_lh <- 0
  for (i in c(1:length(div.branch))) {
    if (div.branch[i] == div.branch_given[i]) {
      time_log_lh <- time_log_lh + dnorm(div.times[i],
                                          mean=div.times_given[i],
                                          sd=1)#sd=max(edges$length[div.branch[i]], 1))
    } else {
      branch_log_lh <- log(1/2)
      #time_log_lh <- time_log_lh - log(edges$length[div.branch[i]])
      if (!all(outgoing[[edges$node.child[div.branch_given[i]]]] != div.branch[i])) {
        branch_log_lh <- branch_log_lh + log(1/2)
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

  K_upd <- K
  rates_upd <- rates

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
  which_half <- pre$which_half

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
    r2 <- runif(1,1,3)

    if (r2 < 2) {
        div.branch_upd[which_idx] <- outgoing[[edges$node.child[div.branch[which_idx]]]][1]
    } else {
        div.branch_upd[which_idx] <- outgoing[[edges$node.child[div.branch[which_idx]]]][2]  
    }
  }

  old_len <- edges$length[div.branch[which_idx]]
  new_len <- edges$length[div.branch_upd[which_idx]]

  rel_pos <- (div.times[which_idx]-nodes$times[edges$node.parent[div.branch[which_idx]]])/old_len
  div.times_upd[which_idx] <- nodes$times[edges$node.parent[div.branch_upd[which_idx]]] + new_len*rel_pos
  
  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- N
  x_next[[4]] <- div.times_upd
  x_next[[5]] <- div.branch_upd

  return(x_next)
}


move_4 <- function(x_prev, pre) {
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  N <- x_prev[[3]]
  div.times <- x_prev[[4]]
  div.branch <- x_prev[[5]]

  N_upd <- rnorm(1, mean = N, sd = 1)

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- N_upd
  x_next[[4]] <- div.times
  x_next[[5]] <- div.branch

  return(x_next)
}