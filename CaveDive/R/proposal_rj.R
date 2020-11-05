###list of: N, div_event_1, ... div_event_N
###div_event list of r, K, x_0, branch
offset <- 1

prop.sampler <- function (x_prev, i_prev, pre, para.initialiser){
  p <- 1 
  r <- runif(1,0,2) 

  if (r < p) { ## transdimensional
    upd <- transdimensional.sampler(x_prev, i_prev, pre, para.initialiser)
    x_next <- upd$x_next
    i_next <- upd$i_next
  } else { ## within-model
    x_next <- within_model.sampler(x_prev, i_prev, pre)
    i_next <- i_prev
  }
  return(list(x_next=x_next,i_next=i_next))
}

prop.cond_lh <- function(x, i, x_given, i_given, initialiser.log_lh, pre) {
  N <- x[[1]]
  N_given <- x_given[[1]]
  out <- dnorm(N, mean=N_given, sd=1, log=TRUE)

  if(i > i_given) {
    out <- out + initialiser.log_lh(x[[i+offset]])
  } else if(i==i_given && i > 0) {
    out <- out + sum(sapply(c(1:i_given), function(j) within_model.cond_log_lh(x[[j+offset]], x_given[[j+offset]], pre)))
  }
  return(out)
}

transdimensional.sampler <- function(x_prev, i_prev, pre, para.initialiser) {
  which_move <- sample.int(2,size=1)
  if (which_move==1) { ### increase dim
    i_next <- i_prev + 1
    x_next <- x_prev
    x_next[[length(x_next) + offset]] <- para.initialiser()
   } else { ### decrease dim
    x_next <- x_prev
    if (i_prev > 0) {
      which_elem <- sample.int(i_prev, size=1)
      x_next <- x_next[-(which_elem+offset)]
      i_next <- i_prev - 1
    } else {
      i_next <- i_prev
    }
  }
  return(list(x_next=x_next, i_next=i_next))
}

within_model.sampler <- function(x_prev, i_prev, pre) {

  if (i_prev > 0) {
    which_move <- sample.int(4, size=1)
  } else {
    which_move <- 4
  }

  if (which_move==1) {
    x_next <- make_move(x_prev, i_prev, pre, move_update_rates)
  } else if(which_move==2) {
    x_next <- make_move(x_prev, i_prev, pre, move_update_time)
  } else if(which_move==3) {
    x_next <- make_move(x_prev, i_prev, pre, move_update_branch)
  } else if(which_move==4) {
    N_prev <- x_prev[[1]]
    N_next <- move_update_N(N_prev, pre)

    x_next <- x_prev
    x_next[[1]] <- N_next 
  }
  return(x_next)
}

make_move <- function(x_prev, i_prev, pre, move) {
    which_model <- offset+sample.int(i_prev, size=1)
    mdl_prev <- x_prev[[which_model]] 
    mdl_next <- move(mdl_prev, pre)
    x_next <- x_prev
    x_next[[which_model]] <- mdl_next
    return(x_next)
}

within_model.cond_log_lh <- function(x, given, pre) {
  rates_given <- given[[1]]
  K_given <- given[[2]]
  div.times_given <- given[[3]]
  div.branch_given <- given[[4]]

  rates <- x[[1]]
  K <- x[[2]]
  div.times <- x[[3]]
  div.branch <- x[[4]]

  outgoing <- pre$outgoing
  edges <- pre$edges.df

  param_log_lh <- sum(dnorm(rates, mean=rates_given, sd=1, log=TRUE)) + 
         sum(dnorm(K, mean=K_given, sd=1, log=TRUE)) 

  time_log_lh <- 0
  branch_log_lh <- 0
  for (i in c(1:length(div.branch))) {
    if (div.branch[i] == div.branch_given[i]) {
      time_log_lh <- time_log_lh + dnorm(div.times[i],
                                          mean=div.times_given[i],
                                          sd=1)#sd=max(edges$length[div.branch[i]], 1))
    } else {
      branch_log_lh <- log(1/2)
      if (!all(outgoing[[edges$node.child[div.branch_given[i]]]] != div.branch[i])) {
        branch_log_lh <- branch_log_lh + log(1/2)
      }
    }
  }
  return(param_log_lh+time_log_lh+branch_log_lh)
}

move_update_rates <- function(x_prev, pre) { ### update rates
  
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

  which_idx <- sample.int(length(rates), size=1)

  K_upd <- K
  rates_upd <- rates

  rates_upd[which_idx] <- rnorm(1, mean=rates[which_idx], sd=1) 
  K_upd[which_idx]  <- rnorm(1, mean=K[which_idx], sd=1)

  x_next[[1]] <- rates_upd 
  x_next[[2]] <- K_upd
  x_next[[3]] <- div.times
  x_next[[4]] <- div.branch

  return(x_next)
}

move_update_time <- function(x_prev, pre) { ### update time
  
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

  which_idx <- sample.int(length(div.times), size=1)
  div.times_upd <- div.times
  div.times_upd <- rnorm(1, mean=div.times[which_idx], sd=1)#sd=max(pre$edges.df$length[div.branch[which_idx]]))

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- div.times_upd
  x_next[[4]] <- div.branch

  return(x_next)
}

move_update_branch <- function(x_prev, pre) { ### update branch
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

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
  r1 <- runif(1,0,2)

  if (r1 < 1) { ###up
    if (edges$node.parent[div.branch[which_idx]]!=root) {
        div.branch_upd[which_idx] <- incoming[[edges$node.parent[div.branch[which_idx]]]]
      } else {
        div.branch_upd[which_idx] <- outgoing[[root]][which(outgoing[[root]]!=div.branch[which_idx])]
      }
  } else { ###down
    r2 <- runif(1,0,2)

    if (r2 < 1) {
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
  x_next[[3]] <- div.times_upd
  x_next[[4]] <- div.branch_upd

  return(x_next)
}


move_update_N <- function(N_prev, pre) {

  N_upd <- rnorm(1, mean = N_prev, sd = 1)
  return(N_upd)

} 