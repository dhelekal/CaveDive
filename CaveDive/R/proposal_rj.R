###list of: N, div_event_1, ... div_event_N
###div_event list of r, K, x_0, branch
offset <- 2

prop.sampler <- function (x_prev, i_prev, pre, para.initialiser, initialiser.log_lh, fn_log_J, fn_log_J_inv){
  p <- 1 
  r <- runif(1,0,2) 

  if (r < p) { ## transdimensional
    x_prop <- transdimensional.sampler(x_prev, i_prev, pre, para.initialiser, initialiser.log_lh, fn_log_J, fn_log_J_inv)
    x_next <- x_prop$x_next
    i_next <- x_prop$i_next
    
    log_J <- x_prop$log_J
    qr <- x_prop$qr 
  } else { ## within-model
    x_prop <- within_model.sampler(x_prev, i_prev, pre)
    x_next <- x_prop$x_next
    i_next <- i_prev
    
    log_J <- 0
    qr <- x_prop$qr 
  }
  return(list(x_next=x_next,i_next=i_next, qr=qr, log_J=log_J))
}

transdimensional.sampler <- function(x_prev, i_prev, pre, para.initialiser, initialiser.log_lh, fn_log_J, fn_log_J_inv) {
  which_move <- sample.int(2, size=1)
  log_J <- 0
  qr <- 0
  if (which_move==1) { ### increase dim
    i_next <- i_prev + 1
    x_next <- x_prev
    
    x_next[[i_next + offset]] <- para.initialiser()

    old_probs <- x_prev[[2]]
    which_split <- sample.int((i_prev+1),1)
    u <- runif(1, 0, old_probs[which_split])
    new_probs <- c(old_probs, u)
    new_probs[which_split] <- new_probs[which_split] - u

    x_next[[2]] <- new_probs

    qr <- qr - log(1/(i_prev+1)) - log(1/old_probs[which_split])
    qr <- qr - initialiser.log_lh(x_next[[i_next + offset]]) 
    qr <- qr + log(1/(i_prev+1)) + log(1/(i_prev+1))

    if(abs(sum(x_next[[2]]) - 1) > 1e-8) warning("prob sum error")

    log_J <- fn_log_J(i_prev, x_prev, x_next)
   } else { ### decrease dim
    x_next <- x_prev
    if (i_prev > 0) {
      which_elem <- sample.int(i_prev, size=1)
      x_next <- x_next[-(which_elem+offset)]

      which_prob <- x_prev[[2]][which_elem+1]
      x_next[[2]] <- x_next[[2]][-(which_elem+1)]

      which_merge <- sample.int(i_prev, size=1)

      x_next[[2]][which_merge] <- x_next[[2]][which_merge] + which_prob
      
      i_next <- i_prev - 1

      qr <- qr - log(1/i_prev) ## proposal remove model 
      qr <- qr - log(1/i_prev) ## proposal merge with this probaility
      qr <- qr + initialiser.log_lh(x_prev[[(which_elem+offset)]]) ## reverse lh of adding that model
      qr <- qr + log(1/i_prev) + log(1/x_next[[2]][which_merge]) ##pick the same split from uniform
      
      if(abs(sum(x_next[[2]]) - 1) > 1e-8) warning("prob sum error")

      log_J <- fn_log_J_inv(i_prev, x_prev, x_next, which_elem)
    } else {
      i_next <- i_prev
    }
  }
  return(list(x_next=x_next, i_next=i_next, log_J=log_J, qr=qr))
}

within_model.sampler <- function(x_prev, i_prev, pre) {

  if (i_prev > 0) {
    which_move <- sample(c(1,2,3,4,5), size=1)
  } else {
    which_move <- 4
  }

  if (which_move==1) {
    x_prop <- make_move(x_prev, i_prev, pre, move_update_rates)
    x_next <- x_prop$x_next
    qr <- x_prop$qr
  } else if(which_move==2) {
    x_prop <- make_move(x_prev, i_prev, pre, move_update_time)
    x_next <- x_prop$x_next
    qr <- x_prop$qr
  } else if(which_move==3) {
    x_prop <- make_move(x_prev, i_prev, pre, move_update_branch)
    x_next <- x_prop$x_next
    qr <- x_prop$qr
  } else if(which_move==4) {
    N_prev <- x_prev[[1]]
    N_prop <- move_update_N(N_prev, pre)
    N_next <- N_prop$N_next
    qr <- N_prop$qr

    x_next <- x_prev
    x_next[[1]] <- N_next 
  } else if(which_move==5) {

    probs_prev <- x_prev[[2]]
    probs_prop <- move_update_probs(probs_prev, pre)
    probs_next <- probs_prop$probs_next
    qr <- probs_prop$qr

    x_next <- x_prev
    x_next[[2]] <- probs_next 
  }
  return(list(x_next=x_next,qr=qr))
}

make_move <- function(x_prev, i_prev, pre, move) {
    which_model <- offset+sample.int(i_prev, size=1)
    mdl_prev <- x_prev[[which_model]] 
    mdl_prop <- move(mdl_prev, pre)
    mdl_next <- mdl_prop$x_next
    qr <- mdl_prop$qr
    x_next <- x_prev
    x_next[[which_model]] <- mdl_next
    return(list(x_next=x_next, qr=qr))
}

move_update_rates <- function(x_prev, pre) { ### update rates
  
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

  if(length(rates) > 1) print(length(rates))

  K_upd <- K
  rates_upd <- rates

  rates_upd <- rnorm(1, mean=rates, sd=1) 
  K_upd  <- rnorm(1, mean=K, sd=1)

  qr <- -dnorm(rates_upd, mean=rates, sd=1, log=TRUE) - dnorm(K_upd, mean=K, sd=1, log=TRUE) ## proposal lh
  qr <- qr + dnorm(rates, mean=rates_upd, sd=1, log=TRUE) + dnorm(K, mean=K_upd, sd=1, log=TRUE) ## reverse lh

  x_next[[1]] <- rates_upd 
  x_next[[2]] <- K_upd
  x_next[[3]] <- div.times
  x_next[[4]] <- div.branch

  return(list(x_next=x_next,qr=qr))
}

move_update_time <- function(x_prev, pre) { ### update time
  
  x_next <- vector(mode = "list", length = length(x_prev))

  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

  div.times_upd <- rnorm(1, mean=div.times, sd=1)

  qr <- -dnorm(div.times_upd, mean=div.times, sd=1, log=TRUE) ##proposal lh
  qr <- qr + dnorm(div.times, mean=div.times_upd, sd=1, log=TRUE) ##reversal lh

  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- div.times_upd
  x_next[[4]] <- div.branch

  return(list(x_next=x_next, qr=qr))
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

  div.branch_upd <- div.branch
  div.times_upd <- div.times

  ## decide whether moving up or down
  r1 <- runif(1,0,2)
  qr <- -log(1/2) ## proposal direction lh
  qr <- qr + log(1/2) ## reverse lh
  if (r1 < 1) { ###up
    if (edges$node.parent[div.branch]!=root) {
        div.branch_upd <- incoming[[edges$node.parent[div.branch]]]
        qr <- qr + log(1/2) ## reverse lh
      } else {
        div.branch_upd <- outgoing[[root]][which(outgoing[[root]]!=div.branch)]
        ##no reverse  lh as direction effectively changes
      }
  } else { ###down
    r2 <- runif(1,0,2)
    qr <- qr - log(1/2) ## proposal which branch lh
    ## no reverse lh as always only one parent
    if (r2 < 1) {
        div.branch_upd <- outgoing[[edges$node.child[div.branch]]][1]
    } else {
        div.branch_upd <- outgoing[[edges$node.child[div.branch]]][2]  
    }
  }

  old_len <- edges$length[div.branch]
  new_len <- edges$length[div.branch_upd]

  qr <- qr - log(1/new_len)
  qr <- qr + log(1/old_len)
  
  div.times_upd <- runif(1, nodes$times[edges$node.parent[div.branch_upd]], nodes$times[edges$node.child[div.branch_upd]])
  
  x_next[[1]] <- rates
  x_next[[2]] <- K
  x_next[[3]] <- div.times_upd
  x_next[[4]] <- div.branch_upd

  return(list(x_next = x_next, qr=qr))
}


move_update_N <- function(N_prev, pre) {

  N_upd <- rnorm(1, mean = N_prev, sd = 1)
  qr <- -dnorm(N_upd, mean = N_prev, sd = 1, log=TRUE) ## proposal lh
  qr <- qr + dnorm(N_prev, mean = N_upd, sd = 1, log=TRUE) ## reverse lh
  return(list(N_next = N_upd, qr=qr))
}


move_update_probs <- function(probs, pre) {
  d <- 1e-2
  delta <- runif(1, -d, d)

  which <- sample.int(length(probs), 2, replace=T)

  probs_upd <- probs
  probs_upd[which] <- probs_upd[which] + c(delta, -delta)

  qr <- -log(1/(2*d)) - log(1/((length(probs)**2)))  ## proposal lh
  qr <- qr + log(1/(2*d)) + log(1/((length(probs)**2)))## reverse lh
  
  return(list(probs_next=probs_upd, qr=qr))
}  