###list of: N, div_event_1, ... div_event_N
###div_event list of r, K, x_0, branch
offset <- 2

prop.sampler <- function (x_prev, i_prev, pre, para.initialiser, initialiser.log_lh, fn_log_J, fn_log_J_inv, pop_scale=10){
  p <- 1 
  r <- runif(1,0,2) 

  jump_scale <- pop_scale/10

  if (r < p) { ## transdimensional
  x_prop <- transdimensional.sampler(x_prev, i_prev, pre, para.initialiser, initialiser.log_lh, fn_log_J, fn_log_J_inv)
  x_next <- x_prop$x_next
  i_next <- x_prop$i_next

  log_J <- x_prop$log_J
  qr <- x_prop$qr 
} else { ## within-model
x_prop <- within_model.sampler(x_prev, i_prev, pre, jump_scale)
x_next <- x_prop$x_next
i_next <- i_prev

log_J <- 0
qr <- x_prop$qr 
}
return(list(x_next=x_next,i_next=i_next, qr=qr, log_J=log_J))
}

transdimensional.sampler <- function(x_prev, i_prev, pre, para.initialiser, initialiser.log_lh, fn_log_J, fn_log_J_inv, fixed_move=NA) {
  if (is.na(fixed_move)){
    which_move <- sample.int(2, size=1)
  } else {
    which_move <- fixed_move
  }

  log_J <- 0
  qr <- 0
  if (which_move==1) { ### increase dim
  i_next <- i_prev + 1
  x_next <- x_prev
  N_prev <- x_prev[[1]]

  x_next[[i_next + offset]] <- para.initialiser(N_prev)

    old_probs <- x_prev[[2]]
    which_split <- sample.int((i_prev+1),1)
    u <- runif(1, 0, old_probs[which_split])

    stopifnot("probs vector lengths not equal i_prev+1"=length(old_probs)==(i_prev+1))

    new_probs <- c(old_probs[-(i_prev+1)], u, old_probs[(i_prev+1)])

    if(which_split < (i_prev+1)) {
      new_probs[which_split] <- new_probs[which_split] - u
    } else {
      new_probs[which_split+1] <- new_probs[which_split+1] - u
    }
  x_next[[2]] <- new_probs

  qr <- qr - log(1/(i_prev+1))
  qr <- qr - log(1/old_probs[which_split])
  qr <- qr - initialiser.log_lh(x_next[[i_next + offset]], N_prev) 
  qr <- qr + log(1/(i_prev+1))
  qr <- qr + log(1/(i_prev+1))

  if(abs(sum(x_next[[2]]) - 1) > 1e-8) warning("prob sum error")

    log_J <- fn_log_J(i_prev, x_prev, x_next)
   } else { ### decrease dim
    x_next <- x_prev
    N_prev <- x_prev[[1]]
    if (i_prev > 0) {
      i_next <- i_prev - 1
      which_elem <- sample.int(i_prev, size=1)
      x_next <- x_next[-(which_elem+offset)]

      stopifnot("x_prev[[2]] vector lengths not equal i_prev+1"=length(x_prev[[2]])==(i_prev+1))

      which_prob <- x_prev[[2]][which_elem]
      x_next[[2]] <- x_next[[2]][-which_elem]

      which_merge <- sample.int(i_prev, size=1)

      x_next[[2]][which_merge] <- x_next[[2]][which_merge] + which_prob
      
      qr <- qr - log(1/i_prev) ## proposal remove model 
      qr <- qr - log(1/i_prev) ## proposal merge with this probaility
      qr <- qr + initialiser.log_lh(x_prev[[(which_elem+offset)]], N_prev) ## reverse lh of adding that model
      qr <- qr + log(1/i_prev) ## reverse pick same prob for split 
      qr <- qr + log(1/x_next[[2]][which_merge]) ##pick the same split from uniform
      
      if(abs(sum(x_next[[2]]) - 1) > 1e-8) warning("prob sum error")
      if(length(x_next[[2]]) != (i_next+1)) warning("prob length error")
      log_J <- fn_log_J_inv(i_prev, x_prev, x_next, which_elem)
  } else {
    i_next <- i_prev
  }
}
return(list(x_next=x_next, i_next=i_next, log_J=log_J, qr=qr))
}

within_model.sampler <- function(x_prev, i_prev, pre, scale, fixed_move=NA, fixed_index=NA) {

  if (is.na(fixed_move)){
    if (i_prev > 0) {
      which_move <- sample(c(1,2,3,4), size=1)
    } else {
      which_move <- 3
    }
  } else {
    which_move <- fixed_move
  }

  if (which_move==1) {
    x_prop <- make_move(x_prev, i_prev, pre, move_update_mid.time, scale, fixed_index)
    x_next <- x_prop$x_next
    qr <- x_prop$qr
  } else if(which_move==2) {
    x_prop <- make_move(x_prev, i_prev, pre, move_update_branch2, scale, fixed_index)
    x_next <- x_prop$x_next
    qr <- x_prop$qr
  } else if(which_move==3) {
    N_prev <- x_prev[[1]]
    N_prop <- move_update_N(N_prev, pre, scale)
    N_next <- N_prop$N_next
    qr <- N_prop$qr

    x_next <- x_prev
    x_next[[1]] <- N_next 
  } else if(which_move==4) {
    probs_prev <- x_prev[[2]]
    probs_prop <- move_update_probs(probs_prev, pre, fixed_index)
    probs_next <- probs_prop$probs_next
    qr <- probs_prop$qr

    x_next <- x_prev
    x_next[[2]] <- probs_next 
  }
  return(list(x_next=x_next,qr=qr))
}

make_move <- function(x_prev, i_prev, pre, move, scale, fixed_index) {

  if(is.na(fixed_index)) {
    which_model <- offset+sample.int(i_prev, size=1)
  } else {
    which_model <- offset+fixed_index
  }
  mdl_prev <- x_prev[[which_model]] 
  mdl_prop <- move(mdl_prev, pre, scale)
  mdl_next <- mdl_prop$x_next
  qr <- mdl_prop$qr
  x_next <- x_prev
  x_next[[which_model]] <- mdl_next
  return(list(x_next=x_next, qr=qr))
}

move_update_mid.time <- function(x_prev, pre, scale) { ### update mid.time

  x_next <- vector(mode = "list", length = length(x_prev))

  mid.time <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

  mid.time_upd <- rnorm(1, mean=mid.time, sd=scale) 
  K_upd  <- rlnorm(1, meanlog=log(K), sdlog=0.15)

  qr <- -dnorm(mid.time_upd, mean=mid.time, sd=scale, log=TRUE) - dlnorm(K_upd, meanlog=log(K), sdlog=0.15, log=TRUE) ## proposal lh
  qr <- qr + dnorm(mid.time, mean=mid.time_upd, sd=scale, log=TRUE) + dlnorm(K, meanlog=log(K_upd), sdlog=0.15, log=TRUE) ## reverse lh
  
  x_next[[1]] <- mid.time_upd 
  x_next[[2]] <- K_upd
  x_next[[3]] <- div.times
  x_next[[4]] <- div.branch
  return(list(x_next=x_next, qr=qr))
}

move_update_branch2 <- function(x_prev, pre, scale) { ### update branch
  x_next <- vector(mode = "list", length = length(x_prev))

  mid.time <- x_prev[[1]]
  K <- x_prev[[2]]
  div.times <- x_prev[[3]]
  div.branch <- x_prev[[4]]

  edges <- pre$edges.df
  nodes <- pre$nodes.df

  outgoing <- pre$outgoing
  incoming <- pre$incoming

  root <- pre$root_idx

  delta_t <- rnorm(1, mean=0, sd=scale)
  qr <- -dnorm(delta_t, mean=0, sd=scale, log=TRUE)
  qr <- qr + dnorm(delta_t, mean=0, sd=scale, log=TRUE)

  ## Root acts as reflecting boundary for time. 
  if ((div.times+delta_t) < nodes$times[root]) {
    div.times_upd <- ((div.times+delta_t) - nodes$times[root])
  } else {
    div.times_upd <- (div.times+delta_t)
  } 

  direction <- sign(delta_t)
  curr_br <- div.branch
  while(!(nodes$times[edges$node.child[curr_br]] > div.times_upd &&
    nodes$times[edges$node.parent[curr_br]] < div.times_upd) && !is.na(curr_br)) {
    if (direction < 0) {
      if (edges$node.parent[curr_br] == root) { ### traversing root, switch direction
        direction <- -direction
        curr_br <- outgoing[[root]][which(outgoing[[root]]!=curr_br)]
        stopifnot("Invalid Branch Selected"=length(curr_br)==1)
      } else {
        curr_br <- incoming[[edges$node.parent[curr_br]]]
        qr <- qr + log(1/2)
      } 
    } else {
     r <- sample(2,1)
     curr_br <- outgoing[[edges$node.child[curr_br]]][r]
     stopifnot("Invalid Branch Selected"=length(curr_br)==1)
     qr <- qr - log(1/2)
    }
  }

  if (is.na(curr_br)) {
    qr <- -Inf
  }

  x_next[[1]] <- mid.time
  x_next[[2]] <- K
  x_next[[3]] <- div.times_upd
  x_next[[4]] <- curr_br

  return(list(x_next = x_next, qr=qr))
}


move_update_N <- function(N_prev, pre, scale) {

  N_upd <- rlnorm(1, meanlog=log(N_prev), sdlog=0.15)
  qr <- -dlnorm(N_upd, meanlog=log(N_prev), sdlog=0.15, log=TRUE) ## proposal lh
  qr <- qr + dlnorm(N_prev, meanlog = log(N_upd), sdlog=0.15, log=TRUE) ## reverse lh
  return(list(N_next = N_upd, qr=qr))
}


move_update_probs <- function(probs, pre, fixed_index = NA) {

  k <- length(probs)
  if (!all(!is.na(fixed_index))){
    which <- sample.int(k, 2, replace=T)
  } else {
    which <- fixed_index
  }

  i1 <- which[1]
  i2 <- which[2]

  probs_upd <- probs
  u <- runif(1,0,probs[i1])

  probs_upd[i1] <- probs_upd[i1] - u
  probs_upd[i2] <- probs_upd[i2] + u

  qr <- -log(1/probs[i1]) - log(1/((k**2)))  ## proposal lh
  qr <- qr + log(1/probs_upd[i2]) + log(1/((k**2)))## reverse lh
  
  return(list(probs_next=probs_upd, qr=qr))
}  


### Proposal likelihood functions used for testing
prop_lh <- function(x_prev, i_prev, x_next, i_next, pre, initialiser.log_lh, scale=10) {
  lh <- 0
    if (i_next == i_prev) {
        lh <- lh + dlnorm(x_next[[1]], meanlog= log(x_prev[[1]]), sdlog= 0.15, log=TRUE)
        p_next <- x_next[[2]]
        p_prev <- x_prev[[2]]
        lh <- lh + log(1/((length(p_next)**2)))

    p_diff <- p_prev - p_next

    p_idx <- which(p_diff > 1e-8)
    if (length(p_idx) > 0) {
      lh <- lh + (log(1/p_prev[p_idx[1]]))
    }

    if (i_next > 0){ 
      mdl_lh <- sum(sapply(c(1:i_next), function (j) model_lh2(x_prev[[offset+j]], x_next[[offset+j]], pre, scale)))
    } else { 
      mdl_lh <- 0
    }
    lh <- lh + mdl_lh
  } else if (i_next > i_prev) {
      # Need to find which event was added -- use that events are uniquely identified by branch + time
      # No need for good performance or efficiency this is a test function
    if (i_prev > 0 ) {
      prev_branches <- sapply(c(1:i_prev), function (j) x_prev[[offset+j]][[4]])
      prev_times <- sapply(c(1:i_prev), function (j) x_prev[[offset+j]][[3]])
      unique_idx <- which(sapply(c(1:i_next), 
        function (j) length(which(x_next[[offset+j]][[4]] == prev_branches & abs(x_next[[offset+j]][[3]] - prev_times) < 1e-6)))==0)
    } else {
      unique_idx <- 1
    }
    lh <- lh + initialiser.log_lh(x_next[[offset+unique_idx]], x_prev[[1]])
    lh <- lh + log(1/i_next)

    p_next <- x_next[[2]]
    p_prev <- x_prev[[2]]

    p_diff <- p_prev - p_next[-(unique_idx)]

    p_idx <- which(p_diff > 1e-8)
    lh <- lh + (log(1/p_prev[p_idx[1]]))
  } else {
    lh <- lh + log((1/i_prev**2))
  }
  return(lh)
}

model_lh2 <- function(mdl_prev, mdl_next, pre, scale) {
  lh <- 0
  lh <- lh + dnorm(mdl_next[[1]], mean=mdl_prev[[1]], sd=scale, log=TRUE) +  ### Rates and Carrying capacity
  dlnorm(mdl_next[[2]], meanlog=log(mdl_prev[[2]]), sdlog=0.2, log=TRUE) 

  br_next <- mdl_next[[4]]
  br_prev <- mdl_prev[[4]]

  if(is.na(br_next) || is.na(br_prev)) {
    lh <- -Inf
  } else {
    t_next <- mdl_next[[3]]
    t_prev <- mdl_prev[[3]]

    edges <- pre$edges.df
    nodes <- pre$nodes.df

    outgoing <- pre$outgoing
    incoming <- pre$incoming

    root <- pre$root_idx
    traversed_1 <- 0
    traversed_2 <- 0

    root_traversed <- FALSE

    if (br_prev != br_next) {
        ### if 1 previous branch comes before next branch in tree, if 2 next branch comes before previous branch in tree
      if (nodes$times[edges$node.parent[br_prev]] > nodes$times[edges$node.parent[br_next]]){
        lower_br <- br_prev
        upper_br <- br_next
        dir <- TRUE
      } else {
        lower_br <- br_next
        upper_br <- br_prev
        dir <- FALSE
      }
      current_br <- lower_br
      while(current_br != upper_br) { 
        if (edges$node.parent[current_br] == root) {
          tmp <- upper_br
          upper_br <- outgoing[[root]][which(outgoing[[root]]!=current_br)]
          current_br <- tmp
          root_traversed <- TRUE
          dir <- !dir
        } else {
          current_br <- incoming[[edges$node.parent[current_br]]]
          if (dir == FALSE) {
            traversed_1 <- traversed_1 + 1
          }
        }
      }
      
      lh <- lh + log(1/2)*traversed_1

    }
    if (root_traversed) {
      lh <- lh + dnorm((nodes$times[root]-t_prev) + (nodes$times[root]-t_next), mean=0, sd=scale, log=TRUE)
    } else {
      lh <- lh + dnorm(t_next-t_prev, mean=0, sd=scale, log=TRUE)
    }
  }
  return(lh)
}