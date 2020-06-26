inv_t_conditional_exp <- function(u,rate,delta_t){
  return((-1/rate)*log(1-u*(1-exp(-rate*delta_t))))
}

cond_exp.lh <- function(rate, u, delta_t){
  return(rate*(1/(1-exp(-rate*delta_t)) - u))
}

exp.lh <- function(rate, t){
  return(rate*exp(-rate*t))
}

exp.prob <- function(rate, t){
  return(1-exp(-rate*t))
}


poi_0.lh <- function(rate, t){
  return(exp(-rate*t))
}

inhomogenous_exp.lh <- function(rate, rate.int, t, s){
  return(rate(t+s)*exp(-rate.int(t,s)))
}

inhomogenous_exp.prob <- function(rate.int, t, s) {
  return(1-exp(-rate.int(t,s)))
}

inv_t_inhomogenous_exp_conditional <- function(rate.int, rate.inv_int , exp.rate, t, s){
  u <- runif(1,0,1)
  Q <- inhomogenous_exp.prob(rate.int, t, s)
  wt <- (-1/exp.rate)*log(1-u*Q)
  return(rate.inv_int(t , wt))
}

inhomogenous_poi_0.lh <- function(rate.int, t, s){
  return(exp(-rate.int(t, s)))
}

build_coal_tree <- function(sampling_times, coalescent_times){
  coal_times_desc <- coalescent_times[order(-coalescent_times)]
  times_desc <- sampling_times[order(-sampling_times)]
  tree_nodes <- seq(1,length(sampling_times))
  tree_nodes <- sapply(tree_nodes, function (x) return (paste0("S", x)))
  extant_entries <- c(1)
  extant_times <- c(times_desc[1])
  
  coal_idx <- 1 
  s_idx <- 1
  t <- times_desc[1]
  
  while (coal_idx <= length(coalescent_times)) {
    if (s_idx<length(times_desc) && (times_desc[s_idx+1] > coal_times_desc[coal_idx])){
      s_idx <- s_idx+1  
      t <- times_desc[s_idx]
      extant_entries <- c(extant_entries, s_idx)
      extant_times <- c(extant_times, t)
    } else {
      t <- coal_times_desc[coal_idx]
      coal_node1_idx <- trunc(runif(1,1,length(extant_entries)+1))
      
      coal_node1 <- extant_entries[coal_node1_idx]
      ct1 <- extant_times[coal_node1_idx]
      
      extant_times <- extant_times[-coal_node1_idx]
      extant_entries <- extant_entries[-coal_node1_idx]
      
      br_len_1 <- ct1-t
      entry1 <- tree_nodes[coal_node1]
      
      coal_node2_idx <- trunc(runif(1,1,length(extant_entries)+1))
      coal_node2 <- extant_entries[coal_node2_idx]
      ct2 <- extant_times[coal_node2_idx]
      br_len_2 <- ct2-t
      entry2<-tree_nodes[coal_node2]
      
      tree_nodes[coal_node2] <- paste0("(",entry1,
                                       ":",br_len_1,
                                       ",",
                                       entry2,
                                       ":",br_len_2,")")
      extant_times[coal_node2_idx] <- t
      coal_idx <- coal_idx+1
    }
  }
  tree_str<-paste0(tree_nodes[extant_entries[1]], ";")
  return(tree_str)
}