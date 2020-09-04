library(CaveDive)
library(ggplot2)
library(viridis)
library(rmutil)
library(ape)

set.seed(123456)

n_tips <- 100
sam <- runif(n_tips, 0, 10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

colours <- trunc(runif(n_tips, 1, 4))

n <- 2

N <- rexp(1, rate = 1/100)
A <- rexp(n, rate = 1/20)
K <- rexp(n, rate = 1/100)
div_times <- c(-1*runif(n,20,80), -Inf)

div_cols <- c(1, 2, 3)
rates <- list(function (s) sat.rate(s, K[1], A[1], div_times[1]), function (s) sat.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
rate.ints <- list(function(t,s) sat.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) sat.rate.int(t, s, K[2], A[2], div_times[2]), function(t,s) constant.rate.int(t,s,N))
co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)
#print(co)

tree.div.str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)
tree.plt <- plot_structured_tree(tree.div, 3)

pdf(file="tree_strucutred.pdf")
plot(tree.plt)
dev.off()

tree.str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
tree <- read.tree(text = tree.str$full)

pre <- structured_coal.preprocess_phylo(tree)

total_branch_len <- sum(pre$edges.df$length)

# x[[1]] rates
# x[[2]] K
# x[[3]] neutral size
# x[[4]] div.times
# x[[5]] div.branch

## Select random branch and a time along the branch 
inner_branches <- pre$edges.df$id[which(pre$edges.df$node.child>n_tips)] 
rand_br <- inner_branches[runif(n, 1, (length(inner_branches)+1))]
rand_times <- sapply(rand_br, function (x) runif(1, pre$nodes.df$times[pre$edges.df$node.parent[x]],
                                                    pre$nodes.df$times[pre$edges.df$node.child[x]]))
x_0 <- list()
x_0[[1]] <- rexp(n, rate = 1/10)
x_0[[2]] <- rexp(n, rate = 1/100)
x_0[[3]] <- rexp(1, rate = 1/100)
x_0[[4]] <- rand_times
x_0[[5]] <- rand_br

# x[[1]] rates
# x[[2]] K
# x[[3]] neutral size
# x[[4]] div.times
# x[[5]] div.branch
log_lh <- function(x){
  rates <- x[[1]]
  K <- x[[2]]
  N <- x[[3]]
  div.times <- x[[4]]
  div.branch <- x[[5]]

  prior_rates <- sum(dexp(rates, rate = 1/10, log = TRUE))
  prior_K <- sum(dexp(K, rate = 1/100, log = TRUE))
  prior_N <- dexp(N, rate = 1/100, log = TRUE)


  if (all(K > 0) && all(rates > 0) && all(N > 0) && all(!is.na(div.branch))){
      MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)  
      if (all(!is.na(MRCAs))) {
        prior_br <- 0 #sum(log(pre$edges.df$length[div.branch]/total_branch_len))
        lh <- 0 #structured_coal.likelihood(pre, MRCAs, div.times, rates, K, N, type="Sat")$log_lh
        prior <- prior_rates + prior_K + prior_N + prior_br
        lh <- lh + prior
      } else {
        lh <- -Inf
      }
  } else {
    lh <- -Inf
  }
  return(lh)
}

prop.sampler <-function (x_prev, it){
  rates <- x_prev[[1]]
  K <- x_prev[[2]]
  N <- x_prev[[3]]
  div.times <- x_prev[[4]]
  div.branch <- x_prev[[5]]

  x_next <- vector(mode = "list", length = length(x_prev))

  if (it %% 2) {
    rates_upd <- rnorm(length(rates), mean=rates, 1) 
    K_upd <- rnorm(length(K), mean=K, 1)
    N_upd <- rnorm(length(N), mean=N, 1)

    x_next[[1]] <- rates_upd
    x_next[[2]] <- K_upd
    x_next[[3]] <- N_upd
    x_next[[4]] <- div.times
    x_next[[5]] <- div.branch

  } else {
    div.times_upd <- rnorm(length(div.times), mean=div.times, 10)
    div.branch_upd <- sapply(c(1:length(div.branch)), function (i) select_br(pre, div.branch[i], div.times[i], div.times_upd[i]))

    x_next[[1]] <- rates
    x_next[[2]] <- K
    x_next[[3]] <- N
    x_next[[4]] <- div.times_upd
    x_next[[5]] <- div.branch_upd
  }
  return(x_next)
}


select_br <- function(pre, div.br, div.time, div.time_upd) {
  edges <- pre$edges.df
  nodes <- pre$nodes.df
  br <- div.br

  if (div.time_upd-div.time > 0) {
    outgoing <- pre$outgoing
    while (!is.na(br) && nodes$times[edges$node.child[br]] < (div.time_upd)) {
      r <- runif(1,1,2)
      br <- outgoing[[edges$node.child[br]]][r]
    }
  } else {
    incoming <- pre$incoming
    while (!is.na(br) && nodes$times[edges$node.parent[br]] > (div.time_upd)) {
      br <- incoming[[edges$node.parent[br]]]
    }
  }
  return(br)
}

branch_log_lh <- function(src, dest) {
  edges <- pre$edges.df
  nodes <- pre$nodes.df
  out <- 0
  if (nodes$times[edges$node.child[src]] < nodes$times[edges$node.child[dest]]){
    ### Check if direction is forwards. If backwards likelihood is one.
    ### iterate backwards exponentiating two at each branch point
    br <- dest 
    incoming <- pre$incoming
    while (!is.na(br) && br != src) {
      br <- incoming[[edges$node.parent[br]]]
      out <- out + log(1/2) ### loops infinitely
    }
  }
  return(out)
}

proposal.cond_lh <- function(x_cand, x_prev, it){

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

  out <- 0

  if (it %% 2) {
    out <- sum(dnorm(rates_cand, mean=rates_prev, 1, log=TRUE)) + 
           sum(dnorm(K_cand, mean=K_prev, 1, log=TRUE)) +
           sum(dnorm(N_cand, mean=N_prev, 1, log=TRUE))
  } else {
    out <- sum(dnorm(div.times_cand, mean=div.times_prev, 10, log=TRUE)) +
           sum(sapply(c(1:length(div.branch_cand)), function (x) branch_log_lh(div.branch_prev[x], div.branch_cand[x])))
           

  }
  return(out)
}

n_it <- 100000
burn_in <- 1

o <- run_mcmc(log_lh, proposal.cond_lh, prop.sampler, x_0, n_it, FALSE)
#o.df <- as.data.frame(t(as.data.frame(o)))
branches <- sapply(c(1:length(o)), function (i) o[[i]][[5]])
h<-table(branches[1,burn_in:n_it])
freq<-sapply(as.integer(unlist(dimnames(h))), function (i) h[toString(i)]/pre$edges.df$length[i]) 