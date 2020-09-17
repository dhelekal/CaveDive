library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

set.seed(1234567)

source("./proposal.R")

n <- 1

n_tips <- 100
sam <- runif(n_tips, 0, 10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

colours <- trunc(runif(n_tips, 1, n+2))

N <- 100#rexp(1, rate = 1/100)
K <- rep(N, n)#rexp(n, rate = 1/100)
A <- c(0.1) #rexp(n, rate = 1/20)
div_times <- c(-1*runif(n,20,50), -Inf)

div_cols <- c(1:(n+1))

rates <- lapply(c(1:n), function(i) return(function (s) sat.rate(s, K[i], A[i], div_times[i])))
rates[[n+1]] <- function (s) constant.rate(s, N)

rate.ints <- lapply(c(1:n), function(i) return(function (t, s) sat.rate.int(t, s, K[i], A[i], div_times[i])))
rate.ints[[n+1]] <- function(t, s) constant.rate.int(t, s, N)

co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)

tree.div.str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)
tree.plt <- plot_structured_tree(tree.div, n+1)

pdf(file="tree_structured.pdf")
plot(tree.plt)
dev.off()

tree.str <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE, aux_root = FALSE)
tree <- read.tree(text = tree.str$full)

pre <- structured_coal.preprocess_phylo(tree)
total_branch_len <- sum(pre$edges.df$length)
root_MRCA <- pre$phy$node.label[pre$nodes.df$id[which.min(pre$nodes.df$times)]-n_tips]
root_div <- -Inf

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
x_0[[1]] <- rexp(n, rate = 1/20)
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

  prior_rates <- sum(dexp(rates, rate = 1, log = TRUE))
  prior_K <- sum(dexp(K, rate = 1/100, log = TRUE))
  prior_N <- dexp(N, rate = 1/100, log = TRUE)

  if (all(K > 0) &&
      all(rates > 0) &&
      all(N > 0) && 
      all(!is.na(div.branch)) && 
      all(div.times>pre$nodes.df$times[pre$edges.df$node.parent[div.branch]]) &&
      all(div.times<pre$nodes.df$times[pre$edges.df$node.child[div.branch]]))
  {
      MRCAs <- sapply(pre$edges.df$node.child[div.branch], function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA)
      MRCAs <- c(MRCAs, root_MRCA)
      div.times <- c(div.times, root_div)
      if (all(!is.na(MRCAs))) {
        prior_br <- 0#sum(log(pre$edges.df$length[div.branch]/total_branch_len))
        lh <- 0#structured_coal.likelihood(pre, MRCAs, div.times, rates, K, N, type="Sat")$log_lh
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

n_it <- 3e6
burn_in <- 1e5

set.seed(1)

o <- run_mcmc(log_lh, 
  function (x, x_given) prop.cond_log_lh(x, x_given, pre), 
  function (x_prev) prop.sampler(x_prev, pre), x_0, n_it, FALSE)

marginals <- list()
names <- list()
offset <- 0
### rate marginals
for(i in c(1:n)) {
  names[[i+offset]] <- paste0("rate_",i)
  marginals[[i+offset]] <- sapply(o, function(x) x[[1]][[i]])
}

offset <- offset + n

### K marginals
for(i in c(1:n)) {
  names[[i+offset]] <- paste0("K_",i)
  marginals[[i+offset]] <- sapply(o, function(x) x[[2]][[i]])
}

offset <- offset + n

### N marginal

names[[1+offset]] <- "N"
marginals[[1+offset]] <- sapply(o, function(x) x[[3]])

offset <- offset + 1

### time marginals

for(i in c(1:n)) {
  names[[i+offset]] <- paste0("times_",i)
  marginals[[i+offset]] <- sapply(o, function(x) x[[4]][[i]])
}

offset <- offset + n

### branch marginals

for(i in c(1:n)) {
  names[[i+offset]] <- paste0("branches_",i)
  marginals[[i+offset]] <- sapply(o, function(x) x[[5]][[i]])
}

o.df <- as.data.frame(marginals)
colnames(o.df) <- names 
o.df <- o.df[burn_in:n_it, ]

### Filter auxilliary branch entries
### First locate auxilliary branch
#r_idx <- which(pre$phy$node.label=="R")+100
#br_idx <- pre$outgoing[[r_idx]][1]
#o.df <- o.df[!(o.df$branches_1==br_idx),]

for(i in c(1:n)) {
  pdf(file=paste0("branch_hist_",i,".pdf"), width = 5, height = 5)
  plt <- ggplot(o.df, aes_string(x=paste0("branches_",i))) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
  plot(plt)
  dev.off()

  pdf(file=paste0("branch_bar_",i,".pdf"), width = 5, height = 5)
  tt <- table(o.df[[paste0("branches_",i)]])
  freq <- sapply(c(1:length(tt)), function (i) tt[i]/pre$edges.df$length[as.integer(names(tt))[i]])

  aux.df <- data.frame(x=names(tt), y = freq)

  plt <- ggplot(aux.df, aes(x=x,y=y)) +
         geom_bar(stat="identity", fill="steelblue") + 
         theme(aspect.ratio=1)
  plot(plt)
  dev.off()

  pdf(file=paste0("time_hist_",i,".pdf"), width = 5, height = 5)
  plt <- ggplot(o.df, aes_string(x=paste0("times_",i))) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 0.01) + 
         geom_vline(xintercept = div_times[i], colour="orange", linetype = "longdash") + 
         theme(aspect.ratio=1)
  plot(plt)
  dev.off()

  pdf(file=paste0("rate_hist_",i,".pdf"), width = 5, height = 5)
  plt <- ggplot(o.df, aes_string(x=paste0("rate_",i))) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 0.01) + 
         geom_vline(xintercept = A[i], colour="orange", linetype = "longdash") + 
         theme(aspect.ratio=1)
  plot(plt)
  dev.off()

  pdf(file=paste0("K_hist_",i,".pdf"), width = 5, height = 5)
  plt <- ggplot(o.df, aes_string(x=paste0("K_",i))) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 0.01) + 
         geom_vline(xintercept = K[i], colour="orange", linetype = "longdash") + 
         theme(aspect.ratio=1)
  plot(plt)
  dev.off()
}

pdf(file="N_hist.pdf", width = 5, height = 5)
plt <- ggplot(o.df, aes(x=N)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 0.01) + 
         geom_vline(xintercept = N, colour="orange", linetype = "longdash") + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()


pdf(file="tree_freq.pdf", width = 5, height = 5)
    labs <- c(tree$node.label, tree$tip.label)
    tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
    ids <- nodeid(tree, labs)
    id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) NA else if (is.na(freq[paste0(pre$incoming[[i]])])) 0 else freq[paste0(pre$incoming[[i]])])

    ldf <- data.frame(node = ids, frequency = id_freq, tip=tip)
    tree.full <- full_join(tree, ldf, by = 'node')

    plt<-ggtree(tree.full, aes(color=frequency, linetype=tip), ladderize=TRUE) +
                          geom_point() +
                          scale_linetype(c("solid","dashed"), na.value = "blank") +
                          scale_size_manual(values=c(1)) +
                          scale_color_viridis() +
                          theme_tree2()
plot(plt)
dev.off()

pdf(file="tree_split.pdf", width = 5, height = 5)
    
    labs <- c(tree$node.label, tree$tip.label)
    ids <- nodeid(tree, labs)
    id_freq <- sapply(ids, function (i) if(is.na(pre$incoming[[i]])) 3 else pre$which_half[pre$incoming[[i]]])

    ldf <- data.frame(node = ids, frequency = unlist(id_freq))
    tree.full <- full_join(tree, ldf, by = 'node')

    plt<-ggtree(tree.full, aes(color=frequency), ladderize=TRUE) +
                          geom_point() +
                          scale_size_manual(values=c(1)) +
                          #scale_color_viridis() +
                          theme_tree2()
plot(plt)
dev.off()