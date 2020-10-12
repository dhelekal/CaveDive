library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

set.seed(12345678)

n <- 2

n_tips <- 100
sam <- runif(n_tips, 0, 0.1)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

colours <- trunc(runif(n_tips, 1, n+2))

N <- 10#rexp(1, rate = 1/100)
K <- rep(N, n)#rexp(n, rate = 1/100)
A <- rexp(n, rate = 0.1)
div_times <- c(-1*runif(n,1,3), -Inf)
div_times <- div_times[order(-div_times)]

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


prior_i <- function(x) dpois(x, 1, log = TRUE)

prior_N <- function(x) dlnorm(x, meanlog = 1, sdlog = 1, log = TRUE)
prior_N.sample <- function() rlnorm(1, meanlog = 1, sdlog = 1) 

prior_r <- function(x) dexp(x, rate=0.1, log=TRUE) 
prior_r.sample <- function(x) rexp(1, rate = 0.1)

prior_K <- function(x) dlnorm(x, meanlog = 1, sdlog = 1, log = TRUE)
prior_K.sample <- function() rlnorm(1, meanlog = 1, sdlog = 1) 

set.seed(0)

o <- infer_outbreaks(tree, prior_i, prior_N, prior_N.sample, prior_r, prior_r.sample, prior_K, prior_K.sample, n_it=1e6, thinning=1, debug=FALSE)

y <- sapply(o$dims, function(x) x)
n <- sapply(o$para, function(x) x[[1]])
r <- sapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) o$para[[x]][[2]][[1]] else NA)
branches <- sapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) o$para[[x]][[2]][[4]] else NA)

x <- c(1:length(y))
df <- data.frame(x=x, y=y, n=n, r=r, branches=branches)
df <- df[1e5:1e6,]

pdf(file=paste0("dim_hist.pdf"), width = 5, height = 5)
plt <- ggplot(df, aes(y)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("N_hist.pdf"), width = 5, height = 5)
plt <- ggplot(df, aes(n)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("r_hist.pdf"), width = 5, height = 5)
plt <- ggplot(df, aes(r)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()


tt <- table(df$branches)
freq <- sapply(c(1:length(tt)), function (i) tt[i]/pre$edges.df$length[as.integer(names(tt))[i]])

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