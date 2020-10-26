library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)
library(logitnorm)

set.seed(123456789)

n_tips <- 30
sam <- runif(n_tips, 0, 0.1)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

N <- rlnorm(1, meanlog = 5, sdlog = 1) 

co <- homogenous_coal.simulate(sam, N)

tree.str <- build_coal_tree(sam, co$coalescent_times)
tree <- read.tree(text = tree.str)
pdf(file="tree_neutral.pdf")
plot(tree)
dev.off()
pre <- structured_coal.preprocess_phylo(tree)

prior_i <- function(x) dpois(x, 1, log = TRUE)

prior_N <- function(x) dlnorm(x, meanlog = 5, sdlog = 1, log = TRUE)
prior_N.sample <- function() rlnorm(1, meanlog = 5, sdlog = 1) 

prior_r <- function(x) dlnorm(x, meanlog = 0, sdlog = 0.5, log = TRUE) 
prior_r.sample <- function(x) rlnorm(1, meanlog = 0, sdlog = 0.5) 

prior_K <- function(x) dlnorm(x, meanlog = 5, sdlog = 1, log = TRUE)
prior_K.sample <- function() rlnorm(1, meanlog = 5, sdlog = 1) 

set.seed(0)

o <- infer_outbreaks(tree, prior_i, prior_N, prior_N.sample, prior_r, prior_r.sample, prior_K, prior_K.sample, n_it=1e6, thinning=5, debug=F)

y <- sapply(o$dims, function(x) x)
n <- sapply(o$para, function(x) x[[1]]) 
lh <- sapply(o$log_lh, function(x) x)
prior <- sapply(o$log_prior, function(x) x)

branches <- lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(o$para[[x]][-1], function(y) list(br=y[[4]], x=x)) else NA)
times <- unlist(lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(o$para[[x]][-1], function(y) y[[3]])))
times.df <- data.frame(times=times)
branches.br <- unlist(lapply(branches, function(x) if(is.na(x)) NA else lapply(x, function(y) y$br)))
branches.x <-  unlist(lapply(branches, function(x) if(is.na(x)) NA else lapply(x, function(y) y$x)))
br.df <- data.frame(x=branches.x, br=branches.br)

evt <- lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(o$para[[x]][-1], function(y) list(r=y[[1]], K=y[[2]], t=y[[3]], br=y[[4]], x=x)) else NA)
evt.r <- unlist(lapply(evt, function(x) if(is.na(x)) NA else lapply(x, function(y) y$r)))
evt.K <- unlist(lapply(evt, function(x) if(is.na(x)) NA else lapply(x, function(y) y$K)))
evt.t <- unlist(lapply(evt, function(x) if(is.na(x)) NA else lapply(x, function(y) y$t)))
evt.br <- unlist(lapply(evt, function(x) if(is.na(x)) NA else lapply(x, function(y) y$br)))
evt.x <-  unlist(lapply(evt, function(x) if(is.na(x)) NA else lapply(x, function(y) y$x)))
evt.df <- data.frame(r=evt.r, K=evt.K, t=evt.t, br=evt.br, x=evt.x)

x <- c(1:length(y))
df <- data.frame(x=x, y=y, n=n, lh=lh, prior=prior)

root_MRCA <- pre$phy$node.label[pre$nodes.df$id[which.min(pre$nodes.df$times)]-n_tips]
root_div <- -Inf
z <- o$para[[2e5]]
z_mrca <- c(sapply(pre$edges.df$node.child[sapply(z[-1], function(x) x[[4]])],
                 function(x) if (x > n_tips) pre$phy$node.label[x-n_tips] else NA), root_MRCA)
z_times <-c(sapply(z[-1], function(x) x[[3]]), root_div)

l_times <- extract_lineage_times(pre, z_mrca, z_times)

sizes <-c(sapply(z[-1], function(x) x[[2]]), z[[1]])

z_lh <- sapply(c(1:length(sizes)), )

save.image()

png(file="trace_lh.png", width=600, height=600)
plt <- ggplot(df, aes(x=x, y=lh)) +
       geom_line(alpha = 0.3)+
       theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()

png(file="trace_prior.png", width=600, height=600)
plt <- ggplot(df, aes(x=x, y=prior)) +
       geom_line(alpha = 0.3)+
       theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()

png(file="trace_branch.png", width=800, height=800)
plt <- ggplot(br.df, aes(x=x, y=br)) +
       geom_point(alpha=0.1,size=0.1) +
       theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()

png(file="trace_dim.png", width=600, height=600)
plt <- ggplot(df, aes(x=x, y=y)) +
       geom_line(alpha = 0.3)+
       theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()

png(file="trace_n.png", width=600, height=600)
plt <- ggplot(df, aes(x=x, y=n)) +
       geom_line(alpha = 0.3)+
       theme_bw() + theme(aspect.ratio=1, legend.position = "bottom")
plot(plt)
dev.off()

pdf(file=paste0("r_hist.pdf"), width = 5, height = 5)
plt <- ggplot(evt.df, aes(r)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("K_hist.pdf"), width = 5, height = 5)
plt <- ggplot(evt.df, aes(K)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         xlim(c(0,200)) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("times_hist.pdf"), width = 5, height = 5)
plt <- ggplot(times.df, aes(times)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("dim_hist.pdf"), width = 5, height = 5)
plt <- ggplot(df, aes(y)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("N_hist.pdf"), width = 5, height = 5)
plt <- ggplot(df, aes(n)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 1) + 
         xlim(c(0,200)) + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

pdf(file=paste0("branch_bar.pdf"), width = 10, height = 10)
tt.br <- table(br.df$br)
freq <- sapply(c(1:length(tt.br)), function (i) tt.br[i])

aux.df <- data.frame(x=names(tt.br), y = freq)

plt <- ggplot(aux.df, aes(x=x,y=y)) +
         geom_bar(stat="identity", fill="steelblue") + 
         theme(aspect.ratio=1)
plot(plt)
dev.off()

tt <- table(branches.br)
freq <- sapply(c(1:length(tt)), function (i) tt[i])

pdf(file="tree_freq.pdf", width = 5, height = 5)
    labs <- c(tree$node.label, tree$tip.label)
    tip <- c(rep("1",length(tree$node.label)), rep("2", length(tree$tip.label)))
    ids <- nodeid(tree, labs)
    id_freq <- sapply(ids, function (i) if(pre$nodes.df$is_tip[i]) NA else if (is.na(freq[paste0(pre$incoming[[i]])])) 0 else freq[paste0(pre$incoming[[i]])])

    ldf <- data.frame(node = ids, frequency = id_freq, tip=tip)
    tree.full <- full_join(tree, ldf, by = 'node')

    plt<-ggtree(tree.full, aes(color=frequency, linetype=tip), ladderize=TRUE) +
                          geom_point() +
                          scale_linetype(c("solid","dashed","dotted"), na.value = "blank") +
                          scale_size_manual(values=c(1)) +
                          scale_color_viridis() +
                          theme_tree2()
plot(plt)
dev.off()
