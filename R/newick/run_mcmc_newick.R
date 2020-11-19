library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

output.dir <- "./grad2014" #Name of the output directory to be created

dir.create(file.path(".", output.dir))
setwd(file.path(".", output.dir))

tree <- read.tree(file = "../grad2014.nwk") #Path to newick string goes here
tree <- makeNodeLabel(tree)

pdf(file= paste0("tree.pdf"), width = 5, height = 5)
    plt<-ggtree(tree, ladderize=TRUE) +
                          geom_point() +
                          theme_tree2()
plot(plt)
dev.off()

pre <- structured_coal.preprocess_phylo(tree)

r_mean <- 0 ## growth rate lognormal prior mean
K_mean <- 6 ## carrying capacity rate lognormal prior mean

r_sd <- 2 ## growth rate lognormal prior sd
K_sd <- 0.5 ## carrying capacity  rate lognormal prior sd

prior_i <- function(x) dpois(x, 0, log = TRUE) ### poisson 0 prior

prior_N <- function(x) dlnorm(x, meanlog = K_mean, sdlog = K_sd, log = TRUE)
prior_N.sample <- function() rlnorm(1, meanlog = K_mean, sdlog = K_sd) 

prior_r <- function(x) dlnorm(x, meanlog = r_mean, sdlog = r_sd, log = TRUE) 
prior_r.sample <- function(x) rlnorm(1, meanlog = r_mean, sdlog = r_sd) 

prior_K <- function(x) dlnorm(x, meanlog = K_mean, sdlog = K_sd, log = TRUE)
prior_K.sample <- function() rlnorm(1, meanlog = K_mean, sdlog = K_sd) 

prior_t <- function(x) {
       return(log(1/abs(max(pre$nodes.df$times)-min(pre$nodes.df$times)))) ### Uniform time prior
} 

prior_t.sample <-function() runif(1, min(pre$nodes.df$times), max(pre$nodes.df$times)) ### Uniform time prior

set.seed(1)

o <- outbreaks_infer(tree, prior_i,  prior_N,  prior_N.sample, 
                     prior_r, prior_r.sample,  prior_K,  prior_K.sample,  prior_t,
                     prior_t.sample, 1, n_it=1e6, thinning=2, debug=F)

y <- sapply(o$dims, function(x) x)
n <- sapply(o$para, function(x) x[[1]])
lh <- sapply(o$log_lh, function(x) x)
prior <- sapply(o$log_prior, function(x) x)

branches <- lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(o$para[[x]][c(-1, -2)], function(y) list(br=y[[4]], x=x)) else NA)
times <- unlist(lapply(c(1:length(o$para)), function(x) if(o$dims[[x]] > 0) lapply(o$para[[x]][c(-1, -2)], function(y) y[[3]])))
times.df <- data.frame(times=times)
branches.br <- unlist(lapply(branches, function(x) if(is.na(x)) NA else lapply(x, function(y) y$br)))
branches.x <-  unlist(lapply(branches, function(x) if(is.na(x)) NA else lapply(x, function(y) y$x)))
br.df <- data.frame(x=branches.x, br=branches.br)

x <- c(1:length(y))
df <- data.frame(x=x, y=y, n=n, lh=lh, prior=prior)


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
       geom_point(alpha=0.1,size=0.1)+
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

pdf(file=paste0("times_hist.pdf"), width = 5, height = 5)
plt <- ggplot(times.df, aes(times)) +
         geom_histogram(colour="darkgreen", fill="white", binwidth = 0.01) + 
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


setwd(file.path(".", ".."))