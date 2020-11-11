library(CaveDive)
library(rmutil)
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(viridis)

set.seed(12345678)

poi_rate <- 3
concentration <- 3

n_tips <- 200
sam <- runif(n_tips, 0, 10)
sam <- sam - max(sam)
sam <- sam[order(-sam)]

r_mean <- 1
r_sd <- 1

K_mean <- 5
K_sd <- 0.5


time_var <- 40
time_mean <- 80

time_shape <- (time_mean**2)/(time_var**2)
time_rate <- time_mean/(time_var**2)

sim <- outbreaks_simulate(poi_rate, concentration, sam, r_mean, r_sd, K_mean, K_sd, time_rate, time_shape)

co <- sim$co
tip_cols <- sim$tip_colours
div_times <- sim$div_times
div_cols <- sim$div_cols


tree.div.str <- build_coal_tree.structured(sam, co$times, tip_cols, co$colours, div_times, div_cols, co$div_from, include_div_nodes = TRUE)
tree.div <- read.tree(text = tree.div.str$full)
tree.plt <- plot_structured_tree(tree.div, sim$n_exp)

pdf(file="tree_structured.pdf")
plot(tree.plt)
dev.off()

tree.str <- build_coal_tree.structured(sam, co$times, tip_cols, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE, aux_root = FALSE)
tree <- read.tree(text = tree.str$full)
pre <- structured_coal.preprocess_phylo(tree)

prior_i <- function(x) dpois(x, 3, log = TRUE)

prior_N <- function(x) dlnorm(x, meanlog = K_mean, sdlog = K_sd, log = TRUE)
prior_N.sample <- function() rlnorm(1, meanlog = K_mean, sdlog = K_sd) 

prior_r <- function(x) dlnorm(x, meanlog = r_mean, sdlog = r_sd, log = TRUE) 
prior_r.sample <- function(x) rlnorm(1, meanlog = r_mean, sdlog = r_sd) 

prior_K <- function(x) dlnorm(x, meanlog = K_mean, sdlog = K_sd, log = TRUE)
prior_K.sample <- function() rlnorm(1, meanlog = K_mean, sdlog = K_sd) 

prior_t <- function(x) {
       out <- dgamma(-(x - min(sam)), shape=time_shape, scale = 1/time_rate, log=TRUE)
       if(!all(out != Inf)) print("Err")
       return(out)
}
prior_t.sample <-function() min(sam)-rgamma(1, shape=time_shape, scale = 1/time_rate)


set.seed(0)

o <- outbreaks_infer(tree, prior_i,  prior_N,  prior_N.sample, 
                     prior_r, prior_r.sample,  prior_K,  prior_K.sample,  prior_t,
                     prior_t.sample, concentration, n_it=1e6, thinning=1, debug=F)

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

B_set <- nodeid(tree,tree$node.label[grep("N_B", tree$node.label)]) 
B_root<- B_set[which.min(pre$nodes.df$times[B_set])]
B_root.edge <- pre$incoming[[B_root]]

A_set <- nodeid(tree,tree$node.label[grep("N_A", tree$node.label)]) 
A_root<- A_set[which.min(pre$nodes.df$times[A_set])]
A_root.edge <- pre$incoming[[A_root]]

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
       geom_point(alpha=0.1,size=0.1)+
       geom_hline(yintercept = B_root.edge, colour="orange", alpha=1, linetype = "longdash") + 
       geom_hline(yintercept = A_root.edge, colour="red", alpha=1, linetype = "longdash")
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
         geom_vline(xintercept = B_root.edge, colour="orange", linetype = "longdash") + 
         geom_vline(xintercept = A_root.edge, colour="red", linetype = "longdash") + 
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
