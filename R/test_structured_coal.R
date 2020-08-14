library(CaveDive)
library(ape)
library(ggtree)
library(ggplot2)
library(treeio)
library(viridis)

set.seed(1)
    
    sam <- runif(100, 0, 10)

    tmax <- max(sam)

    sam <- sam - tmax
    sam <- sam[order(-sam)]
    
    colours <- trunc(runif(100, 1, 4))

    N <- 100
    
    A <- c(1.1, 0.3)
    
    K <- c(100,150)

    div_times <- c(-25, -40, -Inf)
    div_cols <- c(1, 2, 3)

    rates <- list(function (s) half_log.rate(s, K[1], A[1], div_times[1]), function (s) half_log.rate(s, K[2], A[2], div_times[2]), function (s) constant.rate(s, N))
    rate.ints <- list(function(t,s) half_log.rate.int(t, s, K[1], A[1], div_times[1]), function(t,s) half_log.rate.int(t, s, K[2], A[2],div_times[2]), function(t,s) constant.rate.int(t,s,N))

    co <- structured_coal.simulate(sam, colours, div_times, div_cols, rates, rate.ints)
    print(co)

    tr <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from)
    
    tree <-read.tree(text = tr$full)

    labs <- c(tree$node.label, tree$tip.label)
    n_lineages <- 3

    lineages <- lapply(c(1:n_lineages), function (x) labs[grep(paste0("[N,X,S]_",LETTERS[x]), labs)])

    lin_names <- c("expansion 1", "expansion 2", "neutral")

    lineage_labs <- unlist(lineages)
    membership <- unlist(lapply(c(1:length(lineages)), function (x) rep(lin_names[x], length(lineages[[x]]))))
    type <- sapply(lineage_labs, 
        function (x) if (length(grep("X_", x, value="FALSE"))>0) "Divergence" else if(length(grep("N_", x, value="FALSE"))>0) "Coalescent" else "Sampling")

    ldf <- data.frame(node = nodeid(tree, lineage_labs), lineage = membership, type=type)
    ldf$node <- as.numeric(ldf$node)
    ldf$lineage <- as.factor(ldf$lineage)
    ldf$type <- as.factor(ldf$type)

    tree.full <- full_join(tree, ldf, by = 'node')

    pdf("structured_tree.pdf")
    plt<-ggtree(tree.full, aes(color=lineage), ladderize=TRUE) +
                    geom_point(aes(shape=type, size=type)) +
                    scale_size_manual(values=c(1,4,1)) +
                    scale_shape_manual(values=c(1,8,2)) +
                    theme_tree2()#+ 
                    #scale_colour_viridis(discrete = TRUE, option = "plasma")
    plot(plt)
    dev.off()

    tr.nodiv <- build_coal_tree.structured(sam, co$times, colours, co$colours, div_times, div_cols, co$div_from, include_div_nodes = FALSE)
    tree.nodiv <- read.tree(text = tr.nodiv$full)

    times.nodiv <- node.depth.edgelength(tree.nodiv)
    times.nodiv <- times.nodiv-max(times.nodiv)
    times.nodiv <- times.nodiv[101:length(times.nodiv)]
    times.ord <- order(-times.nodiv)
    times.nodiv <- times.nodiv[times.ord]

    MRCAs.idx <- sapply(c(1:3), function (x) (which(co$colours==x)[which.min(times.nodiv[which(co$colours==x)])])) 
    MRCAs<- sapply(MRCAs.idx, function (x) tree.nodiv$node.label[times.ord[x]])

    pre <- structured_coal.preprocess_phylo(tree.nodiv)
    lh.comp <-  structured_coal.likelihood(pre, MRCAs, div_times, A, K, N)

    print(lh.comp)

