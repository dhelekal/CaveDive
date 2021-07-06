#' Simulate expansion structured coalescent
#' Population size must be locally smooth around divergence time and must reach 0 at divergence time !!
#' Note: This simulation relies on inverse transform sampling and numerical inverses of Neg.rate.ints  
#'
#' @param sampling_times Nx1 times of leaves.
#' @param colours Nx1 colour assignment of leaves
#' @param div_times Kx1 divergence times for colours, neutral colour must diverge at -Inf. MUST BE IN DESCENDING ORDER
#' @param div_events Kx1 colour assignment for divergence times, i.e. which colour corresponds to which divergence time
#' @param Neg.rates Kx1 List of Functions (t)->(R+). 1/Neg(t) for each colour.
#' @param Neg.rate.ints Kx1 List of Functions (t, s)->(R+). int_{t}^{s} 1/Neg(tau),d tau for each colour
#' @param div.from Kx1 List of integers denoting parents of diverging lineages. If NA parents will be randomised
#' @return A list consisting of the simulate coalescent times \code{coalescent_times} and the log-likelihood of the simulation \code{log_likelihood}.
#' @export
structured_coal.simulate <- function(sampling_times, colours, div_times, div_events, Neg.rates, Neg.rate.ints, div.from=NA) {

    sam_ord <- order(-sampling_times)
    times_desc <- sampling_times[sam_ord]
    colours_desc <- colours[sam_ord]

    n_col <- length(unique(colours_desc))

    lh_ctr <- 0

    if (n_col != length(div_times)) {
        warning("Number of divergence events does not match number of unique colours")
        return(-1)
    }

    extant_lineages <- as.integer(rep(0, n_col))
    future_lineages <- as.integer(sapply(div_events, function (x) sum(colours_desc==x)))

    future_lineages <- future_lineages[order(div_events)]

    coalescent_times <- rep(0, length(times_desc)-1)
    coalescent_cols <- rep(0, length(times_desc)-1)
    div_from <- rep(0, n_col-1)

    t <- 0
    t0 <- times_desc[1]

    div_idx <- 1
    coal_idx <- 1
    sam_idx <- 2

    extant_lineages[colours_desc[1]] <- 1
    future_lineages[colours_desc[1]] <- future_lineages[colours_desc[1]] - 1

    log_lh <- 0

    while (max(future_lineages) > 0 || sum(extant_lineages) > 1) {
        if (max(extant_lineages) < 2) {
            ### If no lineage has more than 1 member continue to next divergence or sampling event
            if(sam_idx <= length(times_desc) && (div_idx > length(div_times) || div_times[div_idx] < times_desc[sam_idx])) { 
                t <- t0 - times_desc[sam_idx]

                which_colour <- colours_desc[sam_idx]

                extant_lineages[which_colour] <- extant_lineages[which_colour] + 1    
                future_lineages[which_colour] <- future_lineages[which_colour] - 1
                sam_idx <- sam_idx + 1

            } else {
                t <- t0 - div_times[div_idx]

                which_div <- div_events[div_idx]
                if (!all(!is.na(extant_lineages)) || extant_lineages[which_div] != 1) {
                    warning("Divergence event registered at this time but lineage MRCA has not been reached. This is likely a bug")
                    warning(paste0("lineage diverging node count: ", extant_lineages[which_div]))
                    warning(paste0("time: ", t))
                    warning(paste0("time increment: ", s))
                    warning(paste0("lineage diverging: ", which_div))
                    warning(paste0("extant_lineages vector: ", extant_lineages, "\n"))
                    return(-1)
                }
                ### choose which lineage to diverge from
                extant_lineages[which_div] <- extant_lineages[which_div] - 1 
                    
                i <- NA
                if (is.na(div.from)) {
                    non_zero_l <- which(extant_lineages > 0)
                    i <- non_zero_l[runif(1,1,length(non_zero_l)+1)]
                    log_lh <- log_lh + log(1/length(non_zero_l))
                } else {
                    i <- div.from[div_idx]
                }

                extant_lineages[i] <- extant_lineages[i] + 1    
                div_from[div_idx] <- i
 
                div_idx <- div_idx + 1     
            }
        } else {
            comb_ns <- sapply(extant_lineages, function (x) choose(x, 2))
            opt <- -1
            if (sam_idx > length(times_desc) && div_idx > (length(div_times)-1)) {
                opt <- 1
                s <- Inf
            } else if(sam_idx > length(times_desc)){
                opt <- 2
                s <- t0 - t - div_times[div_idx]
            } else if(div_idx > (length(div_times)-1)){
                opt <- 3
                s <- t0 - t - times_desc[sam_idx]
            } else {
                opt <- 4
                s <- t0 - t - max(times_desc[sam_idx], div_times[div_idx])
            }

            if (t > t+s) {
                print(paste0("sam_idx: ",sam_idx))
                print(paste0("div_idx: ",div_idx))
                print(paste0("t0: ",t0))
                print(paste0("t: ",t))
                print(paste0("t+s: ",t+s))
                print(paste0("opt: ", opt))
                warning("next time step must be greater than current time")
                return(-1)
            }

            if (!all(comb_ns >= 0)){
                warning("combination numbers must be positive")
                return(-1)
            }

            ### Decide if an coalescent event happens in this interval
            p_coal <- inhomogenous_exp.prob(rate_int_sum(Neg.rate.ints, comb_ns, n_col),t,s)
            r <- runif(1,0,1)
            if (r<=p_coal) { 
                ### choose which event happens and compute waiting time.

                inv_int <- function(t, Fs) {
                    func <- function(s) rate_int_sum(Neg.rate.ints, comb_ns, n_col)(t,s)
                    dfun <- function(s) sum(sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rates[[x]](t+s)))
                    return(inv_rates(func, Fs, dfun))

                }

                w_t <- inv_t_inhomogenous_exp_conditional(rate_int_sum(Neg.rate.ints, comb_ns, n_col),
                                                            function(t, s) inv_int(t, s),
                                                            1,
                                                            t,
                                                            s)
                if(is.na(w_t)) {
                    warning("Numeric inverse transform failed: failed to find initial bisection interval")
                    warning(paste0("time: ", t))
                    warning(paste0("time increment: ", s))
                    warning(paste0("lineage diverging: ", which_div))
                    warning(paste0("p_coal: ", p_coal))
                    warning(paste0("extant_lineages vector: ", extant_lineages, "\n"))
                    warning(paste0("opt: ", opt))
                    warning(paste0("comb_ns vector: ", comb_ns, "\n"))
                    return(-1)
                } 

                if(is.infinite(w_t)) {
                    warning("Numeric inverse transform failed: infinity returned")
                    warning(paste0("time: ", t))
                    warning(paste0("time increment: ", s))
                    warning(paste0("lineage diverging: ", which_div))
                    warning(paste0("p_coal: ", p_coal))
                    warning(paste0("extant_lineages vector: ", extant_lineages, "\n"))
                    warning(paste0("opt: ", opt))
                    warning(paste0("comb_ns vector: ", comb_ns, "\n"))
                    return(-1)
                } 
                ### Compute rates for each lineage
                rates <- sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rates[[x]](t+w_t))
                i <- choose_reaction(rates)

                log_upd <- inhomogenous_exp.loglh(function(s) (sum(sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*Neg.rates[[x]](s)))),
                                                          rate_int_sum(Neg.rate.ints, comb_ns, n_col), t, w_t) - log(comb_ns[i]) + log(rates[i]) - log(sum(rates))

                log_lh <- log_lh + log_upd
                t <- t + w_t
                extant_lineages[i] <- extant_lineages[i] - 1
                coalescent_times[coal_idx] <- t0 - t
                coalescent_cols[coal_idx] <- i
                coal_idx <- coal_idx + 1
            } else {
                lh_ctr <- lh_ctr + 1
                log_upd <- log(1 - p_coal)
                log_lh <- log_lh + log_upd
                ### determine whether sampling or divergence event comes next
                ### If no lineage has more than 1 member continue to next divergence or sampling event
                if(sam_idx <= length(times_desc) && (div_idx > (length(div_times)-1) || div_times[div_idx] < times_desc[sam_idx])) {
                    t <- t0 - times_desc[sam_idx]

                    which_colour <- colours_desc[sam_idx]

                    extant_lineages[which_colour] <- extant_lineages[which_colour] + 1    
                    future_lineages[which_colour] <- future_lineages[which_colour] - 1
                    sam_idx <- sam_idx + 1
                } else {

                    which_div <- div_events[div_idx]
                    if (!all(!is.na(extant_lineages)) || extant_lineages[which_div] != 1) {
                        warning("Divergence event registered at this time but lineage MRCA has not been reached. This is likely a bug")
                        warning(paste0("lineage diverging node count: ", extant_lineages[which_div]))
                        warning(paste0("time: ", t))
                        warning(paste0("time increment: ", s))
                        warning(paste0("lineage diverging: ", which_div))
                        warning(paste0("p_coal: ", p_coal))
                        warning(paste0("extant_lineages vector: ", extant_lineages, "\n"))
                        warning(paste0("comb_ns vector: ", comb_ns, "\n"))
                        return(-1)
                    }
                    t <- t0 - div_times[div_idx]
                    ### choose which lineage to diverge from
                    extant_lineages[which_div] <- extant_lineages[which_div] - 1 
                    
                    i <- NA
                    if (all(is.na(div.from))) {
                        non_zero_l <- which(extant_lineages > 0)
                        i <- non_zero_l[runif(1,1,length(non_zero_l)+1)]
                        log_lh <- log_lh + log(1/length(non_zero_l))
                    } else {
                        i <- div.from[div_idx]
                    }

                    extant_lineages[i] <- extant_lineages[i] + 1    

                    div_from[div_idx] <- i

                    div_idx <- div_idx + 1     
                } 
            }
        }
    }
    return(list(times=coalescent_times, colours=coalescent_cols, div_from=div_from, log_lh=log_lh))
}

#' Preprocess phylogeny to speed up likelihood computation
#'
#' @param phy phylogeny
#' @param order_edges_by_node_label should edge indexing be ordered by child node label. Default: TRUE
#' @return preprocessed phylogeny
#' @export
preprocess_phylo <- function(phy, order_edges_by_node_label=TRUE){
    stopifnot("phy must be an ape phylogeny"= class(phy) == "phylo")
    stopifnot("phylogeny must be labeled" = all(!is.na(phy$node.label))&&all(!is.na(phy$tip.label)))

    labs <- c(phy$node.label, phy$tip.label)
    nodes <- nodeid(phy, labs)
    is_tip <- c(rep(FALSE, length(phy$node.label)), rep(TRUE, length(phy$tip.label)))

    times <- node.depth.edgelength(phy)
    times <- times - max(times)
    times <- times[nodes]

    t_min <- min(times)
    t_max <- max(times)

    edges.parent <- phy$edge[,1]
    edges.child <- phy$edge[,2]
    edges.len <- phy$edge.length

    nodes.df <- data.frame(id=nodes, times=times, is_tip=is_tip, lab=labs)
    nodes.df <- nodes.df[order(nodes.df$id), ]

    if (!order_edges_by_node_label) {
        edges.df <- data.frame(node.parent=edges.parent, node.child=edges.child, length=edges.len, id=c(1:length(edges.parent)))
        edges.df <- edges.df[order(edges.df$id), ]
    } else {
        edges.df <- data.frame(node.parent=edges.parent, node.child=edges.child, length=edges.len)
        edges.df <- edges.df[order(nodes.df$lab[edges.df$node.child]), ]
        edges.df$id <- c(1:length(edges.parent))
    }

    stopifnot("branch lengths must be nonegative" = all(edges.df$length >= 0))

    edges.outgoing <- lapply(nodes.df$id, function (x) edges.df$id[which(edges.df$node.parent==x)])
    edges.outgoing <- lapply(edges.outgoing, function (x) if (length(x) > 0) x else NA)
    edges.incoming <- lapply(nodes.df$id, function (x) edges.df$id[which(edges.df$node.child==x)])
    edges.incoming <- lapply(edges.incoming, function (x) if (length(x) > 0) x else NA)

    root <- which(is.na(edges.incoming))

    clades.list <- lapply(nodes.df$id[which(nodes.df$is_tip==FALSE)], function(x) extract.clade(phy, x))

    cl <- clade_lookup(phy, root, nodes.df, edges.df, edges.outgoing)

    node_times_clade_lookup<-nodes.df$times[cl$clade_node_lookup$node_list] 
    tip_times_clade_lookup<-nodes.df$times[cl$clade_tip_lookup$tip_list] 

    return(structure(list(phy=phy,
                nodes.df = nodes.df,
                edges.df = edges.df,
                clades.list = clades.list, 
                n_tips = length(phy$tip.label), 
                incoming = edges.incoming,
                outgoing = edges.outgoing,
                t_min = t_min, 
                t_max = t_max,
                root_idx = root,
                clade_tip_lookup=cl$clade_tip_lookup,
                clade_node_lookup=cl$clade_node_lookup,
                node_times_clade_lookup=node_times_clade_lookup,
                tip_times_clade_lookup=tip_times_clade_lookup), class="preprocessedPhy"))
}

clade_lookup <- function(phy, root_idx, nodes.df, edges.df, edges.outgoing) {
    n_node <- length(which(nodes.df$is_tip==FALSE))
    n_tips <- length(which(nodes.df$is_tip==TRUE))

    clade_bounds_tips <- data.frame(lo=rep(NA,n_node+n_tips), hi=rep(NA,n_node+n_tips))
    tip_list <- rep(NA, n_tips)

    tl <- clade_tip_recursion(1, n_tips, root_idx, tip_list, clade_bounds_tips, nodes.df, edges.df, edges.outgoing, n_tips, phy)

    clade_bounds_nodes <- data.frame(lo=rep(NA,n_node+n_tips), hi=rep(NA,n_node+n_tips))
    node_list <- rep(NA, n_node)

    nl <- clade_node_recursion(1, n_node, root_idx, node_list, clade_bounds_nodes, nodes.df, edges.df, edges.outgoing, n_tips, phy)

    stopifnot(all(!is.na(nl$node_list)))
    stopifnot(all(!is.na(tl$tip_list)))

    return(list(clade_tip_lookup=list(bounds=tl$clade_bounds, tip_list=tl$tip_list),
                clade_node_lookup=list(bounds=nl$clade_bounds, node_list=nl$node_list)))
}

clade_tip_recursion <- function(index_lo, index_hi, id, tip_list, clade_bounds, nodes.df, edges.df, edges.outgoing, n_tips, phy) {
    cb <- clade_bounds
    cb$hi[id] <- index_hi
    cb$lo[id] <- index_lo
    if (nodes.df$is_tip[id]) {
        stopifnot("index error" = index_hi==index_lo)
        tl <- tip_list
        tl[index_lo] <- id
        return(list(tip_list=tl, clade_bounds=cb))
    } else {
        next_v <- sapply(edges.outgoing[[id]], function (i) edges.df$node.child[i])
        next_v <- next_v[order(next_v)]

        if (!nodes.df$is_tip[next_v[1]]) {
            count_left <- length(extract.clade(phy, next_v[1])$tip.label)
        } else {
            count_left <- 1
        }

        if (!nodes.df$is_tip[next_v[2]]) {
            count_right <- length(extract.clade(phy, next_v[2])$tip.label)
        } else {
            count_right <- 1
        }

        c <- clade_tip_recursion(index_lo, index_hi-count_right, next_v[1], tip_list, cb, nodes.df, edges.df, edges.outgoing, n_tips, phy)
        c <- clade_tip_recursion(index_lo+count_left, index_hi, next_v[2], c$tip_list, c$clade_bounds, nodes.df, edges.df, edges.outgoing, n_tips, phy)
        return(c)
    }
}

clade_node_recursion <- function(index_lo, index_hi, id, node_list, clade_bounds, nodes.df, edges.df, edges.outgoing, n_tips, phy) {
    stopifnot("index error" = !nodes.df$is_tip[id])
    cb <- clade_bounds
    cb$hi[id] <- index_hi
    cb$lo[id] <- index_lo

    nl <- node_list
    nl[index_lo] <- id

    next_v <- sapply(edges.outgoing[[id]], function (i) edges.df$node.child[i])
    next_v <- next_v[order(next_v)]

    if (!nodes.df$is_tip[next_v[1]]) {
        count_left <- length(extract.clade(phy, next_v[1])$node.label)
    } else {
        count_left <- 0
    }

    if (!nodes.df$is_tip[next_v[2]]) {
        count_right <- length(extract.clade(phy, next_v[2])$node.label)
    } else {
        count_right <- 0
    }

    if (!nodes.df$is_tip[next_v[1]]) {
        c <- clade_node_recursion(index_lo+1, index_hi-count_right, next_v[1], nl, cb, nodes.df, edges.df, edges.outgoing, n_tips, phy)
        nl <- c$node_list
        cb <- c$clade_bounds
    }

    if (!nodes.df$is_tip[next_v[2]]){
        c <- clade_node_recursion(index_lo+count_left+1, index_hi, next_v[2], nl, cb, nodes.df, edges.df, edges.outgoing, n_tips, phy)
        nl <- c$node_list
        cb <- c$clade_bounds
    }

    return(list(node_list=nl, clade_bounds=cb))
}

#' @export
print.preprocessedPhy <- function(x, ...) {
    cat(paste("\nPreprocessed phylogeny with", x$n_tips, "tips\n\n"))
    cat("  ",paste("\nThe most recent sampling t_max time is: ", x$t_max), "\n", sep="")
    cat("  ",paste("\nThe time of the MRCA t_min is: ", x$t_min), "\n", sep="")
    cat("  ",paste("\nRoot index root_idx: ", x$root_idx), "\n", sep="")
    cat("  ",paste("\nincoming: a list of length: ", length(x$incoming)), "\n", sep="")
    cat("  ",paste("\nedges.df: dataframe with names: ", colnames(x$edges.df)), sep="")
    cat("  ",paste("\nnodes.df: dataframe with names: ", colnames(x$nodes.df)), sep="")
    cat("  ",paste("\nclades.list: a list of length: ", length(x$clades.list)), "\n", sep="")
    cat("  ",paste("\nincoming: a list of length: ", length(x$incoming)), "\n", sep="")
    cat("  ",paste("\noutgoing: a list of length: ", length(x$outgoing)), "\n", sep="")
}


#' Compute likelihood for preprocessed phylogeny 
#'
#' @param phylo.preprocessed preprocessed phylogeny
#' @param div.MRCA.nodes labels of MRCA nodes of diverging lineages
#' @param div_times absolute divergence times
#' @param diverging.rates growth rates for diverging lineages
#' @param diverging.sizes asymptotic size for diverging lineages
#' @param neutral.size size of neutral phylogeny
#' @param times.extracted list returned by extract_lineage_times. Used instead of div.MRCA.nodes to avoid recomputing the tree partition.
#' @return Log-likelihood
#' @export
structured_coal.likelihood <- function(phylo.preprocessed, div.MRCA.nodes, div_times, diverging.rates, diverging.sizes, neutral.size, times.extracted=NA, type="Sat"){
    n_tips <- phylo.preprocessed$n_tips

    MRCA.idx <- nodeid(phylo.preprocessed$phy, div.MRCA.nodes)-n_tips
    k_div <- length(div.MRCA.nodes)
    log_lh <- 0

    if(is.na(times.extracted)) {
        times.extracted <- extract_lineage_times_native(phylo.preprocessed, div.MRCA.nodes, div_times)
    }
    partition_counts <- times.extracted$partition_counts

    if(times.extracted$empty_tips) {
        #warning("MRCA selection produced subtrees with empty tips.")
        log_lh <- -Inf
    } else {

        t_max <- max(sapply(c(1:length(times.extracted$sam.times)), function (x) max(times.extracted$sam.times[[x]])))

        if (k_div > 1) range <- c(1:(k_div-1)) else range <- c()

        for (i in range){
            if (type == "Log-Exp") {
                    log_lh <- log_lh + logexp_coalescent_loglh(times.extracted$sam.times[[i]], times.extracted$coal.times[[i]], div_times[i], diverging.rates[i], diverging.sizes[i], t_max)
                } else if (type=="Sat") {
                    log_lh <- log_lh + sat_coalescent_loglh(times.extracted$sam.times[[i]], times.extracted$coal.times[[i]], div_times[i], diverging.rates[i], diverging.sizes[i], t_max)
                } else {
                    warning(paste0("Unrecognised likelihood option: ", type))
                    return(NA)
                }
            
        }
        log_lh <- log_lh + coalescent_loglh(times.extracted$sam.times[[k_div]], times.extracted$coal.times[[k_div]], neutral.size, t_max)
        #log_lh <- log_lh -lgamma(length(div_times))#+ log(1/factorial(length(div_times)-1))
    }

    return(list(log_lh = log_lh, sam.times = times.extracted$sam.times, coal.times = times.extracted$coal.times, partition_counts=partition_counts))
}

extract_lineage_times_native <- function(phylo.preprocessed, div.MRCA.nodes, div_times) {
        stopifnot("div_times must be in ascending order"=order(-div_times)==c(1:length(div_times)))

        res <- extract_partition_times_fast(nodeid(phylo.preprocessed$phy, div.MRCA.nodes), div_times, 
                                           phylo.preprocessed$node_times_clade_lookup, 
                                           phylo.preprocessed$tip_times_clade_lookup, 
                                           phylo.preprocessed$clade_node_lookup$bounds, 
                                           phylo.preprocessed$clade_tip_lookup$bounds)
        return(res)
}

#' Extracts Lineage times (but faster) from a parent phylogeny based on provided divergence MRCA nodes and times. 
#'
#' @param phylo.preprocessed preprocessed phylogeny
#' @param div.MRCA.nodes labels of MRCA nodes of diverging lineages
#' @param div_times absolute divergence times
#' @param return_partitions if true, return a list named 'partitions' containing tip labels for each partition
#' @return A list containing lists of sampling time vectors sam.times and coalescent time vectors coal.times for each lineage.
#' @export
extract_lineage_times <- function(phylo.preprocessed, div.MRCA.nodes, div_times, return_partitions=FALSE) {
    times.ord <- order(-div_times)
    k_div <- length(div_times)

    div_from <- rep(NA, k_div)

    coal.times <- list()
    sam.times <- list()

    if (return_partitions) {
        partitions <- list()
    } else {
        partitions <- NA
    }

    empty_tips <- FALSE
    partition_counts <- rep(0, k_div)

    n_tips <- phylo.preprocessed$n_tips
    nodes.df <- phylo.preprocessed$nodes.df

    MRCA.idx <- nodeid(phylo.preprocessed$phy, div.MRCA.nodes)

    clade_node_lookup <- phylo.preprocessed$clade_node_lookup
    clade_tip_lookup <- phylo.preprocessed$clade_tip_lookup
    
    for (i in c(1:k_div)) {
        ii <- times.ord[i]

        node_bounds <- list(lo=clade_node_lookup$bounds$lo[MRCA.idx[ii]], hi=clade_node_lookup$bounds$hi[MRCA.idx[ii]])
        tip_bounds <-list(lo=clade_tip_lookup$bounds$lo[MRCA.idx[ii]],hi=clade_tip_lookup$bounds$hi[MRCA.idx[ii]])

        subset_bounds_lo_n <- c()
        subset_bounds_hi_n <- c()

        subset_bounds_lo_t <- c()
        subset_bounds_hi_t <- c()

        leaf.times <- c()
        for (j in which(c(1:k_div) < i)) {
            jj <- times.ord[j]
            bounds_sub_tip <- list(lo=clade_tip_lookup$bounds$lo[MRCA.idx[jj]],hi=clade_tip_lookup$bounds$hi[MRCA.idx[jj]])

            if (bounds_sub_tip$lo >= tip_bounds$lo && 
                bounds_sub_tip$hi <= tip_bounds$hi && is.na(div_from[j])){

                div_from[j] <- i 
                leaf.times <- c(leaf.times, div_times[jj])

                bounds_sub_node <- list(lo=clade_node_lookup$bounds$lo[MRCA.idx[jj]], hi=clade_node_lookup$bounds$hi[MRCA.idx[jj]])

                subset_bounds_lo_t<-c(subset_bounds_lo_t, bounds_sub_tip$lo)
                subset_bounds_hi_t<-c(subset_bounds_hi_t, bounds_sub_tip$hi)

                subset_bounds_lo_n<-c(subset_bounds_lo_n, bounds_sub_node$lo)
                subset_bounds_hi_n<-c(subset_bounds_hi_n, bounds_sub_node$hi)

                partition_counts[ii] <- partition_counts[ii] - 1
            }
        }

        total_tip_len <- tip_bounds$hi - tip_bounds$lo - 
                         sum(subset_bounds_hi_t-subset_bounds_lo_t) + 1 - length(subset_bounds_hi_t)
        total_node_len <- node_bounds$hi - node_bounds$lo - 
                         sum(subset_bounds_hi_n-subset_bounds_lo_n) + 1 - length(subset_bounds_hi_n)

        tip_idx <- rep(Inf,total_tip_len)
        node_idx <- rep(Inf,total_node_len)

        arr_idx <- 1
        s_idx <- tip_bounds$lo

        if (length(subset_bounds_lo_t)>0) {
            ord.t <- order(subset_bounds_lo_t)
            subset_bounds_hi_t <- subset_bounds_hi_t[ord.t]
            subset_bounds_lo_t <- subset_bounds_lo_t[ord.t]
            for (k in c(1:length(subset_bounds_lo_t))) {
                seg_len <- subset_bounds_lo_t[k]-s_idx
                if (seg_len > 0) {
                    tip_idx[arr_idx:(arr_idx+seg_len-1)] <- clade_tip_lookup$tip_list[s_idx:(subset_bounds_lo_t[k]-1)]
                    arr_idx<-arr_idx+seg_len
                }
                s_idx <- subset_bounds_hi_t[k]+1
            }
            seg_len <- tip_bounds$hi-s_idx
            if (seg_len > 0) {
                tip_idx[arr_idx:(arr_idx+seg_len)] <- clade_tip_lookup$tip_list[s_idx:tip_bounds$hi]
            }
        } else {
            tip_idx[arr_idx:total_tip_len] <- clade_tip_lookup$tip_list[tip_bounds$lo:tip_bounds$hi]
        }

        arr_idx <- 1
        s_idx <- node_bounds$lo
        
        if (length(subset_bounds_lo_t)>0){
            ord.n <- order(subset_bounds_lo_n)
            subset_bounds_hi_n <- subset_bounds_hi_n[ord.n]
            subset_bounds_lo_n <- subset_bounds_lo_n[ord.n]
        
            for (k in c(1:length(subset_bounds_lo_n))) {
                seg_len <-  subset_bounds_lo_n[k]-s_idx
                if (seg_len > 0) {
                    node_idx[arr_idx:(arr_idx+seg_len-1)] <- clade_node_lookup$node_list[s_idx:(subset_bounds_lo_n[k]-1)]
                    arr_idx<-arr_idx+seg_len
                }
                s_idx <- subset_bounds_hi_n[k]+1
            }
            seg_len <- node_bounds$hi-s_idx
            if (seg_len > 0) {
                node_idx[arr_idx:(arr_idx+seg_len)] <- clade_node_lookup$node_list[s_idx:node_bounds$hi]
            }
        } else {
            node_idx[arr_idx:total_node_len] <- clade_node_lookup$node_list[node_bounds$lo:node_bounds$hi]
        }
        tip.times <- nodes.df$times[tip_idx]
        node.times <- nodes.df$times[node_idx]

        if (return_partitions){
            partitions[[ii]] <- nodes.df$lab[tip_idx][order(-tip.times)]
        }
        if (length(tip.times)==0) {
            empty_tips <- TRUE
        }

        leaf.times <- c(leaf.times, tip.times)

        sam.times[[ii]] <- leaf.times[order(-leaf.times)]
        coal.times[[ii]] <- node.times[order(-node.times)]

        partition_counts[ii] <- partition_counts[ii] + length(sam.times[[ii]])
    }
    return(list(sam.times=sam.times, 
                coal.times=coal.times, 
                empty_tips=empty_tips, 
                partition_counts=partition_counts, 
                partitions=partitions))
}

choose_reaction <- function(rates) {

    if (sum(rates==Inf) > 1) {
        warning("More than one entry equal to infinity. This should not happen")
    } else if (!all(rates != Inf)) {
        which_inf <- which(rates==Inf)
        return(which_inf[runif(1,1,length(which_inf+1))])
    } else {

        total_rate <- sum(rates)
        r <- runif(1, 0, total_rate)

        i <- 0 
        rate_sum <- total_rate
        while (rate_sum > r) {
            i <- i+1
            rate_sum <- rate_sum-rates[i]
        }
        return(i)
    }
}

rate_int_sum <- function(rate_ints, comb_ns, n_col) {
    out <- function(t,s){
        rate_vec <- sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else comb_ns[x]*rate_ints[[x]](t, s))
        rate_vec2 <- sapply(c(1:n_col), function (x) if(comb_ns[x] == 0) 0 else rate_ints[[x]](t, s))
        if(!all(rate_vec >= 0)) {
            warning("Negative rates encountered")
            warning(paste0("rate vec: ", rate_vec2, "\n"))
            warning(paste0("t0 at error: ", t))
            warning(paste0("s at error: ", s))
        }
        return(sum(rate_vec))
    }
    return(out)
}

