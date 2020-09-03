#'Generate a Newick string representation of a tree corresponding to a realisation of a coalescent process.
#' 
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events / internal nodes.
#' @return a Newick string corresponding to the tree.
#' @export

build_coal_tree <- function(sampling_times, coalescent_times, leaf_names=NULL,node_name_prefix="N", terminate_string = TRUE) {
  coal_times_desc <- coalescent_times[order(-coalescent_times)]
  sam_ord <- order(-sampling_times)
  times_desc <- sampling_times[sam_ord]

  if (!is.null(leaf_names)) {
    leaf_names <- leaf_names[sam_ord]
  }

  if (is.null(leaf_names)) {
    tree_nodes <- seq(1, length(sampling_times))
    tree_nodes <-
    sapply(tree_nodes, function (x)
      return (paste0("S", x)))
  } else {
    tree_nodes <- leaf_names
  }
  extant_entries <- c(1)
  extant_times <- c(times_desc[1])
  
  coal_idx <- 1
  s_idx <- 1
  t <- times_desc[1]
  
  while (coal_idx <= length(coalescent_times)) {
    if (s_idx < length(times_desc) &&
        (times_desc[s_idx + 1] > coal_times_desc[coal_idx])) {
      s_idx <- s_idx + 1
      t <- times_desc[s_idx]
      extant_entries <- c(extant_entries, s_idx)
      extant_times <- c(extant_times, t)
    } else {
      t <- coal_times_desc[coal_idx]
      coal_node1_idx <-
        trunc(runif(1, 1, length(extant_entries) + 1))
      
      coal_node1 <- extant_entries[coal_node1_idx]
      ct1 <- extant_times[coal_node1_idx]
      
      extant_times <- extant_times[-coal_node1_idx]
      extant_entries <- extant_entries[-coal_node1_idx]
      
      br_len_1 <- ct1 - t
      entry1 <- tree_nodes[coal_node1]
      
      coal_node2_idx <-
        trunc(runif(1, 1, length(extant_entries) + 1))
      coal_node2 <- extant_entries[coal_node2_idx]
      ct2 <- extant_times[coal_node2_idx]
      br_len_2 <- ct2 - t
      entry2 <- tree_nodes[coal_node2]
      
      tree_nodes[coal_node2] <- paste0("(", entry1,
                                       ":", br_len_1,
                                       ",",
                                       entry2,
                                       ":", br_len_2, ")",node_name_prefix,coal_idx)
      extant_times[coal_node2_idx] <- t
      coal_idx <- coal_idx + 1
    }
  }
  if (terminate_string) {
    tree_str <- paste0(tree_nodes[extant_entries[1]], ";")
  } else{
     tree_str <- tree_nodes[extant_entries[1]]
  }

  return(tree_str)
}

#'Generate a Newick string representation of a tree corresponding to a realisation of a structured,linage-based coalescent process.
#'
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events / internal nodes.
#' @param leaf_colours leaf colouring
#' @param coalescent_colours a numeric vector where number at j-th position indicates that the j-th coalescent event corresponds to j-th lineage.
#' @param div_times times of lineage divergence events. 
#' @param div_events a numeric vector where number at j-th position indicates that the j-th divergence time corresponds to j-th lineage
#' @param div_from numeric vector identifying the parent of the j-th diverging lineage
#' @return a list with the Newick string corresponding to the tree, and a list of leaf colouring assignments
#' @export

build_coal_tree.structured <- function(sampling_times, coalescent_times, leaf_colours, coalescent_colours, div_times, div_events, div_from, include_div_nodes = TRUE) {

  sam_ord <- order(-sampling_times)
  coal_ord <- order(-coalescent_times)

  coal_times_desc <- coalescent_times[coal_ord]
  times_desc <- sampling_times[sam_ord]

  leaf_colours_desc <- leaf_colours[sam_ord]
  coalescent_colours_desc <- coalescent_colours[coal_ord]

  subtrees <- rep("", length(div_times))

  leaf_subs <- lapply(c(1:length(div_times)), function(x) which(leaf_colours==x))
  node_subs <- lapply(c(1:length(div_times)), function(x) which(coalescent_colours==x))

  if (!include_div_nodes) {
    subtree_MRCA <- sapply(c(1:length(div_times)), function(x) min(coalescent_times[node_subs[[x]]]))
  }

  ### First generate appropriate subtrees.
  for (i in c(1:length(div_times))) {
    
    ### Subset leafs and coalescent times

    sam_subs <- sampling_times[leaf_subs[[i]]]
    coal_subs <- coalescent_times[node_subs[[i]]]

    leaf_names <- sapply(c(1:length(leaf_subs[[i]])), function (x) paste0("S_",LETTERS[i],x)) 
    node_name_prefix <- paste0("N_", LETTERS[i])

    ### add any divergence event times as sampling times to this lineage, mark them with D[j] where j is the number of lineage diverging

    for (j in c(1:length(div_from))) {
      if (div_from[j] == i) {

        if (include_div_nodes) {
          sam_subs <- c(sam_subs, div_times[j])
        } else {
          sam_subs <- c(sam_subs, subtree_MRCA[j])
        }

        leaf_names <- c(leaf_names, paste0("#D_", j))
      } 
    }

    tree <- build_coal_tree(sam_subs, coal_subs, leaf_names=leaf_names, node_name_prefix=node_name_prefix, terminate_string = FALSE)

    if (include_div_nodes && div_times[i] > -Inf) {
        branch_len <- min(coal_subs)-div_times[i]
        tree <- paste0("(",tree,":",branch_len,")","X_",LETTERS[i])
    }

    subtrees[i] <- tree
  }
  
  subtrees.ret <- sapply(subtrees, function (x) paste0(x,";"))

  ### Next build the combined tree
  for(i in c(1:(length(div_events)-1))){
    child <- div_events[i]
    parent <- div_from[i]

    diverging_tree <- subtrees[child]
    parent_tree <- subtrees[parent]

    subtrees[parent] <- gsub(paste0("#D_", i), diverging_tree, parent_tree)
  }

  tree_str <- paste0(subtrees[length(subtrees)], ";")
  return(list(full=tree_str, subtrees=subtrees.ret))
}

#'Plots a structured tree with different lineages distinguished by different colours, and divergence events highlighted.
#'
#' @param sampling_times a tree including divergence events
#' @param n_lineages number of lineages
#' @return a list with the Newick string corresponding to the tree, and a list of leaf colouring assignments
#' @export
plot_structured_tree <- function(tree, n_lineages){

    labs <- c(tree$node.label, tree$tip.label)

    lineages <- lapply(c(1:n_lineages), function (x) labs[grep(paste0("[N,X,S]_",LETTERS[x]), labs)])

    lin_names <- c(sapply(c(1:n_lineages), function (x) paste0("expansion ",x)), "neutral")

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
                    theme_tree2()
    return(plt)
}