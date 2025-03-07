library(phyclust)

get_trees_from_ms <- function(ms_output) {
  i <- 3
  tree_index <- c()
  while(i<=length(ms_output)) {
    tree_index <- c(tree_index,i)
    i <- i+2
  }
  trees <- read.tree(text = ms_output[tree_index])
  return(trees)
}

make_ultrametric <- function(tree) {
  max_depth <- max(node.depth.edgelength(tree))
  for (i in 1:length(tree$tip.label)) {
    tip_depth <- node.depth.edgelength(tree)[i]
    tree$edge.length[tree$edge[,2] == i] <- tree$edge.length[tree$edge[,2] == i] + (max_depth - tip_depth)
  }
  return(tree)
}

# Function to simulate gene trees using ms
simulate_gene_trees <- function(species_tree, num_trees = 100) {
  # set node labels to null so that branching.times names the node numbers
  species_tree$node.label <- NULL
  if(any(branching.times(species_tree)<0)) {
    species_tree <- make_ultrametric(species_tree)
  }
  
  # Extract branching times and order of coalescent events
  branching_times <- branching.times(species_tree)
  sorted_events <- sort(branching_times)
  
  # Define population splits for ms
  num_tips <- length(species_tree$tip.label)
  ms_commands <- character()
  
  node_lookup <- c()
  
  for (i in seq_along(sorted_events)) {
    node <- as.numeric(names(sorted_events)[i])
    time <- sorted_events[i]
    
    # Find descendant nodes of the current internal node
    descendants <- species_tree$edge[species_tree$edge[, 1] == node, 2]
    
    # Convert node numbers to population IDs
    desc1 <- ifelse(descendants[1] <= num_tips, descendants[1], node_lookup[descendants[1]])
    desc2 <- ifelse(descendants[2] <= num_tips, descendants[2], node_lookup[descendants[2]])
    
    # Merge the two descendant populations into a new population ID
    ms_commands <- c(ms_commands, sprintf("-ej %f %d %d", time, desc1, desc2))
    node_lookup[node] <- desc2
    # getMRCA(species_tree,c(desc1,desc2))

  }
  
  # Construct ms command
  ms_command <- sprintf("ms %d %d -T -I %d %s %s",
                        num_tips, num_trees, num_tips,
                        paste(rep(1, num_tips), collapse = " "),
                        paste(ms_commands, collapse = " "))
  
  ms_opts <- paste("-T -I",num_tips,paste(rep(1, num_tips), collapse = " "),paste(ms_commands, collapse = " "))
  ms_out <- ms(nsam=num_tips,nreps = num_trees,opts = ms_opts)
  trees <- get_trees_from_ms(ms_out)

  # put original tip labels onto tree
  for (i in seq_along(trees)) {
    tip_labels <- trees[[i]]$tip.label
    trees[[i]]$tip.label <- species_tree$tip.label[as.integer(gsub("s([0-9]+)","\\1",tip_labels))]
  }
  return(trees)
  # return(ms_opts)
}


# Example usage
# species_tree_newick <- "(A:1, ((B:1,C:1):1, (D:1,E:1):1):1);"
# species_tree1 <- read.tree(text = species_tree_newick)
# 
# gene_trees1 <- simulate_gene_trees(species_tree1)
# 
# species_tree_2 <- read.tree("~/Research/genome_assembly/manuscript/trees/astral4_BS10.tre")
# 
# # species_tree_2$edge.length
# 
# gene_trees2 <- simulate_gene_trees(species_tree_2)
# plot(species_tree_2)
# edgelabels()
# 
# gt2 <- read.tree("~/Research/genome_assembly/manuscript/trees/merged.treefile")

# species_tree$tip.label


