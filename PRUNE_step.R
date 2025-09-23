PRUNE_step = function(T, ds) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## function to PRUNE a singly internal node of a tree T
  
  ## INPUT:
  ## T: current tree, with only root
  ## ds: terminal node whose parent is pruned
  
  D = T$leafCount
  
  tmp_parent = T$leaves[[ds]]$parent
  
  tmp_node_L = tmp_parent$children[[1]]$name
  tmp_node_R = tmp_parent$children[[2]]$name
  
  tmp_parent$RemoveChild(tmp_node_L)
  tmp_parent$RemoveChild(tmp_node_R)
  
  Di_star_new = rep(0, length(y))
  for (d in 1:T$leafCount) {
    for (i in T$leaves[[d]]$ix) {
      Di_star_new[i] = d
    }
  }
  
  return(list(T, Di_star_new))
}
