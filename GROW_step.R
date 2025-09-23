GROW_step = function(T, y, x_j, ds, j, threshold) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## function to GROW node of a tree T
  
  ## INPUT:
  ## T: current tree
  ## y: response variable of the type y = (1, 2, .., K)^n
  ## x_j: n-dim vector containing variable x_j to split on
  ## ds: terminal node to split
  ## j: id of variable to split on (integer in 1:p)
  ## threshold: value to split at
  ## p_Adj: an integer representing the number of available variable to split on
  ## n_adj: an integer representing the number of available unique values to split at
  ## alpha, beta: hyperparameters for the probability of eta being a terminal node
  ## sigma2: prior variance of parameter theta
  
  ###########################
  ## create T*
  
  D = T$leafCount
  # ds = sample(1:D, 1)
  
  tmp_parent = T$leaves[[ds]]
  depth = tmp_parent$level - 1
  K = length(tmp_parent$res)
  
  tmp_node_L = paste("d", T$totalCount+1, sep = "")
  tmp_node_R = paste("d", T$totalCount+2, sep = "")
  
  # tmp_parent$Climb()$AddChild(tmp_node_L)
  # tmp_parent$Climb()$AddChild(tmp_node_R)
  tmp_parent$AddChild(tmp_node_L)
  tmp_parent$AddChild(tmp_node_R)
  
  tmp_parent$children[[1]]$j = c(tmp_parent$j, j)
  tmp_parent$children[[1]]$sign = c(tmp_parent$sign, "L")
  tmp_parent$children[[1]]$threshold = c(tmp_parent$threshold, threshold)
  
  tmp_parent$children[[2]]$j = c(tmp_parent$j, j)
  tmp_parent$children[[2]]$sign = c(tmp_parent$sign, "R")
  tmp_parent$children[[2]]$threshold = c(tmp_parent$threshold, threshold)
  
  tmp_ix = tmp_parent$ix
  tmp_ix_L = intersect(which(x_j <= threshold), tmp_ix)
  tmp_ix_R = setdiff(tmp_ix, tmp_ix_L)
  
  tmp_parent$children[[1]]$ix = tmp_ix_L
  s_L = get_s(y[tmp_ix_L], K)
  # s_L = table(y[tmp_ix_L])
  tmp_parent$children[[1]]$res = s_L
  
  tmp_parent$children[[2]]$ix = tmp_ix_R
  s_R = get_s(y[tmp_ix_R], K)
  # s_R = table(y[tmp_ix_R])
  tmp_parent$children[[2]]$res = s_R
  
  Di_star_new = rep(0, length(y))
  for (d in 1:T$leafCount) {
    for (i in T$leaves[[d]]$ix) {
      Di_star_new[i] = d
    }
  }
  
  return(list(T, Di_star_new, depth))
}
