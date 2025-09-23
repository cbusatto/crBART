CHANGE_step = function(T, y, x_j, ds, j, threshold) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## function to CHANGE a singly internal node of a tree T
  
  ## INPUT:
  ## T: current tree
  ## y: response variable of the type y = (1, 2, .., K)^n
  ## x_j: n-dim vector containing variable x_j to split on
  ## ds: terminal node whose parent is pruned
  ## j: id of variable to split on (integer in 1:p)
  ## threshold: value to split at
  ## p_Adj: an integer representing the number of available variable to split on
  ## n_adj: an integer representing the number of available unique values to split at
  ## alpha, beta: hyperparameters for the probability of eta being a terminal node
  ## epsilon: tuning parameter representing the variance of proposal q \sim N(m_omega, epsilon)
  ## sigma2: prior variance of parameter theta
  
  D = T$leafCount
  
  tmp_parent = T$leaves[[ds]]$parent
  K = length(tmp_parent$res)
  
  depth = length(tmp_parent$children[[1]]$j)
  tmp_parent$children[[1]]$j[depth] = j
  # tmp_parent$children[[1]]$sign = c(tmp_parent$sign, "L")
  tmp_parent$children[[1]]$threshold[depth] = threshold
  
  tmp_parent$children[[2]]$j[depth] = j
  # tmp_parent$children[[2]]$sign = c(tmp_parent$sign, "R")
  tmp_parent$children[[2]]$threshold[depth] = threshold
  
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

  return(list(T, Di_star_new))
}
