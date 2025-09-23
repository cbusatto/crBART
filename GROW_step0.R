GROW_step0 = function(T, y, x_j, j, threshold) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## function to GROW first node of a tree T
  
  ## INPUT:
  ## T: current tree, with only root
  ## y: response variable of the type y = (1, 2, .., K)^n
  ## x_j: n-dim vector containing variable x_j to split on
  ## j: id of variable to split on (integer in 1:p)
  ## threshold: value to split at
  
  K = length(T$res)
    
  ## create T*
  T$AddChild("d1")
  T$AddChild("d2")
  
  T$d1$j = c(j)
  T$d1$sign = "L"
  T$d1$threshold = threshold
  
  T$d2$j = c(j)
  T$d2$sign = "R"
  T$d2$threshold = threshold
  
  tmp_ix_L = which(x_j <= threshold)
  # tmp_ix_R = which(x_j > threshold)
  tmp_ix_R = setdiff(1:length(y), tmp_ix_L)
  
  T$d1$ix = tmp_ix_L
  s_L = get_s(y[tmp_ix_L], K)
  T$d1$res = s_L
  # T$d1$res = table(y[tmp_ix_L])
  
  T$d2$ix = tmp_ix_R
  s_R = get_s(y[tmp_ix_R], K)
  T$d2$res = s_R
  # T$d2$res = table(y[tmp_ix_R])
  
  Di_star_new = rep(0, length(y))
  for (d in 1:T$leafCount) {
    for (i in T$leaves[[d]]$ix) {
      Di_star_new[i] = d
    }
  }
  
  return(list(T, Di_star_new))
}


