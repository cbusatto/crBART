var_importance_crBART = function(RES, Bt, p) {
  L = length(RES[[1]]$T)
  
  count_j = rep(0, p)
  count_tot = 0
  
  for (b in 1:Bt) {
    for (l in 1:L) {
      v = RES[[b]]$T[[l]]$rules[-1]
      ind_v = as.integer(sub("^x([0-9]+)\\..*", "\\1", v))
      
      count_j[ind_v] = count_j[ind_v] + 1
      count_tot = count_tot + length(ind_v)
    }
  }
  
  out = count_j / count_tot
  return (out)
}