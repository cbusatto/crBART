pred_crBART = function(RES, Xtest, Bt) {
  
  n = nrow(Xtest)
  L = length(RES[[1]]$T)
  K = length(RES[[1]]$T[[1]]$T_l$res)
  
  # posterior mean probabilities
  pred_mat = matrix(0, n, K)
  
  for (b in 1:Bt) {
    
    # latent cumulative logits mu_i,k
    mu_mat = matrix(RES[[b]]$alpha, n, K-1, byrow = TRUE)
    
    for (l in 1:L) {
      # print("---------------------------------")
      
      Dl = RES[[b]]$T[[l]]$T_l$leafCount
      
      if (Dl == 1) {
        v = RES[[b]]$T[[l]]$mu
        mu_mat = mu_mat + matrix(v, n, K-1, byrow = TRUE)
        next
      }
      
      mu_l = RES[[b]]$T[[l]]$mu
      leaf_id = integer(n)
      for (d in 1:Dl) {
        tmp = RES[[b]]$T[[l]]$T_l$leaves[[d]]
        j  = tmp$j
        ts = tmp$threshold
        st = tmp$sign
        
        idx = 1:n
        if (d > 1) {
          idx = which(leaf_id == 0)
        }
        
        mask = rep(TRUE, length(idx))
        for (k in seq_along(j)) {
          if (st[k] == "L") {
            mask = mask & (Xtest[idx, j[k]] <= ts[k])
          } else {
            mask = mask & (Xtest[idx, j[k]] > ts[k])
          }
        }
        
        assigned = idx[mask]
        if (length(assigned) != 0) {
          leaf_id[assigned] = d
          if (length(assigned) > 0) {
            v = mu_l[((d-1)*(K-1)+1):(d*(K-1))]
            mu_mat[assigned, ] = mu_mat[assigned, ] + matrix(v, length(assigned), K-1, byrow=TRUE)
          }
        }
      }
    }
    
    # ---- Convert cumulative logits to probabilities ----
    
    log_prod_term = rep(0, n)
    log_exp_prob = matrix(-Inf, n, K)
    for (k in K:2) {
      prob_k = plogis(mu_mat[, k-1])
      log_exp_prob[, k] = log_prod_term + log(prob_k)
      log_prod_term = log_prod_term + log1p(-prob_k)
    }
    log_exp_prob[, 1] = log_prod_term
    
    maxlog = apply(log_exp_prob, 1, max)
    exp_prob = exp(log_exp_prob - matrix(maxlog, n, K))
    prob_mat = exp_prob / rowSums(exp_prob)
    pred_mat = pred_mat + prob_mat
  }
  
  pred_mat = pred_mat / Bt
  
  return(pred_mat)
}

