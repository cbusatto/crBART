pred_ordBART = function(RES, xtest, Bt) {
  L = length(RES[[1]]$T)
  K = length(RES[[1]]$T[[1]]$T_l$res)
  
  ypred = rep(0, length(RES))
  pred = rep(0, K)
  for (b in 1:Bt) {

    mu = RES[[b]]$alpha
    
    for (l in 1:L) {
      Dl = RES[[b]]$T[[l]]$T_l$leafCount
      if (Dl == 1) {
        v = RES[[b]]$T[[l]]$mu
        mu = mu + v
      } else {
        
        ds = 0
        for (d in 1:Dl) {
          tmp = RES[[b]]$T[[l]]$T_l$leaves[[d]]
          
          j = tmp$j
          ts = tmp$threshold
          stmp = tmp$sign
          
          key = 0
          js = 1
          while ((key == 0) & (js <= length(j))) {
            if (stmp[js] == "L") {
              key = ifelse(xtest[j[js]] <= ts[js], 0, 1)
            } else {
              key = ifelse(xtest[j[js]] > ts[js], 0, 1)
            }
            
            js = js + 1
          }
          
          if (key == 0) {
            ds = d
          }
        }
        
        if (ds == 0) print("Warning")
        tmp = RES[[b]]$T[[l]]$T_l$leaves[[ds]]
          
        mu_l = RES[[b]]$T[[l]]$mu
        v = mu_l[((ds-1)*(K-1)+1):(ds*(K-1))]
        mu = mu + v
      }
    }

    # prod_term = 1
    # exp_prob = rep(0, K)
    # for (k in K:2) {
    #   prob_k = plogis(mu[k-1])
    #   exp_prob[k] = prod_term * prob_k
    #   prod_term = prod_term * (1-prob_k)
    # }
    # exp_prob[1] = prod_term
    log_prod_term = 0
    log_exp_prob = rep(-Inf, K)
    
    for (k in K:2) {
      prob_k = plogis(mu[k-1])              # P(y=k | not in 1,...,k-1)
      log_exp_prob[k] = log_prod_term + log(prob_k)
      log_prod_term = log_prod_term + log1p(-prob_k)  # log(1 - prob_k)
    }
    log_exp_prob[1] = log_prod_term
    
    # exponentiate & normalize
    maxlog = max(log_exp_prob)              # for stability
    exp_prob = exp(log_exp_prob - maxlog)
    exp_prob = exp_prob / sum(exp_prob)

    pred = pred + exp_prob
    ypred[b] = which.max(exp_prob)
  }

  return(list(pred / Bt, ypred))
}


pred_ordBART_mat = function(RES, Xtest, Bt) {
  
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









pred_ordBART_mat2 = function(RES, Xtest, Bt) {
  
  n = nrow(Xtest)
  L = length(RES[[1]]$T)
  K = length(RES[[1]]$T[[1]]$T_l$res)
  
  # posterior mean probabilities
  pred_mat = matrix(0, n, K)
  
  for (b in 1:Bt) {
    
    # latent cumulative logits mu_i,k
    mu_mat = matrix(RES[[b]]$alpha, n, K-1, byrow = TRUE)
    
    for (l in 1:L) {
      Dl = RES[[b]]$T[[l]]$T_l$leafCount
      
      # stump (single leaf)
      if (Dl == 1) {
        v = RES[[b]]$T[[l]]$mu
        mu_mat = mu_mat + matrix(v, n, K-1, byrow = TRUE)
        next
      }
      
      # multiple leaves
      for (i in 1:n) {
        ds = 0
        
        for (d in 1:Dl) {
          tmp = RES[[b]]$T[[l]]$T_l$leaves[[d]]
          
          j  = tmp$j
          ts = tmp$threshold
          stmp = tmp$sign
          
          key = 0
          js = 1
          while ((key == 0) & (js <= length(j))) {
            if (stmp[js] == "L") {
              key = ifelse(Xtest[i, j[js]] <= ts[js], 0, 1)
            } else {
              key = ifelse(Xtest[i, j[js]] > ts[js], 0, 1)
            }
            js = js + 1
          }
          
          if (key == 0) ds = d
        }
        
        if (ds == 0) stop("Tree traversal error")
        
        mu_l = RES[[b]]$T[[l]]$mu
        v = mu_l[((ds-1)*(K-1)+1):(ds*(K-1))]
        mu_mat[i, ] = mu_mat[i, ] + v
      }
    }
    
    # ---- Convert cumulative logits to probabilities ----
    
    prob_mat = matrix(0, n, K)

    for (i in 1:n) {
      log_prod_term = 0
      log_exp_prob = rep(-Inf, K)

      for (k in K:2) {
        prob_k = plogis(mu_mat[i, k-1])
        log_exp_prob[k] = log_prod_term + log(prob_k)
        log_prod_term = log_prod_term + log1p(-prob_k)
      }
      log_exp_prob[1] = log_prod_term

      maxlog = max(log_exp_prob)
      exp_prob = exp(log_exp_prob - maxlog)
      prob_mat[i, ] = exp_prob / sum(exp_prob)
    }

    pred_mat = pred_mat + prob_mat
  }
  
  pred_mat = pred_mat / Bt
  
  return(pred_mat)
}
