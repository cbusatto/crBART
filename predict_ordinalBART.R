pred_ordBART = function(RES, xtest, Bt, y) {
  # B = length(RES)
  L = length(RES[[1]]$T)
  K = length(RES[[1]]$T[[1]]$T_l$res)
  
  ypred = rep(0, length(RES))
  pred = rep(0, K)
  
  for (b in 1:Bt) {
    print(":::::::: iter ")
    print(b)

    mu = RES[[b]]$alpha
    # mu = alpha0_mean
    
    for (l in 1:L) {
      # print(" :: l: ")
      # print(l)
      # print(xtest)

      # print(RES[[b]]$T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf))
      # print(RES[[b]]$T[[l]]$mu)
      
      Dl = RES[[b]]$T[[l]]$T_l$leafCount
      if (Dl == 1) {
        
        v = RES[[b]]$T[[l]]$mu
        mu = mu + v
        # for (k in 1:(K-1)) {
        #   mu[k] = mu[k] + v[k]
        # }
        
      } else {
        
        ds = 0
        for (d in 1:Dl) {
          tmp = RES[[b]]$T[[l]]$T_l$leaves[[d]]
          
          j = tmp$j
          ts = tmp$threshold
          stmp = tmp$sign

          # print(tmp$j)
          # print(tmp$sign)
          # print(tmp$threshold)
          
          key = 0
          js = 1
          while ((key == 0) & (js <= length(j))) {
            # for (js in 1:length(j)) {
            
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
        # print(ds)
        tmp = RES[[b]]$T[[l]]$T_l$leaves[[ds]]
          
        mu_l = RES[[b]]$T[[l]]$mu
        v = mu_l[((ds-1)*(K-1)+1):(ds*(K-1))]
        mu = mu + v
        # print("mu and ds")
        # print(v)
        # print(tmp$res)
        # print(tmp$ix)
        # print(y[tmp$ix])
        # print(table(y[tmp$ix]))
        # 
        # tmp = v
        # print(exp(tmp[1]) / (1 + exp(tmp[1])))
        # print(exp(sum(tmp[1:2])) / (1 + exp(sum(tmp[1:2]))) - exp(tmp[1]) / (1 + exp(tmp[1])))
        # print(exp(sum(tmp[1:3])) / (1 + exp(sum(tmp[1:3]))) - exp(sum(tmp[1:2])) / (1 + exp(sum(tmp[1:2]))))
        # print(1 / (1 + exp(sum(tmp[1:3]))))
        
        # print(ds)
        # 
        # print(tmp$j)
        # print(tmp$sign)
        # print(tmp$threshold)
          
        # for (k in 1:(K-1)) {
        #   # print(((ds-1)*(K-1)+1):(ds*(K-1)))
        #   mu[k] = mu[k] + v[k]
        # }
      }
    }

    # print("mu")
    # print(mu)
    
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
    
    # 
    # # print("mu")
    # # print(mu)
    # # mu = mu / L
    # exp_prob = rep(0, K)
    # # exp_prob[1] = exp(mu[1]) / (1 + exp(mu[1]))
    # exp_prob[1] = plogis(mu[1])
    # # exp_prob[K] = 1 / (1 + exp(mu[K-1]))
    # exp_prob[K] = 1 - plogis(mu[K-1])
    # for (k in 2:(K-1)) {
    #   # exp_prob[k] = exp(mu[k]) / (1 + exp(mu[k])) - exp(mu[k-1]) / (1 + exp(mu[k-1]))
    #   exp_prob[k] = plogis(mu[k]) - plogis(mu[k-1])
    # }
    # # for (k in 2:(K-1)) {
    # #   exp_prob[k] = exp(mu[k]) / (1 + exp(mu[k]))
    # # }
    # print(exp_prob)
    

    # print(exp_prob)
    pred = pred + exp_prob
    ypred[b] = which.max(exp_prob)
  }

  return(list(pred / Bt, ypred))
}

