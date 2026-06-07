crBART = function(y, X, L, w0 = 0.5, s0 = 0.3, mcmc_iter = 1000, mcmc_bn = 500, alpha = 0.95, beta = 2, prob_move = rep(1/3, 3)) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Ordinal BART
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
  ## INPUT:
  ## y: response variable of the type y = (1, 2, .., K)^n
  ## X: data matrix n x p
  ## L: number of trees T_l
  ## w0: positive value related to prior variance sigma2_l (w0 = 0.5 default)
  ## s0: positive value related to prior variance sigma2_alpha (s0 = 0.3 default)
  ## mcmc_iter: number of post-burnin iterations
  ## mcmc_bn: number of burnin iterations
  ## alpha, beta: hyperparameters for the probability of being a terminal node
  ## prob_move: probability of GROW, PRUNE and CHANGE steps
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## INITIALIZATION
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  n = length(y)
  p = ncol(X)
  s = table(y)
  K = length(s)
  
  D = rep(1, L)
  depth = rep(0, L)
  Di_star = matrix(1, n, L)
  
  ## ::::: compute latent factor z_tilde
  
  ms = K - pmax(1, y-1) 
  offset = c(0, cumsum(ms))[1:n] + 1  
  
  z0 = c()
  for (i in 1:n) {
    if (y[i] == 1) {
      z0 = c(z0, rep(-0.5, K-1))
    } else if (y[i] == K) {
      z0 = c(z0, 0.5)
    } else {
      z0 = c(z0, 0.5, rep(-0.5, K-y[i]))
    }
  }
  ns = length(z0)
  
  ## ::::: initialize T1, .., T_L, sigma2
  
  sigma2 = numeric(L)
  T = list()
  
  for (l in 1:L) {
    sigma2[l] = (3 / (w0 * sqrt(l)))^2
    # sigma2[l] = (3 / (w0 * sqrt(L)))^2
    
    T[[l]] = list()
    T[[l]]$rules = "T"
    T[[l]]$T_l = Node$new("T")
    T[[l]]$T_l$ix = 1:n
    T[[l]]$T_l$res = s
  }

  ## ::::: initialiaze alpha and mu
  
  # prior mean
  # alpha0_mean = get_prior_mean(K, sigma2 = sigma2[1], L = L)
  alpha0_mean = get_prior_mean(K, sigma2 = s0, L = 1)
  alpha0 = alpha0_mean
  
  mu = rep(0, L*(K-1))
  for (l in 1:L)  T[[l]]$mu = rep(0, K-1)
  
  mu_x = matrix(0, n, K-1)
  for (i in 1:n) {
    for (k in max(2, y[i]):K)  mu_x[i, k-1] = alpha0_mean[k-1]
  }
  
  T_out = list()
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## MAIN CYCLE
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  cat(":::::::::::: START", "\n")
  
  for (b in 1:(mcmc_bn+mcmc_iter)) {
    
    if (b %% 500 == 0) {
      cat(":::::::::::: iteration:", b, "\n")
    }

    ## ::::: update latent factors omega and the intercept
    
    c_tmp = 1
    w = numeric(ns)
    eta_nol = numeric(ns)
    omega = matrix(0, n, K-1)
    for(i in 1:n) {
      for (k in max(2, y[i]):K) {
        tmp = mu_x[i, k-1]
        w[c_tmp] = pgdraw::pgdraw(1, tmp)
        omega[i, k-1] = w[c_tmp]
        eta_nol[c_tmp] = tmp - alpha0[k-1]
        c_tmp = c_tmp + 1
      }
    }
    r = z0 / w - eta_nol
    
    var_post = c()
    mean_post = c()
    for (k in 1:(K-1)) {
      idx = which(y <= (k+1))  
      ind_pos = c(offset[idx] + k - pmax(1, y[idx]-1))
      w_sum = sum(omega[idx, k])
      r_sum = sum(omega[idx, k] * r[ind_pos]) 
      
      # posterior var and mean
      tmp_var = 1 / (w_sum + 1 / s0)
      var_post = c(var_post, tmp_var)   
      mean_post = c(mean_post, tmp_var * (r_sum + alpha0_mean[k] / s0)) 
    }
    
    tmp_alpha0 = alpha0
    alpha0 = rnorm(K-1, mean_post, sqrt(var_post))
    
    for (i in 1:n) {
      for (k in max(2, y[i]):K) {
        mu_x[i, k-1] = mu_x[i, k-1] + alpha0[k-1] - tmp_alpha0[k-1]
      } 
    }
    
    ## ::::::::::::::::::::: 
    ## update T and M
    
    for (l in 1:L) {
      # print("l")
      # print(l)
      # print(T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf))
      
      if (depth[l] == 0) {
        
        ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
        ## GROW STEP 0: there is only root
        ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
       
        check = 1
        while (check == 1) {
          j = sample(1:p, 1)
          id_n_adj = sample(1:n, 1)
          threshold = X[id_n_adj, j]
          check = ifelse(((threshold == max(X[, j])) || (threshold == min(X[, j]))) == TRUE, 1, 0)
        }
        tmp_rule = paste("x", j, ".", id_n_adj, sep = "")
        
        T_star = data.tree::Clone(T[[l]]$T_l)
        TMP = GROW_step0(T_star, y, X[, j], j, threshold)
        T_star = TMP[[1]]
        
        D_star = D
        D_star[l] = D[l] + 1
        Di_star_new = Di_star
        Di_star_new[, l] = TMP[[2]]
      
      } else {
        
        check = apply(T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf), 2, sum)
        if (max(check) <= 5) {
          move = 3
        } else {
          move = sample(1:3, 1, prob = prob_move)
        }
        # move = sample(1:3, 1, prob = prob_move)
        
        if (move == 1) {
          
          ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
          ## GROW STEP
          ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
          
          check = 1
          while (check == 1) {
            ds = sample(1:D[l], 1)
            if (length(T[[l]]$T_l$leaves[[ds]]$ix) > 5) {
              check = 0
            }
          }
          
          check = 1
          while (check == 1) {
            j = sample(1:p, 1)
            ind_tmp = T[[l]]$T_l$leaves[[ds]]$ix
            id_n_adj = sample(ind_tmp, 1)
            threshold = X[id_n_adj, j]
            check = ifelse(((threshold == min(X[ind_tmp, j])) || (threshold == max(X[ind_tmp, j]))) == TRUE, 1, 0) 
          }
          
          tmp_rule = paste("x", j, ".", id_n_adj, sep = "")
          
          T_star = data.tree::Clone(T[[l]]$T_l)
          TMP = GROW_step(T_star, y, X[, j], ds, j, threshold)
          T_star = TMP[[1]]
          
          depth_eta = TMP[[3]]
          
          D_star = D
          D_star[l] = D[l] + 1
          Di_star_new = Di_star
          Di_star_new[, l] = TMP[[2]]
          
          # print(T_star$Get(function(node) node$res, filterFun = isLeaf))
          
        } else if (move == 2) {
          
          ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
          ## PRUNE STEP
          ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
          
          check = FALSE
          while (check == FALSE) {
            ds = sample(1:D[l], 1)
            tmp_parent = T[[l]]$T_l$leaves[[ds]]$parent
            check = length(tmp_parent$children) == 2 && all(sapply(tmp_parent$children, function(x) x$isLeaf))
          }
          
          depth_eta = tmp_parent$level - 1
          j = T[[l]]$T_l$leaves[[ds]]$j
          ts = T[[l]]$T_l$leaves[[ds]]$threshold
          if (length(j) > 1) {
            ts = ts[length(j)]
            j = j[length(j)]
          }
          id_n_adj = which(X[, j] == ts)
          tmp_rule = paste("x", j, ".", id_n_adj, sep = "")
          
          T_star = data.tree::Clone(T[[l]]$T_l)
          TMP = PRUNE_step(T_star, ds, n)
          T_star = TMP[[1]]
          
          D_star = D
          D_star[l] = D[l] - 1
          Di_star_new = Di_star
          Di_star_new[, l] = TMP[[2]]
          
          # print(T_star$Get(function(node) node$res, filterFun = isLeaf))
          
        } else {
          
          ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
          ## CHANGE STEP
          ## ::::::::::::::::::::::::::::::::::::::::::::::::::::
          
          ## remove children
          check = FALSE
          while (check == FALSE) {
            ds = sample(1:D[l], 1)
            tmp_parent = T[[l]]$T_l$leaves[[ds]]$parent
            check = (length(tmp_parent$children) == 2) && all(sapply(tmp_parent$children, function(x) x$isLeaf))
          }
          
          j = T[[l]]$T_l$leaves[[ds]]$j
          ts = T[[l]]$T_l$leaves[[ds]]$threshold
          if (length(j) > 1) {
            ts = ts[length(j)]
            j = j[length(j)]
          }
          id_n_adj = which(X[, j] == ts)
          tmp_rule = paste("x", j, ".", id_n_adj, sep = "")
          
          ## add new split
          check = 1
          while (check == 1) {
            j2 = sample(1:p, 1)
            ind_tmp = tmp_parent$ix
            id_n_adj = sample(ind_tmp, 1)
            threshold = X[id_n_adj, j2]
            check = ifelse(((threshold == min(X[ind_tmp, j2])) || (threshold == max(X[ind_tmp, j2]))) == TRUE, 1, 0) 
          }
          tmp_rule2 = paste("x", j2, ".", id_n_adj, sep = "")
          
          T_star = data.tree::Clone(T[[l]]$T_l)
          TMP = CHANGE_step(T_star, y, X[, j2], ds, j2, threshold)
          T_star = TMP[[1]]
          
          D_star = D
          Di_star_new = Di_star
          Di_star_new[, l] = TMP[[2]]
          
          # print(T_star$Get(function(node) node$res, filterFun = isLeaf))
        }
      }

      # print(T_star$Get(function(node) node$res, filterFun = isLeaf))
      
      ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ## compute quantities related to T_l and T*_l
      
      c_tmp = 1
      for(i in 1:n) {
        dl = Di_star[i, l]
        
        if (l > 1) {
          l_start = sum(D[1:(l-1)])
        } else {
          l_start = 0
        }
        l_start = l_start*(K-1) + (dl-1)*(K-1) + 1
        
        for (k in max(2, y[i]):K) {
          eta_nol[c_tmp] = mu_x[i, k-1] - mu[l_start+k-2]
          c_tmp = c_tmp + 1
        }
      }
      r = z0 / w - eta_nol
      
      var_post = c()
      mean_post = c()
      for (dl in 1:D[l]) {
        for (k in 1:(K-1)) {
          idx = which(Di_star[, l] == dl & y <= (k+1))  
          
          if (length(idx) > 0) {
            ind_pos = c(offset[idx] + k - pmax(1, y[idx]-1))
            w_sum = sum(omega[idx, k])
            r_sum = sum(omega[idx, k] * r[ind_pos]) 
          } else {
            w_sum = 0
            r_sum = 0
          }
          
          # posterior var and mean
          tmp_var = 1 / (w_sum + 1 / sigma2[l])
          var_post = c(var_post, tmp_var)   
          mean_post = c(mean_post, tmp_var * r_sum) 
        }
      }
      
      var_post_new = c()
      mean_post_new = c()
      for (dl in 1:D_star[l]) {
        for (k in 1:(K-1)) {
          idx_new = which(Di_star_new[, l] == dl & y <= (k+1))  
          if (length(idx_new) > 0) {
            ind_pos_new = c(offset[idx_new] + k - pmax(1, y[idx_new]-1))
            w_sum_new = sum(omega[idx_new, k])
            r_sum_new = sum(omega[idx_new, k] * r[ind_pos_new])
          } else {
            w_sum_new = 0
            r_sum_new = 0
          }
          
          # posterior var and mean
          tmp_var = 1 / (w_sum_new + 1 / sigma2[l])
          var_post_new = c(var_post_new, tmp_var)   
          mean_post_new = c(mean_post_new, tmp_var * r_sum_new) 
        }
      }
      
      ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ## compute alpha_MH
      
      ## ratios pi(T*) / pi(T) and q(T* -> T) / q(T -> T*)
      if (depth[l] == 0) {
        depth_eta = 0
        log_pqT_r = log(alpha) + 2*log(1-(alpha/2^beta)) - log(1-alpha)
      } else {
        if (move == 1) {
          D2_star = get_nsin(T_star)
          log_pqT_r = log(D[l]) - log(D2_star) + log(alpha) + 2*log(1-alpha/(2+depth_eta)^beta) - log((1+depth_eta)^beta-alpha)
        } else if (move == 2) {
          D2_star = get_nsin(T[[l]]$T_l)
          log_pqT_r = -log(D[l]-1) + log(D2_star) - log(alpha) - 2*log(1-alpha/(2+depth_eta)^beta) + log((1+depth_eta)^beta-alpha)
        } else {
          log_pqT_r = 0
        }
      }
      
      log_target = -0.5 * D[l] * (K-1) * log(2*pi) - 0.5 * sum(log(var_post)) - 0.5 * D[l] * (K-1) * log(sigma2[l]) + 0.5 * sum(mean_post^2 / var_post) 
      log_target_star = -0.5 * D_star[l] * (K-1) * log(2*pi) - 0.5 * sum(log(var_post_new)) - 0.5 * D_star[l] * (K-1) * log(sigma2[l]) + 0.5 * sum(mean_post_new^2 / var_post_new) 
      log_alpha_MH = log_pqT_r + log_target_star - log_target
      
      D_tmp = D[l]
      Di_star_old = Di_star
      check = log_alpha_MH > log(runif(1))
      if (check == TRUE) {
        T[[l]]$T_l = data.tree::Clone(T_star)
        
        D = D_star
        Di_star = Di_star_new
        mean_post = mean_post_new
        var_post = var_post_new
        
        if (depth[l] == 0 || move == 1) {
          T[[l]]$rules = c(T[[l]]$rules, tmp_rule)
          depth[l] = depth[l] + 1
        } else if (move == 2) {
          T[[l]]$rules = T[[l]]$rules[-which(T[[l]]$rules == tmp_rule)]
          depth[l] = depth[l] - 1
        } else {
          T[[l]]$rules = T[[l]]$rules[-which(T[[l]]$rules == tmp_rule)]
          T[[l]]$rules = c(T[[l]]$rules, tmp_rule2)
        }
      }
         
      ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ## sample mu
      
      tmp_mu = T[[l]]$mu
      tmp_mu_new = rnorm(D[l]*(K-1), mean_post, sqrt(var_post))
      T[[l]]$mu = tmp_mu_new
      
      for (i in 1:n) {
        dl = Di_star[i, l]
        tmp_dl = Di_star_old[i, l]
        ind_pos_new = (dl-1)*(K-1) + 1
        ind_pos_old = (tmp_dl-1)*(K-1) + 1
        for (k in max(2, y[i]):K)  {
          mu_x[i, k-1] = mu_x[i, k-1] + tmp_mu_new[ind_pos_new+k-2] - tmp_mu[ind_pos_old+k-2]
        }
      }
      
      
      if (l == 1) {
        mu = c(tmp_mu_new, mu[-c(1:(D_tmp*(K-1)))])
      } else if (l == L) {
        mu = c(mu[-c((length(mu)-D_tmp*(K-1)+1):length(mu))], tmp_mu_new)
      } else {
        l_start = sum(D[1:(l-1)]) * (K-1) + 1
        tmp = c(mu[1:(l_start-1)], tmp_mu_new)
        mu = c(tmp, mu[-c(1:(l_start+D_tmp*(K-1)-1))])
      }
    }
    
    if (b > mcmc_bn) {
      T_out[[b-mcmc_bn]] = list()
      T_out[[b-mcmc_bn]]$mu = mu
      T_out[[b-mcmc_bn]]$alpha = alpha0
      T_out[[b-mcmc_bn]]$T = T
    }
  }
  
  return(T_out)
}
