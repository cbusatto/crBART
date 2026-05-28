sim_data_multinom = function(beta0, n, p, K) {
  
  ## Simulate predictors (independent variables)
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  # X = cbind(1, X)
  colnames(X) = paste0("X", 1:p)
  
  ## Simulate true coefficients (including an intercept for each class)
  # beta = matrix(rnorm((p+1) * (K-1)), nrow = p+1, ncol = K-1)
  # beta0 = matrix(0, nrow = p, ncol = K-1)
  # for (k in 1:(K-1)) {
  #   u = rbinom(p0[k], 1, 0.4)
  #   beta0[1:p0[k], k] = (-1)^(u) * (0.75 * log(n) / sqrt(n) + abs(rnorm(p0[k], 0, 1)))
  #   # beta0[1:p0[k], k] = (-1)^(u) * abs(rnorm(p0[k], 0, 1))
  # }
  
  ## Compute linear predictor for each class
  linear_predictor = X %*% beta0
  
  ## Convert to probabilities using softmax function
  exp_lp = exp(cbind(0, linear_predictor))  # Add baseline class (logit = 0)
  probs = exp_lp / rowSums(exp_lp)
  
  ## Simulate response variable based on probabilities
  y = apply(probs, 1, function(prob) sample(1:K, size = 1, prob = prob))
  
  return(list(y = y, X = X))
}


sim_data_ordinal = function(beta0, n, p, K) {
  
  ## Simulate predictors (independent variables)
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  colnames(X) = paste0("X", 1:p)
  
  ## Compute linear predictor for each class
  linear_predictor = X %*% beta0
  linear_predictor[, K-1] = linear_predictor[, K-1]-runif(n, 1, 2.5)
  linear_predictor[, K-2] = linear_predictor[, K-2]-runif(n, 0, 0.5)
  
  ## Sample response variable based on ordered probabilities
  y = numeric(n)
  probs = matrix(0, n, K)
  for (i in 1:n) {
    mu = linear_predictor[i, ]
    
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
    
    y[i] = sample(1:K, size = 1, prob = exp_prob)
    probs[i, ] = exp_prob
  }
  
  return(list(y = y, X = X, probs = probs))
}

sim_data_ordinal_unbalanced = function(beta0, n, p, K) {
    
    ## Simulate predictors (independent variables)
    X = matrix(rnorm(n*p), nrow = n, ncol = p)
    colnames(X) = paste0("X", 1:p)
    
    ## Compute linear predictor for each class
    linear_predictor = X %*% beta0
    linear_predictor[, K-1] = linear_predictor[, K-1]-runif(n, 3, 5)
    linear_predictor[, K-2] = linear_predictor[, K-2]-runif(n, 1, 3)
    
    ## Sample response variable based on ordered probabilities
    y = numeric(n)
    probs = matrix(0, n, K)
    for (i in 1:n) {
      mu = linear_predictor[i, ]
      
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
      
      y[i] = sample(1:K, size = 1, prob = exp_prob)
      probs[i, ] = exp_prob
    }
    
    return(list(y = y, X = X, probs = probs))
  }
  


