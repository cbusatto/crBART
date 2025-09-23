sim_data_multinom = function(beta0, n, p, p0, K) {
  
  if (length(p0) != (K-1)) return(0)
  
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
  
  return(list(y = y, X = X, beta0 = beta0, prob = probs))
}


sim_data_multinom_test = function(n, p, K, beta0) {
  
  ## Simulate predictors (independent variables)
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  # X = cbind(1, X)
  colnames(X) = paste0("X", 1:p)

  ## Compute linear predictor for each class
  linear_predictor = X %*% beta0
  
  ## Convert to probabilities using softmax function
  exp_lp = exp(cbind(0, linear_predictor))  # Add baseline class (logit = 0)
  probs = exp_lp / rowSums(exp_lp)
  
  ## Simulate response variable based on probabilities
  y = apply(probs, 1, function(prob) sample(1:K, size = 1, prob = prob))
  
  return(list(ytest = y, Xtest = X))
}

sim_data_ordinal = function(n, p, p0, K) {
  
  ## Simulate predictors (independent variables)
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  # X = cbind(1, X)
  colnames(X) = paste0("X", 1:p)
  
  ## Simulate true coefficients (including an intercept for each class)
  # beta = matrix(rnorm((p+1) * (K-1)), nrow = p+1, ncol = K-1)
  beta0 = rep(0, p)
  u = rbinom(p0[1], 1, 0.4)
  beta0[1:p0[1]] = (-1)^(u) * (0.75 * log(n) / sqrt(n) + abs(rnorm(p0[1], 0, sqrt(0.75))))
  
  # Compute the latent variable
  latent = X %*% beta0 + rnorm(n)  # Linear predictor + noise
  
  # Define cutpoints using quantiles
  cuts = quantile(latent, probs = seq(0, 1, length.out = K+1))
  
  # Assign ordinal levels based on cutpoints
  ordinal_response = cut(latent, breaks = cuts, labels = FALSE, include.lowest = TRUE)
  
  # Convert to an ordered factor
  ordinal_response = ordered(ordinal_response, levels = 1:K)
  
  # Create a data frame
  # data = data.frame(X, y = ordinal_response)
  
  return(list(y = ordinal_response, X = X, beta0 = beta0))
}


sim_data_ordinal_test = function(n, p, K, beta0) {
  
  ## Simulate predictors (independent variables)
  X = matrix(rnorm(n*p), nrow = n, ncol = p)
  # X = cbind(1, X)
  colnames(X) = paste0("X", 1:p)
  
  # Compute the latent variable
  latent = X %*% beta0 + rnorm(n)  # Linear predictor + noise
  
  # Define cutpoints using quantiles
  cuts = quantile(latent, probs = seq(0, 1, length.out = K+1))
  
  # Assign ordinal levels based on cutpoints
  ordinal_response = cut(latent, breaks = cuts, labels = FALSE, include.lowest = TRUE)
  
  # Convert to an ordered factor
  ordinal_response = ordered(ordinal_response, levels = 1:K)
  
  # Create a data frame
  # data = data.frame(X, y = ordinal_response)
  
  return(list(ytest = ordinal_response, Xtest = X))
}


sample_data_multinom = function(n, K, p, p0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (p0 > p) stop("p0 must be ≤ p")
  
  # Generate β_k vectors (one for each class k)
  beta_list = vector("list", K)
  for (k in 1:K) {
    beta_k = rep(0, p)
    nonzero_idx = sample(1:p, p0, replace = FALSE)
    beta_k[nonzero_idx] = rnorm(p0)
    beta_list[[k]] = beta_k
  }
  
  # Generate design matrix X: n x p
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Compute η_{ik} = x_i^T β_k and π_{ik} via softmax
  eta_mat = matrix(0, nrow = n, ncol = K)
  for (k in 1:K) {
    eta_mat[, k] = X %*% beta_list[[k]]
  }
  
  # Softmax function (row-wise)
  softmax = function(eta_row) {
    exp_eta = exp(eta_row - max(eta_row))  # stabilize
    exp_eta / sum(exp_eta)
  }
  
  prob_mat = t(apply(eta_mat, 1, softmax))  # n x K
  
  # Sample multinomial outcomes y_i in 1..K
  y = apply(prob_mat, 1, function(p) sample(1:K, size = 1, prob = p))
  
  return(list(
    y = y,
    X = X,
    beta = beta_list,
    probs = prob_mat
  ))
}
