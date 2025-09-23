get_prior_mean = function(K, sigma2 = 1, L = 1) {
  scale = sqrt(1 + (sigma2*pi^2)/8)
  
  mu_means = numeric(K-1)  
  alpha = numeric(K-1)     
  mu_means[K-1] = -log(K-1)
  alpha[K-1] = exp(mu_means[K-1]) / (1 + exp(mu_means[K-1]))  
  
  # Recursively compute alpha_k and mu_k backwards
  for (k in (K-2):1) {
    prod_term = prod(1 - alpha[(k+1):(K-1)])
    alpha[k] = 1 / (K*prod_term)
    alpha[k] = min(max(alpha[k], 1e-10), 1 - 1e-10)
    mu_means[k] = qnorm(alpha[k]) * scale
  }
  
  # Scale means for the number of trees L
  mu_means = mu_means / L
  # names(mu_means) = paste0("k=", 2:K)
  
  return(mu_means)
}