# get_prior_mean_prob = function(probs, sigma2 = 1, L = 1, eps = 1e-10) {
#   # method = match.arg(method)
#   # Ensure probs is a numeric vector and K >= 2
#   if (!is.numeric(probs) || length(probs) < 2) stop("probs must be numeric of length >= 2")
#   # normalize if needed
#   probs = probs / sum(probs)
#   K = length(probs)
#   scale = sqrt(1 + (sigma2 * pi^2) / 8)  # same scale you used
#   
#   # compute sequential alphas using the same factorization as your original code:
#   # for k = 1..K-1: p_k = alpha_k * prod_{j=k+1}^{K-1} (1 - alpha_j)
#   alpha = numeric(K-1)
#   # start from k = K-1 (empty product -> alpha_{K-1} = p_{K-1})
#   alpha[K-1] = pmin(pmax(probs[K-1], eps), 1-eps)
#   
#   if (K-2 >= 1) {
#     for (k in (K-2):1) {
#       prod_term = prod(1 - alpha[(k+1):(K-1)])
#       alpha[k] = probs[k] / prod_term
#       alpha[k] = pmin(pmax(alpha[k], eps), 1 - eps)
#     }
#   }
#   
#   # optional sanity check: implied p_K = prod_{j=1}^{K-1} (1 - alpha_j)
#   implied_pK = prod(1 - alpha)
#   # If you want, you can check closeness:
#   # if (abs(implied_pK - probs[K]) > 1e-6) warning("Implied p_K differs from supplied probs[K]")
#   
#   mu_means = log(alpha / (1 - alpha)) 
#   
#   # scale across L trees as in your original function
#   mu_means = mu_means / L
#   return(mu_means)
# }


# # target marginals for K = 4:
# probs = c(0.21, 0.22, 0.13, 0.44)
# a = get_prior_mean(probs, sigma2 = 1, L = 1, method = "logit")
# exp(a[3]) / (1 + exp(a[3]))
# exp(a[2]) / (1 + exp(a[2])) / (1 + exp(a[3]))
# exp(a[1]) / (1 + exp(a[1])) / (1 + exp(a[2])) / (1 + exp(a[3]))
# 
# 
# 
# # if you prefer logistic-scale mapping (mu such that logistic(mu) = alpha)
# get_prior_mean_from_probs(probs, sigma2 = 1, L = 1, method = "logit")


get_prior_mean_prob = function(probs) {
  K = length(probs)
  eps = 1e-10
  
  alpha = numeric(K-1)
  alpha[K-1] = pmin(pmax(probs[K-1], eps), 1-eps)
  for (k in (K-2):1) {
    prod_term = prod(1 - alpha[(k+1):(K-1)])
    alpha[k] = probs[k] / prod_term
    alpha[k] = pmin(pmax(alpha[k], eps), 1-eps)
  }
  a = log(alpha / (1 - alpha))   # logit mapping
  
  return(a)
}

# probs = table(y) / n
# probs = c(probs[-1], probs[1])
# probs2 = probs
# probs2[length(probs)] = probs[1]
# probs2[1] = probs[length(probs)]
# a = get_prior_mean_prob(probs)
# 
# exp(a[3]) / (1 + exp(a[3]))
# exp(a[2]) / (1 + exp(a[2])) / (1 + exp(a[3]))
# exp(a[1]) / (1 + exp(a[1])) / (1 + exp(a[2])) / (1 + exp(a[3]))
# 
# alpha
# a
