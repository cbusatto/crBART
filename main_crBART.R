# ========================================================================
# LIBRARIES

library(data.tree)
library(nnet)
library(pROC)

# ========================================================================
# SOURCES

setwd("C:/Users/HP/Dropbox/Assegno_unifi")

source("sample_data.R")
source("get_nsin.R")
source("get_s.R")
source("get_res.R")
source("get_prior_mean.R")
source("GROW_step0.R")
source("GROW_step.R")
source("PRUNE_step.R")
source("CHANGE_step.R")
source("crBART.R")
source("predict_crBART.R")
source("var_importance_crBART.R")

# ========================================================================
# PARAMETERS

n     = 50             
p     = 20            
K     = 4                # Number of categories (classes)
p0    = c(5, 5, 7)
ntest = 100

# set hyperparameters
B         = 1000
alpha     = 0.95
beta      = 2
L         = 5
w0        = 0.5
s0        = 0.3
prob_move = rep(1/3, 3)

# ========================================================================
# DATA

set.seed(12)
beta0 = matrix(0, nrow = p, ncol = K-1)
for (k in 1:(K-1)) {
  u = rbinom(p0[k], 1, 0.4)
  beta0[1:p0[k], k] = (-1)^(u) * (0.75 * log(n) / sqrt(n) + abs(rnorm(p0[k], 0, 1)))
}

datai = sim_data_multinom(beta0, n, p, K)
y = as.numeric(datai$y)
X = datai$X

datai = sim_data_multinom(beta0, ntest, p, K)
ytest = as.numeric(datai$y)
Xtest = datai$X

# ========================================================================
# RUN crBART

set.seed(10)
res_crbart = crBART(y = as.numeric(y), X = X, 
                    L = L, w0 = w0, s0 = s0, 
                    mcmc_iter = B/2, mcmc_bn = B/2, 
                    alpha = alpha, beta = beta, 
                    prob_move = prob_move)

# variable importance
var_importance_crBART(res_crbart, B/2, p)

# prediction
yprob_crbart = pred_crBART(res_crbart, Xtest, B/2)
ypred_crbart = apply(yprob_crbart, 1, function(x) which.max(x))
