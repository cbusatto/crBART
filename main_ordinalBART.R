library(data.tree)
library(nnet)
library(pROC)

source("C:/Users/HP/Dropbox/Assegno_unifi/sample_data.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_nsin.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_s.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_res.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_prior_mean.R")
# source("C:/Users/HP/Dropbox/Assegno_unifi/get_prior_mean_prob.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/GROW_step0.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/GROW_step.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/PRUNE_step.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/CHANGE_step.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/ordinal_BART.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/predict_ordinalBART.R")


## Parameters
n = 100              # Number of observations
p = 200              # Number of predictors
K = 4                # Number of categories (classes)
p0 = c(5, 5, 7)
ntest = 100

## set hyperparameters
B = 1000
alpha = 0.95
beta = 2
L = 30
w0 = 0.5
s0 = 0.3
prob_move = rep(1/3, 3)

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## my data 

## beta0
set.seed(12)
K = 4                # Number of categories (classes)
beta0 = matrix(0, nrow = p, ncol = K-1)
for (k in 1:(K-1)) {
  u = rbinom(p0[k], 1, 0.4)
  beta0[1:p0[k], k] = (-1)^(u) * (0.75 * log(n) / sqrt(n) + abs(rnorm(p0[k], 0, 1)))
  # beta0[1:p0[k], k] = (-1)^(u) * abs(rnorm(p0[k], 0, 1))
}

## sample data
datai = sim_data_multinom(beta0, n, p, p0, K)
y = as.numeric(datai$y)
X = datai$X

datai = sim_data_multinom(beta0, ntest, p, p0, K)
ytest = as.numeric(datai$y)
Xtest = datai$X

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set.seed(10)
RES = ordinal_BART(y = as.numeric(y), X = X, 
                   L = L, w0 = w0, s0 = s0, 
                   mcmc_iter = B/2, mcmc_bn = B/2, 
                   alpha = alpha, beta = beta, 
                   prob_move = prob_move)

b = sample(1:round(B/2), 1)
l = sample(1:L, 1)

RES[[b]]$T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf)
round(RES[[b]]$T[[l]]$mu, 3)
round(RES[[b]]$alpha, 3)
round(RES[[b]]$mu, 3)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## PREDICTION

yprob = matrix(0, ntest, K)
tmp = matrix(0, B/2, ntest)
for (i in 1:ntest) {
  print("::::::::: iter:")
  print(i)
  g = pred_ordBART(RES, Xtest[i, ], B/2, y)
  
  yprob[i, ] = g[[1]]
  tmp[, i] = g[[2]]
}

ypred = apply(yprob, 1, function(x) which.max(x))
cbind(round(yprob, 3), ytest)
round(apply(round(yprob, 3), 2, mean), 3)
table(apply(yprob, 1, function(x) which.max(x)), ytest)
