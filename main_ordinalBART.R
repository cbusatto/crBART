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
source("C:/Users/HP/Dropbox/Assegno_unifi/var_importance_ordinalBART.R")


## Parameters
n = 100              # Number of observations
p = 20            # Number of predictors
K = 4                # Number of categories (classes)
p0 = c(5, 5, 7)
ntest = 100

## set hyperparameters
B = 1000
alpha = 0.95
beta = 2
L = 5
# v_narisetty = max(p^2.1/(100*n), log(n))
# w0 = sqrt(9/v_narisetty)
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
datai = sim_data_multinom(beta0, n, p, K)
y = as.numeric(datai$y)
X = datai$X

datai = sim_data_multinom(beta0, ntest, p, K)
ytest = as.numeric(datai$y)
Xtest = datai$X

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## data from https://rawcdn.githack.com/jaeyonggy/OPBART/main/Demo/demonstration.nb.html#opbart

# fx = function(x) -7 + 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]
# 
# gen_data = function(n, us) {
#   X = matrix(runif(n*5), nrow = n)
#   mu = fx(X)
#   Z = rnorm(n, mean = mu, sd = 1)
#   Y = ifelse(-Inf < Z & Z <= us[1], 1, ifelse(-us[1] < Z & Z <= us[2], 2, ifelse(-us[2] < Z & Z <= us[3], 3, ifelse(-us[3] < Z & Z <= us[4], 4, 5))))
#   Y = factor(Y, ordered=TRUE)
#   X = cbind(X, matrix(runif(n*(p-5)), nrow = n))
#   
#   return(data.frame(X, mu, Z, Y))
# }
# 
# K = 5 
# 
# us = c(0, 4, 8, 13)  # (us0 = -Inf), us1 = 0, us2 = 4, us3 = 8, us4 = 13, (us5 = Inf)
# 
# ## sample data
# datai = gen_data(n+ntest, us)
# y = as.numeric(datai$Y[1:n])
# X = datai[1:n, 1:p]
# 
# ytest = as.numeric(datai$Y[-c(1:n)])
# Xtest = datai[-c(1:n), 1:p]

## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set.seed(10)
RES = ordinal_BART(y = as.numeric(y), X = X, 
                   L = L, w0 = w0, s0 = s0, 
                   mcmc_iter = B/2, mcmc_bn = B/2, 
                   alpha = alpha, beta = beta, 
                   prob_move = prob_move)

var_importance_ordBART(RES, B/2, p)

b = sample(1:round(B/2), 1)
l = sample(1:L, 1)
RES[[b]]$T[[1]]$rules

RES[[b]]$T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf)
round(RES[[b]]$T[[l]]$mu, 3)
round(RES[[b]]$alpha, 3)
round(RES[[b]]$mu, 3)

v = RES[[b]]$T[[1]]$rules[-1]
as.integer(sub("^x([0-9]+)\\..*", "\\1", v))

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## PREDICTION

yprob_ordbart = matrix(0, ntest, K)
tmp = matrix(0, B/2, ntest)
for (i in 1:ntest) {
  print("::::::::: iter:")
  print(i)
  g = pred_ordBART(RES, Xtest[i, ], B/2, y)
  
  yprob_ordbart[i, ] = g[[1]]
  tmp[, i] = g[[2]]
}

# yprob_ordbart2 = apply(Xtest, 1, function(x) pred_ordBART(RES, x, B/2, y)[[1]])
ypred_ordbart = apply(yprob_ordbart, 1, function(x) which.max(x))

cbind(round(yprob_ordbart, 3), ytest)
round(apply(round(yprob_ordbart, 3), 2, mean), 3)
table(apply(yprob_ordbart, 1, function(x) which.max(x)), ytest)

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit multinomial regression model

res_nnet = multinom(as.factor(y) ~ ., data = as.data.frame(X))
summary(res_nnet)
yprob_nnet = predict(res_nnet, newdata = as.data.frame(Xtest), type = "probs")
ypred_nnet = apply(yprob_nnet, 1, function(x) which.max(x))

round(yprob_nnet, 3)
round(yprob_ordbart, 3)

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit ordinal regression model

# res_polr = MASS::polr(as.factor(y) ~ ., data = as.data.frame(X), method = "probit")
# yprob_polr = predict(res_polr, newdata = as.data.frame(Xtest), type = "probs")
# ypred_polr = predict(res_polr, newdata = as.data.frame(Xtest), type = "class")


# res_clm = ordinal::clm(ordered(y) ~ ., data = as.data.frame(X), link = "logit")
# summary(res_clm)
# ypred_clm = predict(res_clm, newdata = as.data.frame(Xtest), type = "class")

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit Bayesian multinomial regression model

# res_mnl = MCMCpack::MCMCmnl(y ~ ., data = data.frame(y, X),
#                             baseline = 1,          
#                             burnin = 1000, mcmc = 5000, thin = 5,
#                             verbose = 1000)

library(brms)
res_brms = brm(as.factor(y) ~ ., data = data.frame(y, X), family = categorical())
yprob_brms = fitted(res_brms, newdata = data.frame(Xtest))[, 1, ]
# yprob_brms = cbind(yprob_brms[, 1, 1], yprob_brms[, 1, 2], yprob_brms[, 1, 3], yprob_brms[, 1, 4])
ypred_brms = apply(yprob_brms, 1, function(x) which.max(x))

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit Bayesian ordinal regression model

library(brms)
res_brms_ord = brm(ordered(y) ~ ., data = data.frame(y, X), family = cumulative("logit"))
yprob_brms_ord = fitted(res_brms_ord, newdata = data.frame(Xtest))[, 1, ]
# yprob_brms_ord = cbind(yprob_brms_ord[, 1, 1], yprob_brms_ord[, 1, 2], yprob_brms_ord[, 1, 3], yprob_brms_ord[, 1, 4])
ypred_brms_ord = apply(yprob_brms_ord, 1, function(x) which.max(x))


source("C:/Users/HP/Dropbox/Assegno_unifi/BOPR.R")

library(SoftBart)  # for BART
opts = Opts(num_burn = 500, num_save = 500)

set.seed(7)
fitted_bopr = bopr(y ~ ., as.data.frame(Y = y, x = X), as.data.frame(Y = ytest, x = X_test), N = (opts$num_burn + opts$num_save), burn.in = opts$num_burn)

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit BART multinomial regression model

res_mbart = BART::mbart(x.train = X,        # Training predictors
                        y.train = y,        # Training response (factor with 2+ levels)
                        x.test = Xtest,     # Optional: test predictors
                        type = 'lbart',     # Model type: 'pbart' (probit), 'wbart' (weighted probit), 'lbart' (logit)
                        ntree = 50,         # Number of trees
                        ndpost = 500,      # Number of posterior samples
                        nskip = 500,        # Number of burn-in iterations
                        keepevery = 1)      # Thinning interval
                        # mc.cores = 2,       # Number of CPU cores for parallel processing
                        # seed = 123)         # Random seed for reproducibility
tmp = apply(res_mbart$prob.test, 2, mean)
yprob_mbart = matrix(0, ntest, K)
for (ix in 1:ntest)  yprob_mbart[ix, ] = tmp[(K*(ix-1)+1):(K*ix)]
ypred_mbart = apply(yprob_mbart, 1, function(x) which.max(x))

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit BART ordinal regression model

library(caret)
library(SoftBart)  # for BART
library(truncnorm) # for truncated normal dstn
library(progress)  # for progress bar

source("C:/Users/HP/Dropbox/Assegno_unifi/OPBART.R")

set.seed(7)
opts = Opts(num_burn = 500, num_save = 500)
res_opbart = opbart(ordered(as.factor(y)) ~ ., data.frame(y, X), data.frame(y = ytest, Xtest), opts = opts)
yprob_opbart = res_opbart$test_probs
ypred_opbart = apply(yprob_opbart, 1, function(x) which.max(x))

## :::::::::::::::::::::::::::::::::::::::::::::::
## Fit ordinal Random Forest

library(ordinalForest) 

set.seed(7)
res_orf = ordfor(depvar = "y", data = data.frame(y = ordered(as.factor(y)), X), perffunction = "probability")
yprob_orf = predict(res_orf, data.frame(y = ordered(as.factor(ytest)), Xtest))$classprobs
ypred_orf = apply(yprob_orf, 1, function(x) which.max(x))

## :::::::::::::::::::::::::::::::::::::::::::::::
## RESULTS

r1 = get_res(ypred_ordbart, yprob_ordbart, ytest)
r2 = get_res(ypred_nnet, yprob_nnet, ytest)
# r3 = get_res(ypred_polr, yprob_polr, ytest)
r4 = get_res(ypred_brms, yprob_brms, ytest)
r5 = get_res(ypred_brms_ord, yprob_brms_ord, ytest)
r6 = get_res(ypred_mbart, yprob_mbart, ytest)
r7 = get_res(ypred_opbart, yprob_opbart, ytest)
r8 = get_res(ypred_orf, yprob_orf, ytest)

rps(ytest, yprob_ordbart)
rps(ytest, yprob_opbart)
rps(ytest, yprob_orf)

r1$accuracy
r2$accuracy
r3$accuracy

r1$cross_entropy
r2$cross_entropy
r3$cross_entropy

r1$multi_auc
r2$multi_auc
r3$multi_auc

r1$F1_macro
r2$F1_macro
r3$F1_macro

mean(r1$TPR)
mean(r2$TPR)
mean(r3$TPR)

mean(r1$FPR)
mean(r2$FPR)
mean(r3$FPR)

