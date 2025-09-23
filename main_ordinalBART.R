# library(Rcpp)
library(data.tree)
# library(msos)
library(nnet)
library(pROC)

source("C:/Users/HP/Dropbox/Assegno_unifi/sample_data.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_Precm.R")
# sourceCpp("C:/Users/HP/Dropbox/Assegno_unifi/get_Precm.cpp")
# sourceCpp("C:/Users/HP/Dropbox/Assegno_unifi/solve_covariance.cpp")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_nsin.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_s.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/get_prior_mean.R")
# source("C:/Users/HP/Dropbox/Assegno_unifi/get_prior_mean_prob.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/GROW_step0.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/GROW_step.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/PRUNE_step.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/CHANGE_step.R")
# source("C:/Users/HP/Dropbox/Assegno_unifi/Gcopula.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/ordinal_BART.R")
# source("C:/Users/HP/Dropbox/Assegno_unifi/ordinal_BART_mod.R")
source("C:/Users/HP/Dropbox/Assegno_unifi/predict_ordinalBART.R")

# alpha = 0.95, beta = 2
## n = 100, p = 30, K = 4, B = 500, L = 50, w0 = 1, s0 = 0.3, seed 12 quasi ok
## n = 100, p = 30, K = 4, B = 500, L = 50, w0 = 0.5, s0 = 0.3, seed 12 ok
## n = 200, p = 30, K = 4, B = 500, L = 50, w0 = 1, s0 = 0.3, seed 12 ok
## n = 200, p = 30, K = 4, B = 500, L = 50, w0 = 0.5, s0 = 0.3, seed 12 cosi cosi


## Set seed for reproducibility
set.seed(12)

## Parameters
n = 100          # Number of observations
p = 150                  # Number of predictors
K = 4                # Number of categories (classes)
# p0 = rep(round(p/3), K-1)       # Number of true effects for the (K-1) classes
p0 = c(5, 5, 7)
# p0 = c(5, 7)


data = sim_data_multinom(n, p, p0, K)
table(data$y)

# data2 = sample_data_multinom(n, K, p, 5)
# table(data2$y)
# round(data2$probs, 3)
# data = sim_data_ordinal(n, p, p0, K)

y = as.numeric(data$y)
table(y)
X = data$X
beta0 = data$beta0
# prob = data$prob
table(data$y)

## set hyperparameters
B = 500
# w0 = sqrt(10) * 3
# ws = c(1, 1, 1)
alpha = 0.95
beta = 2
L = 100
w0 = 0.5


b = 0
mcmc_bn = 2
# prob_move = c(0.25, 0.25, 0.5)
prob_move = rep(1/3, 3)

# sample(1:3, 1, prob = c(0.5, 0.5, 0))

# set.seed(10)
# RES = ordinal_BART(y = y, X = X, L = L, w0 = w0, mcmc_iter = 77, mcmc_bn = B/s, alpha = alpha, beta = beta, prob_move = rep(1/3, 3))
# RES = ordinal_BART_mod(y = as.numeric(y), X = X, L = L, w0 = w0, ws = ws, mcmc_iter = B/2, mcmc_bn = B/2, alpha = alpha, beta = beta, prob_move = rep(1/3, 3))

set.seed(10)
RES = ordinal_BART(y = as.numeric(y), X = X, L = L, w0 = w0, s0 = 0.3, mcmc_iter = B/2, mcmc_bn = B/2, alpha = alpha, beta = beta, prob_move = prob_move)

# set.seed(10)
# RES = ordinal_BART(y = y, X = X, L = L, w0 = w0, mcmc_iter = 77, mcmc_bn = B/s, alpha = alpha, beta = beta, prob_move = rep(1/3, 3))
# RES2 = ordinal_BART_chol(y = as.numeric(y), X = X, L = L, w0 = w0, mcmc_iter = B/2, mcmc_bn = B/2, alpha = alpha, beta = beta, prob_move = rep(1/3, 3))

b = sample(1:round(B/2), 1)
l = sample(1:L, 1)
b
round(RES[[b]]$alpha, 3)

RES[[b]]$T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf)
round(RES[[b]]$T[[l]]$mu, 3)
round(RES[[b]]$alpha, 3)
round(RES[[b]]$mu, 3)

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## PREDICTION

## Simulate predictors 
set.seed(10)

ntest = 100
data_test = sim_data_multinom_test(ntest, p, K, beta0)
# data_test = sim_data_ordinal_test(ntest, p, K, beta0)
ytest = as.numeric(data_test$ytest)
Xtest = data_test$X


# ## Compute linear predictor for each class
# linear_predictor = Xtest %*% beta0
# 
# ## Convert to probabilities using softmax function
# exp_lp = exp(cbind(0, linear_predictor))  # Add baseline class (logit = 0)
# probs = exp_lp / rowSums(exp_lp)
# 
# ## Simulate response variable based on probabilities
# ytest = apply(probs, 1, function(prob) sample(1:K, size = 1, prob = prob))
# table(ytest)

# prediction = function(RES, Xtest, B, y) {
#   ntest = nrow(Xtest)
#   K = length(unique(y))  # or set K explicitly if known
#   
#   # run prediction for each test point
#   res_list = lapply(seq_len(ntest), function(i) {
#     g = pred_ordBART(RES, Xtest[i, ], B/2, y)
#     list(prob = g[[1]], tmp = g[[2]])
#   })
#   
#   # extract into matrices
#   yprob = do.call(rbind, lapply(res_list, `[[`, "prob"))
#   tmp   = do.call(cbind, lapply(res_list, `[[`, "tmp"))
#   
#   list(yprob = yprob, tmp = tmp)
# }
# 
# out = prediction(RES, Xtest, B, y)
# 
# out$yprob  # ntest x K
# out$tmp 


yprob = yprob_maj = matrix(0, ntest, K)
# ypred = rep(0, ntest)
tmp = tmp_maj = matrix(0, B/2, ntest)
for (i in 1:ntest) {
  print("::::::::: iter:")
  print(i)
  g = pred_ordBART(RES, Xtest[i, ], B/2, y)
  # g_maj = pred_majority_ordBART(RES, Xtest[i, ], B/2)
  
  yprob[i, ] = g[[1]]
  # yprob_maj[i, ] = g_maj[[1]]
  
  # ypred[i] = which.max(yprob[i, ])
  tmp[, i] = g[[2]]
  # tmp_maj[, i] = g_maj[[2]]
  
  # print(pred_ordBART(RES, Xtest[i, ]))
}

ypred = apply(yprob, 1, function(x) which.max(x))
cbind(round(yprob, 3), ytest)
round(apply(round(yprob, 3), 2, mean), 3)
table(apply(yprob, 1, function(x) which.max(x)), ytest)

round(yprob_nnet, 3)
round(yprob, 3)


# save(y, X, beta0, n, p, p0, K, w0, B, alpha, beta, L, RES, ytest, Xtest, yprob, ypred, file = "try1.RData")

# ypred = apply(tmp, 2, function(x) which.max(table(x)))
# ypred_maj = apply(tmp_maj, 2, function(x) which.max(table(x)))

b = 1
l = 5
b = sample(1:round(B/2), 1)
l = sample(1:L, 1)
b
round(RES[[b]]$alpha, 3)

RES[[b]]$T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf)
round(RES[[b]]$T[[l]]$mu, 3)
round(RES[[b]]$alpha, 3)
round(RES[[b]]$mu, 3)


RES[[b]]$T[[l]]$T_l
RES[[b]]$T[[l]]$T_l$"d1"$j
RES[[b]]$T[[l]]$T_l$leaves[[1]]$j
RES[[b]]$T[[l]]$T_l$leaves[[2]]$j
RES[[b]]$T[[l]]$T_l$leaves[[3]]$j


T[[10]]$T_l$Get(function(node) node$res, filterFun = isLeaf)
t(vs)
theta

T[[l]]$T_l$Get(function(node) node$res, filterFun = isLeaf)
T_star$Get(function(node) node$res, filterFun = isLeaf)
RES[[b]]$T[[l]]$theta

# Fit multinomial regression model
res_nnet = multinom(as.factor(y) ~ ., data = as.data.frame(X))
summary(res_nnet)
yprob_nnet = predict(res_nnet, newdata = as.data.frame(Xtest), type = "probs")
ypred_nnet = apply(yprob_nnet, 1, function(x) which.max(x))
yprob_nnet2 = predict(res_nnet, newdata = as.data.frame(Xtest), type = "class")


res_polr = MASS::polr(as.factor(y) ~ ., data = as.data.frame(X), method = "logistic")
yprob_polr = predict(res_polr, newdata = as.data.frame(Xtest), type = "probs")
ypred_polr = predict(res_polr, newdata = as.data.frame(Xtest), type = "class")


res_clm = ordinal::clm(as.factor(y) ~ ., data = as.data.frame(X), link = "logit")
summary(res_clm)
# summary(res_clm)
ypred_clm = predict(res_clm, newdata = as.data.frame(Xtest), type = "class")
print(pred_class)





cbind(cbind(ytest, ypred), ypred_nnet)
cbind(cbind(cbind(ytest, ypred), ypred_nnet), ypred_polr)








get_res = function(ypred, yprob, ytest) {
  
  l1 = length(table(ypred))
  l2 = length(table(ytest))
  
  # Accuracy
  accuracy = mean(ypred == ytest)
  cat("Accuracy:", accuracy, "\n")
  
  # Multiclass AUC
  TMP = as.matrix(yprob)
  colnames(TMP) = levels(as.factor(ytest))
  multi_auc = multiclass.roc(response = as.factor(ytest), predictor = TMP)
  cat("Multiclass AUC:", multi_auc$auc, "\n")
  
  # Cross-Entropy Loss
  # cross_entropy = -mean(rowSums(model.matrix(~ as.factor(ytest)-1) * log(yprob)))
  cross_entropy = 0
  for (i in 1:length(ytest)) {
    cross_entropy = cross_entropy - log(yprob[i, ytest[i]])
  }
  cross_entropy = cross_entropy / length(ytest)
  cat("Cross-Entropy Loss:", cross_entropy, "\n")
  
  # Confusion Matrix
  conf_matrix = table(ytest, ypred)
  cat("Confusion Matrix:\n")
  print(conf_matrix)
  
  # F1 Score (Macro and Micro Average)
  precision = diag(conf_matrix) / colSums(conf_matrix)
  recall = diag(conf_matrix) / rowSums(conf_matrix)
  f1_class = 2 * precision * recall / (precision + recall)
  f1_macro = mean(f1_class, na.rm = TRUE)
  f1_micro = sum(diag(conf_matrix)) / sum(conf_matrix)
  cat("F1 Score (Macro):", f1_macro, "\n")
  cat("F1 Score (Micro):", f1_micro, "\n")
  
  # TPR and FPR for each class
  tpr = diag(conf_matrix) / rowSums(conf_matrix)  # Sensitivity
  fpr = (colSums(conf_matrix) - diag(conf_matrix)) / (sum(conf_matrix) - rowSums(conf_matrix)) # False Positive Rate
  cat("TPR per class:\n")
  print(tpr)
  cat("FPR per class:\n")
  print(fpr)
  
  # return(list(accuracy = accuracy, multi_auc = multi_auc, cross_entropy = cross_entropy, 
  return(list(accuracy = accuracy, cross_entropy = cross_entropy, 
              F1_micro = f1_micro, F1_macro = f1_macro, TPR = tpr, FPR = fpr))
}

r1 = get_res(ypred, yprob, ytest)
r2 = get_res(ypred_nnet, yprob_nnet, ytest)
r3 = get_res(ypred_polr, yprob_polr, ytest)
r3 = get_res(ypred_polr, yprob_polr, ytest)

r1$accuracy
r2$accuracy
r3$accuracy

r1$cross_entropy
r2$cross_entropy
r3$cross_entropy

r1$F1_micro
r2$F1_micro
r3$F1_micro

r1$F1_macro
r2$F1_macro
r3$F1_macro

mean(r1$TPR)
mean(r2$TPR)
mean(r3$TPR)

mean(r1$FPR)
mean(r2$FPR)
mean(r3$FPR)

