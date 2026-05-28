# Ranked Probability Score (RPS)
RPS = function(ytest, P) {
  
  n = length(y); K = ncol(P)
  F = t(apply(P, 1, cumsum))        # n x K, row cumulative; last col = 1
  
  # keep only k=1..K-1
  Fk = F[, 1:(K-1), drop = FALSE]
  
  # binary truth indicators for thresholds
  O = sapply(1:(K-1), function(k) as.integer(y <= k))
  
  # squared differences, averaged
  mean(rowSums((Fk - O)^2) / (K-1))
}

get_res = function(ypred, yprob, ytest) {
  
  l1 = length(table(ypred))
  l2 = length(table(ytest))
  K = length(table(ytest))
  
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
  # Define the full set of possible classes
  if (length(table(ypred)) != K) {
    tmp = sort(unique(ytest))
    
    # Build the full confusion matrix with all levels
    conf_matrix = table(factor(ytest, levels = tmp),
                        factor(ypred, levels = tmp))
  } else {
    conf_matrix = table(ytest, ypred)
  }
  
  cat("Confusion Matrix:\n")
  print(conf_matrix)
  
  # F1 Score (Macro and Micro Average)
  precision = diag(conf_matrix) / colSums(conf_matrix)
  precision[is.na(precision)] = 0
  
  recall = diag(conf_matrix) / rowSums(conf_matrix)
  
  f1_class = 2 * precision * recall / (precision + recall)
  f1_class[is.na(f1_class)] = 0  # handle 0/0 cases
  
  f1_macro = mean(f1_class, na.rm = TRUE)
  # f1_micro = sum(diag(conf_matrix)) / sum(conf_matrix)
  cat("F1 Score (Macro):", f1_macro, "\n")
  # cat("F1 Score (Micro):", f1_micro, "\n")
  
  # TPR and FPR for each class
  tpr = diag(conf_matrix) / rowSums(conf_matrix)  # Sensitivity
  fpr = (colSums(conf_matrix) - diag(conf_matrix)) / (sum(conf_matrix) - rowSums(conf_matrix)) # False Positive Rate
  cat("TPR per class:\n")
  print(tpr)
  cat("FPR per class:\n")
  print(fpr)
  
  # # RPS (Ranked Probability Score)
  # rps = RPS(ytest, yprob)
  # cat("RPS:\n")
  # print(rps)
  
  return(list(accuracy = accuracy, multi_auc = multi_auc, cross_entropy = cross_entropy,
              F1_macro = f1_macro, TPR = tpr, FPR = fpr))
              # , RPS = rps))
}
