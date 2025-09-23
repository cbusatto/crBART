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
