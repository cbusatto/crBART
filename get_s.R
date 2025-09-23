get_s = function(y, K) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## function to create table(y)
  
  ## INPUT:
  ## y: vector of observed values (1,...,K)^n
  ## K: number of possible outcomes in y
  
  if (length(y) != 0) {
    ym = matrix(0, length(y), K)
    # colnames(ym) = 1:K
    for (i in 1:length(y)) ym[i, y[i]] = 1
    out = apply(ym, 2, sum)
  } else {
    out = rep(0, K)
  }
  
  return(out)
}