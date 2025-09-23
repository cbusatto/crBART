get_nsin = function(T) {
  
  ## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## function to get the number of singly internal nodes in a tree T
  
  ## INPUT:
  ## T: current tree
  
  D = T$leafCount
  
  d = 1
  nsin = 0
  
  while (d < D) {
    tmp_parent = T$leaves[[d]]$parent
    check1 = tmp_parent$children[[1]]$isLeaf
    check2 = tmp_parent$children[[2]]$isLeaf

    if ((check1 == TRUE) & (check2 == TRUE)) {
      nsin = nsin + 1
      d = d + 1
    }

    d = d + 1
  }
  
  
  # while (d < D) {
  #   
  #   str1 = as.numeric(stringr::str_sub(T$leaves[[d]]$name, 2, stringr::str_length(T$leaves[[d]]$name)))
  #   str2 = as.numeric(stringr::str_sub(T$leaves[[d+1]]$name, 2, stringr::str_length(T$leaves[[d+1]]$name)))
  #   
  #   if ((str2 - str1) == 1) {
  #     
  #     nsin = nsin + 1
  #     d = d + 1
  #   }
  #   
  #   d = d + 1
  # }
  
  return(nsin)
}
