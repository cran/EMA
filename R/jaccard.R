`jaccard` <-
function(mat){
  out <- matrix(NA, ncol=ncol(mat), nrow=ncol(mat))
  colnames(out) <- colnames(mat)
  rownames(out) <- colnames(mat)

  for (i in 1:(ncol(mat))){
    for (j in 1:(ncol(mat))){
      a <- mat[,i]
      b <- mat[,j]
      out[i,j] <- (sum(b, na.rm=TRUE) + sum(a, na.rm=TRUE) -2*sum(a & b, na.rm=TRUE))/(sum(a, na.rm=TRUE) + sum(b, na.rm=TRUE) - sum(a & b, na.rm=TRUE))
    }
  }
  return (out)
}

