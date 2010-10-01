`dice` <-
function(mat) {
  out <- matrix(NA, ncol=ncol(mat), nrow=ncol(mat))
  colnames(out) <- colnames(mat)
  rownames(out) <- colnames(mat)

  for (i in 1:(ncol(mat))){
    for (j in 1:(ncol(mat))){
      a <- mat[,i]
      b <- mat[,j]
      out[i,j] <- (sum(abs(b), na.rm=TRUE) + sum(abs(a), na.rm=TRUE) -2*sum(a & b, na.rm=TRUE))/(sum(abs(a), na.rm=TRUE) + sum(abs(b), na.rm=TRUE))
    }
  }
  return (out)
}

