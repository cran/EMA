qualitySample <-
function(acp, axes = c(1:3)){

  ## Goal
  ## Sample representation quality : compute cos2 between samples and each combination of 2 axes
  ## Representation quality for 2 axes is the sum of representation quality for each axe

  ## Input
  ## acp : result from PCA function
  ## axes: vector of integer, axes number for quality computation

  ## Output
  ## cos2 between samples and axes

  qual <- data.frame(acp$ind$cos2)[, axes]
  name <- colnames(qual)
  ncomp <- length(axes)
  k <- ncomp+1

  for (i in 1:(ncomp-1)) {
    for (j in (i+1):ncomp) {
      qual[,k] <- qual[,i]+qual[,j]
      name <- c(name, paste(name[i], "-", name[j], sep=""))
      k <- k+1
    }
  }
  names(qual) <- name  
  return(qual)
}

