FDR.BH<-
  function (pval, TST=FALSE,q){

  m <- length(pval)
  signif <- rep(FALSE,m)

  ##  LSU de BH & BY (1995)
  ## m0/m under-estimation of the FDR
  resFDR<-mt.rawp2adjp(pval, proc=c("BH"))
  adjp <- resFDR$adjp[order(resFDR$index),]

  if (TST){
    ## Second step (TST)
    ## Estimation of m0 to control FDR at q level
    r<-mt.reject(adjp[,"BH"],(q/(1+q)))$r
    outFDR<-adjp[,"BH"]*((m-r)/m)
  }
  else{
    outFDR<-adjp[,"BH"]
  }

  return (outFDR)
}
