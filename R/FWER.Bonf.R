FWER.Bonf <-
  function (pval){

  m <- length(pval)
  signif <- rep(FALSE,m)
  
  ## Bonferroni
  resFDR<-mt.rawp2adjp(pval, proc=c("Bonferroni"))
  adjp <- resFDR$adjp[order(resFDR$index),]

  return (adjp[,"Bonferroni"])
}
