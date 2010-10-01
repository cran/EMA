multiple.correction <-
  function (pval,typeFDR,q){
  if (missing(pval))
    stop("Error : raw pvalues")
  if (missing(typeFDR))
    stop("Error : no typeFDR")
  if(!typeFDR %in% c("FWER","FDR-BH","FDR-TST","qvalue") )  
    stop("Error : unknown typeFDR")

  print (paste("typeFDR=", typeFDR))
  if (typeFDR == "FWER")
    padj <- FWER.Bonf(pval)

  else if (typeFDR == "FDR-BH")
    padj <- FDR.BH(pval)
  
  else if (typeFDR == "FDR-TST")
    padj <- FDR.BH(pval, TST=TRUE,q) 

  else if (typeFDR == "qvalue"){
    ##siggenes qvalue (version=2 for robust=TRUE)
    pi0 <- pi0.est(pval)$p0
    padj <- qvalue.cal(pval, pi0, version=2) 
  }

  return (padj)
}
