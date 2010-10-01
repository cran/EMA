runWilcox <-
  function(data,labels,typeFDR="FDR-BH",q=0.05,plot=TRUE){

  list.exp <- colnames(data)
  list.probe <-rownames(data)
  if (is.null(list.probe)){
    list.probe<-1:nrow(data)
  }
  print ("Launch Wilcoxon test")

  test<-mt.teststat.num.denum(data, labels, test="wilcoxon")
  n<-length(which(labels==0))
  m<-length(which(labels==1))
  
  print ("Calculate pval")
  espU<-(n*m)/2
  pval.test <- 2*pwilcox(espU-abs(test$teststat.num), n, m, lower.tail=TRUE)
  
  print ("Adjusted pval")

  W<-espU-test$teststat.num
  padj<-multiple.correction(pval.test,typeFDR)
  indexsignif<-which(padj<=q)

  ## Result
  out <- data.frame(list.probe, W, pval.test,padj)
  colnames(out)<-c("probeID", "Stat", "RawpValue", "AdjpValue")

  ## Graphical plot
  if (plot){
    col <- rep("black", length(W))
    col[indexsignif] <- "red"
    ## Empirical distribution
    qqplot(W,rwilcox(length(W), m=m, n=n), col=col[order(W)], main=paste("QQplot. n=",n,"m=",m), xlab="X")
  }

  return (out)
}
