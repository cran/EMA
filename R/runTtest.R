runTtest <-
  function(data,labels,typeFDR="FDR-BH",algo="t", q=0.05, plot=TRUE){

  list.exp <- colnames(data)
  list.probe <-rownames(data)
  if (is.null(list.probe)){
    list.probe<-1:nrow(data)
  }
  
  ## Student / Welch
  print (paste("Launch ",algo," test"))
  test<-mt.teststat(data,labels, test=algo)

  ## ddl
  if (algo == "t.equalvar"){
    s1<-apply(data[,which(labels==0)],1,function(x){return (length(which(!is.na(x))))})	
    s2<-apply(data[,which(labels==1)],1,function(x){return (length(which(!is.na(x))))})
    ddl <- s1+s2-2
  }
  else if(algo == "pairt"){
      s1<-apply(data[,which(labels==0)],1,function(x){return (length(which(!is.na(x))))})
      ddl <- s1-1  
  }
  else if (algo == "t"){
    s1<-apply(data[,which(labels==0)],1,function(x){return (length(which(!is.na(x))))})	
    s2<-apply(data[,which(labels==1)],1,function(x){return (length(which(!is.na(x))))})
    w1<-apply(data[,which(labels==0)],1,var,na.rm=TRUE)/s1
    w2<-apply(data[,which(labels==1)],1,var,na.rm=TRUE)/s2
    
    ##Satterthwaite approximation
    ddl<-(w1+w2)^2/((w1^2/(s1-1))+(w2^2/(s2-1)))
  }

  ## pval
  print ("Calculate pval")
  pval.test <- 2*(pt(-abs(test),ddl))

  ## Multiple test correction
  print ("Adjusted pval")
  padj<-multiple.correction(pval.test,typeFDR)
  indexsignif<-which(padj<=q)
  print (paste(length(indexsignif), "significant genes."))

  out <- data.frame(list.probe, test, pval.test,padj)
  colnames(out)<-c("probeID",  "Stat", "RawpValue", "AdjpValue")

  if (plot){
    ## Graphical plot
    col <- rep("black", length(test))
    col[indexsignif] <- "red"
    qqnorm(test, col=col)
    qqline(test)
  }
  return (out)
}
