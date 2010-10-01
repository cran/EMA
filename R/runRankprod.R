`runRankprod` <-
function(data,labels,q=0.05, plot=TRUE){

  if (is.null(rownames(data))){
    rownames(data)<-1:nrow(data)
  }
  
  ori<-rep(1,ncol(data))
  RP.out<-RPadvance(data, as.numeric(as.vector(labels)), ori, num.perm=100, logged=TRUE, gene.names=rownames(data), plot=FALSE)
  res<-topGene(RP.out, cutoff=q)

  ##Table1
  p=as.vector(res$Table1[,"gene.index"])
  out1<-NULL

  if (length(p)>1){
    out1<-data.frame(res$Table1[,2],res$Table1[,4:5],res$Table1[,3])
    colnames(out1)<-c("RP/Rsum","AdjpValue","RawpValue","FC(class1/class2)")

  }
  else if (p==1){
    out1<-as.data.frame(matrix(c(res$Table1[1,2],res$Table1[1,4:5],res$Table1[1,3]),nrow=1))
    colnames(out1)<-c("RP/Rsum","AdjpValue","RawpValue","FC(class1/class2)")
    
  }   

  ##Table2
  p=as.vector(res$Table2[,"gene.index"])
  out2<-NULL
  
  if (length(p)>1){
    out2<-data.frame(res$Table2[,2],res$Table2[,4:5],res$Table2[,3])
    colnames(out2)<-c("RP/Rsum","AdjpValue","RawpValue","FC(class1/class2)")
  }
  else if (length(p)==1){
    out2<-as.data.frame(matrix(c(res$Table2[1,2],res$Table2[1,4:5],res$Table2[1,3]),nrow=1))
    colnames(out2)<-c("RP/Rsum","AdjpValue","RawpValue","FC(class1/class2)") 
  }
 
  out<-list()
  out$table1<-out1
  out$table2<-out2

  #Graphical plot
  if (plot){
    plotRP(RP.out, cutoff=q)
  }
 
return(out)
}

