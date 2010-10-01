`sample.plot` <-
function(data, labels=NULL, plot=TRUE, ...){

  ##Have to add the name for apply function
  data2<-data.frame(matrix(NA, ncol=ncol(data)+1, nrow=nrow(data)))
  data2[,1]<-rownames(data)
  data2[,2:ncol(data2)]<-data
  colnames(data2)<-c("probe", colnames(data))
  data<-data2
  rm(data2)

  ##For each variable
  out<-apply(data,1,exp_grp<-function(x){
    
    id<-x[1]
    x<-as.numeric(x[-1])
    s.name<-colnames(data)[-1]
    if (!is.null(labels)){
    
      x.g<-x.name<-x.col<-x.leg<-c()
      col<-rainbow(n=length(unique(labels)))
      cpt<-1

      for (i in unique(labels)){
        tmp<-x[which(labels==i)]
        x.g<-c(x.g,tmp)
        x.name<-c(x.name, s.name[which(labels==i)])
        x.col<-c(x.col, rep(col[cpt],length(tmp)))
        x.leg<-c(x.leg,i)
        cpt<-cpt+1
      }

      if (plot){
        x11()
        barplot(x.g,names.arg=x.name,axis.lty=1,las=2,col=x.col, ylab='Expression Level',main = paste("Expression level of",id), cex.axis=0.7, cex.lab=0.7, cex.main=0.8, cex.names=0.6, ylim=c(0,(max(x,na.rm=TRUE)+1)),...)
        legend(x=-1,y=(max(x,na.rm=TRUE)+1),x.leg,fill=col, bty="n", cex=0.7, text.col="gray50")
      }
    }
    else{
      if (plot){
        x11()
        barplot(x,names.arg=s.name, axisnames=TRUE,las=2,col="blue", ylab='Expression Level',main = paste("Expression level of",id), cex.axis=0.7, cex.lab=0.7, cex.main=0.8, cex.names=0.6, ...)
      }
    }
  })

  return (out)
}

