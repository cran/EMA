`distrib.plot` <-
    function(data, labels=NULL, plot=TRUE, ...){
        if (is.vector(data)){
            data <- matrix(data, nrow=1)
        }
        
        ##Have to add the name for apply function
        data2<-data.frame(matrix(NA, ncol=ncol(data)+1, nrow=nrow(data)))
        data2[,1]<-rownames(data)
        data2[,2:ncol(data2)]<-data
        colnames(data2)<-c("probe", colnames(data))
        data<-data2
        rm(data2)
        
        ##For each variable
        apply(data,1,exp_grp<-function(x){
            x11()
            id<-x[1]
            x<-as.numeric(x[-1])

            if (!is.null(labels)){
                col<-rainbow(n=length(unique(labels)))
                cpt<-1
                
                maxcounts<-NULL
                for (i in unique(labels)){
                    x.tmp<-x[which(labels==i)]
                    if(cpt==1){
                        plot(density(x.tmp, na.rm=TRUE), col=col[cpt], type="l", lwd=2, cex.axis=0.7, cex.lab=0.7, cex.main=0.8,main = paste("Expression level of",id),xlab='Expression Level',ylab='Density', ...)
                    }
                    else{
                        points(density(x.tmp, na.rm=TRUE), col=col[cpt], type="l", lwd=2)
                    }

                    cpt<-cpt+1
                }
                legend(x=min(x),y=maxcounts,unique(labels),fill=col, bty="n", cex=0.7, text.col="gray50", ...)
            }
            else{
                plot(density(x, na.rm=TRUE), cex.lab=0.7, col="blue",cex.axis=0.7,cex.main=0.8,main = paste("Expression level of",id),xlab='Expression Level',ylab='Density', ...)
            }
        })
    }
