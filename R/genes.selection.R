`genes.selection` <-
  function(data, thres.IQR, thres.num){
    ## check the parameters
    if (missing(thres.IQR) && missing(thres.num)){
        stop("** Stop. No method found.")
    }
    if (!missing(thres.IQR) && !missing(thres.num)){
        stop("** Stop. Choose one of the two options - thres.IQR or thres.num")
    }

    
    iqrValues<-apply(data,1,IQR,na.rm=TRUE)
    
    if (!missing(thres.IQR)){
        ind<-which(iqrValues>=thres.IQR)
        genesList<-rownames(data[ind,])    
    }
    else if (!missing(thres.num)){
        genesList<-names(sort(iqrValues, decreasing=TRUE)[1:thres.num])
    }
    return (genesList)
}

