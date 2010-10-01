`genes.selection` <-
  function(data, thres.diff, thres.num, probs=0.25){
    ## check the parameters
    if (missing(thres.diff) && missing(thres.num)){
      stop("** Stop. No method found.")
    }
    if (!missing(thres.diff) && !missing(thres.num)){
      stop("** Stop. Choose one of the two options - thres.diff or thres.num")
    }

    diff.quantile <- function(x, probs=1/length(x))
      {
        vv <- quantile(x, probs=c(probs, 1-probs), na.rm = TRUE)
        return(vv[2] - vv[1])
        
      }
    
    
    rangeValues<-apply(data, 1, diff.quantile, probs=probs)
    
    if (!missing(thres.diff)){
      ind<-which(rangeValues>=thres.diff)
      genesList<-rownames(data[ind,])    
    }
    else if (!missing(thres.num)){
      genesList<-names(sort(rangeValues, decreasing=TRUE)[1:thres.num])
    }
    return (genesList)
  }

