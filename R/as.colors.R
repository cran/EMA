as.colors <-
function(x, col.na="#E6E6E6", palette="rainbow", ...){
  v2m=FALSE

  if (is.na(palette) || is.null(palette)){
    stop("Palette must be specified")
  }
  
  if (is.vector(x)){
    v2m=TRUE
    x <- matrix(x, nrow=1)
  }
 
  xfactor <- factor(as.vector(as.matrix(x)))
  col.sel <- do.call(palette,list(n=length(levels(xfactor)), ...))

  colorlab<-apply(x,1,function(xx){
    for (i in levels(factor(xx))){
      ind<-which(levels(xfactor)==i)
      xx[which(xx==i)] <- col.sel[ind]
    }
    return (xx)
  })

  
  colorlab[which(is.na(colorlab))] <- col.na
  
  if (v2m){
    return(as.vector(colorlab))
  }
  else if(is.matrix(x) && is.vector(colorlab)){
    colorlab<-as.matrix(colorlab,ncol=ncol(x), nrow=nrow(x))
  }
  else{
    return(t(colorlab))
  }
}
