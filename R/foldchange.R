`foldchange` <-
function (data, labels, unlog=TRUE){
    if(missing(data))
        stop("** FAILURE : 'data' is missing **")
    if(missing(labels))
        stop("** FAILURE : 'labels' is missing **")

    data1.ave = apply(data[,which(labels==0)], 1, mean, na.rm=TRUE)
    data2.ave = apply(data[,which(labels==1)], 1, mean, na.rm=TRUE)
    fold.change = data2.ave - data1.ave
    
    if(unlog){
        fold.change <- 2^fold.change
    }
    
    return(fold.change)
}

