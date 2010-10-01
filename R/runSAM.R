runSAM <-
function(data, labels, nbpermut=500, q=0.05, plot=TRUE, method="d.stat", var.equal=TRUE, include.zero=FALSE, paired=FALSE, seed=123){
  
    ## pour mettre labels en 0 et 1 si non paired
    if(paired == FALSE){
        labels <- as.integer(factor(labels))-1
    }
    
    
    print("SAM analysis is running...")
    if(method=="d.stat")
        output <- sam(data, cl=labels, B=nbpermut, method="d.stat", var.equal=var.equal, include.zero=include.zero, rand=seed, control=samControl(delta=seq(0.1,10,0.05), lambda=0.5), med=TRUE)
    if(method=="wilc.stat")
        output <- sam(data, cl=labels, method="wilc.stat")
    if(method=="chisq.stat")
        output <- sam(data, cl=labels, method="chisq.stat", B=nbpermut, rand=seed)
    if(method=="trend.stat")
        output <- sam(data, cl=labels, method="trend.stat", B=nbpermut, rand=seed)
    
        
    print("Create table for all genes...")
    if(length(unique(labels))==2 & (method=="d.stat" || method=="wilc.stat")){
        alllist <- data.frame(rownames(data), output@d, output@p.value, output@fold)## output@q.value, output@fold)
        colnames(alllist)<-c("probeID", "Stat", "RawpValue", "FoldChange")## "AdjpValue", "FoldChange")
    }else{
        alllist <- data.frame(rownames(data), output@d, output@p.value)##, output@q.value)
        colnames(alllist) <- c("probeID", "Stat", "RawpValue")##, "AdjpValue")
    }
    alllist$probeID=as.character(alllist$probeID)

    print(output@mat.fdr[which(output@mat.fdr[,'FDR']>0),c("Delta","p0","False","Called","FDR")])

    print("Find delta...")
    outdelta <- try(findDelta(output, fdr=q))
    if(class(outdelta) == "matrix"){
        delta <- outdelta[2,1] 
        print(paste("Delta : ", delta, sep=""))
        
        print("Create SAM plot...")
        if(plot){
            #plot(output)
            plot(output, delta)
        }
        
        alllist$Significant <- rep(FALSE, nrow(alllist))  
        if(outdelta[2,1]!=0){
            delta.sum <- summary(output, delta)
            ds.mat.sig <- delta.sum@mat.sig
            rows <- ds.mat.sig$Row
            print(paste("Find",length(rows),"significant genes ..."))
            alllist[rows, "Significant"] <- TRUE
        }
    }
    else
        print("Can't find delta with chosen FDR... No SAM plot... No Significant slots in out data.frame...")

    ##warning("The adjusted pvalue is a qvalue, and is not related to the significant genes found with the SAM's FDR")
    return(alllist)
}

