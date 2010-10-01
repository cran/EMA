
clustering.plot <- function(tree, tree.sup, data=NULL, lab, lab.sup, dendro=TRUE, dendro.sup=TRUE, title="", scale="row", heatcol, pdfname, names=TRUE, names.sup=TRUE, names.dist=TRUE, trim.heatmap=1, palette="rainbow", legend=TRUE, legend.pos="topright", ...){

  if (missing(tree)){
    stop("At least, one tree is required.")
  }

  if ((!missing(tree.sup) && missing(data))){
    stop("Two-ways clustering need two trees and data matrix.")
  }
  ##check that label exist
  if(! "order.lab" %in% names(tree))
      stop("order.lab missing in tree object. Check you data colnames.") 

  if (!missing(lab) && is.vector(lab)){
    lab <- matrix(lab, ncol=1)
  }
  if (!missing(lab) && is.data.frame(lab)){
      for (i in 1:ncol(lab)) lab[, i] <- as.character(lab[, i])
      lab <- as.matrix(lab)
  }
  
  if (!missing(lab.sup) && is.vector(lab.sup)){
    lab.sup <- matrix(lab.sup, ncol=1)
  }
  if (!missing(lab.sup) && is.data.frame(lab.sup)){
      for (i in 1:ncol(lab.sup)) lab.sup[, i] <- as.character(lab.sup[, i])
      lab.sup <- as.matrix(lab.sup)
  }
  
  if (!missing(pdfname)){
    pdf(pdfname)
  }
  
  op <- par(no.readonly = TRUE)
  #op<-par(mar=c(4,4,3,3)) 
  on.exit(par(op))

  ## Heatmap display
  if (!is.null(data)){
    if (!is.matrix(data)){
      data <- as.matrix(data)
    }

    #par(mar=c(5,4,3,9))
    mean.ori<-mean(data, na.rm=TRUE)
    
    ##Always center (to -1/1 scaling)
    data <- data-mean(data, na.rm=TRUE)
    
    ##Heatmap
    if (scale=="row"){
      data=t(scale(t(data)))
      mean.ori<-0
      if (tree$method=="euclidean"){
        warnings("We do not advice to scale the data, when you use an euclidean distance ...")
      }
    }
    else if (scale=="column"){
      data=scale(data)
      mean.ori<-0
    }
    
    ##Trimming
    q <- quantile(data, c((1-trim.heatmap), trim.heatmap), na.rm=TRUE)
    data[data < q[1]] = q[1]
    data[data > q[2]] = q[2]
    
    ##-1/1
    maxi <- max(data, na.rm=TRUE)
    mini <- min(data, na.rm=TRUE)
    data[!is.na(data) & data > 0] <- data[!is.na(data) & data > 0]/maxi
    data[!is.na(data) & data < 0] <- -data[!is.na(data) & data < 0]/mini
    maxi.real <- q[2]+mean.ori
    mini.real <- q[1]+mean.ori

    ##RowSideColors
    if (!missing(lab.sup)){
      rowcolorlab<-t(as.colors(t(lab.sup), palette=palette))
      if (ncol(lab.sup)==1){
        rowcolorlab<-cbind(rowcolorlab,rowcolorlab)
      }
      rownames(rowcolorlab)<-NULL
    }
    ##ColSideColors
    if (!missing(lab)){
      colcolorlab<-as.colors(lab, palette=palette)
      if (ncol(lab)==1){
        colcolorlab<-cbind(colcolorlab,colcolorlab)
      }
      colnames(colcolorlab)<-NULL
    }

    ##Dendrograms
    d.t <- d.ts <-NA
    if (dendro){
      d.t <- as.dendrogram(as.hclust(tree))
    }
    if (dendro.sup & !missing(tree.sup)){
      d.ts <- as.dendrogram(as.hclust(tree.sup))
    }

    ##Margins
    bmar=floor(max(nchar(colnames(data))/2))
    rmar=floor(max(nchar(rownames(data))/2))

    ##Text & legend
    n <- n.sup <-  NULL
    if (!names){
      n <- rep("", length(tree$order))
    }
    if (!names.sup){
        if (!missing(tree.sup)){
            n.sup <- rep("", length(tree.sup$order))
        }
        else{
            n.sup <- rep("", nrow(data)) 
        }
    }
   
    if (names.dist){
      xtitle <- paste("Hierarchical Clustering :",tree$method,attr(tree$diss,"Metric"))
      bmar <- bmar+3
      if (dendro.sup & !missing(tree.sup)){
        ytitle <- paste("Hierarchical Clustering :",tree.sup$method,attr(tree.sup$diss,"Metric"), sep=" ")
        rmar <- rmar+3
      }else{
        ytitle <- NA
      }
    }
    else{
      xtitle <- NA
      ytitle <- NA
    }
    if (missing(heatcol)) 
      heatcol <- myPalette(low = "green", high = "red", mid ="black", k=50)
    
    ##Heatmap
    if (!missing(lab) && !missing(lab.sup)){
        heatmap.plus(data, Rowv=d.ts, Colv=d.t, margins=c(bmar,rmar), cexCol=0.9, cexRow=0.9, scale="none", cex.main=0.7, main=title, ylab=ytitle, xlab=xtitle, RowSideColors=rowcolorlab,ColSideColors=colcolorlab, labRow=n.sup, labCol=n, col=heatcol,...)
    }
    else if (!missing(lab.sup)){
        heatmap.plus(data, Rowv=d.ts, Colv=d.t, margins=c(bmar,rmar), cexCol=0.9, cexRow=0.9, scale="none",  cex.main=0.7, main=title, ylab=ytitle, xlab=xtitle, RowSideColors=rowcolorlab, labRow=n.sup, labCol=n, col=heatcol,...)
    }
    else if (!missing(lab)){
      heatmap.plus(data, Rowv=d.ts, Colv=d.t, margins=c(bmar,rmar), cexCol=0.9, cexRow=0.9, scale="none", cex.main=0.7,  main=title, ylab=ytitle, xlab=xtitle, ColSideColors=colcolorlab, labRow=n.sup, labCol=n, col=heatcol,...)
    }
    else{
        heatmap.plus(data, Rowv=d.ts, Colv=d.t,  scale="none", cexCol=0.9, cexRow=0.9, cex.main=0.7, main=title,  ylab=ytitle, xlab=xtitle, labRow=n.sup, labCol=n, col=heatcol, margins=c(bmar,rmar),...)
    }

    ##Label legend
    if (legend && (!missing(lab) || !missing(lab.sup))){
      op1 <- par(fig = c(0.05,0.2,0.78,0.99), new=TRUE)
      op2 <- par(mar=c(0,0,0,0))
      lall<-c()
      if (!missing(lab)){
         lall<-c(lall, as.vector(lab))
      }
      if (!missing(lab.sup)){
         lall<-c(lall, as.vector(lab.sup))
      }
      leg <- lall
      if (is.element(NA,leg)){
          leg[which(is.na(leg))] <- "NA"
      }
         
      legend(legend.pos,legend=unique(leg), fill=unique(as.colors(lall, palette=palette)), bty="n", cex=0.7, text.col="gray50")
      par(op1)
      par(op2)
    }
    ##Heatmap legend
    if (legend){
      op1 <- par(fig = c(0.93,0.95,0.8,0.95), new=TRUE)
      op2 <- par(mar=c(0,0,0,0))
      iml<-matrix(seq(from=mini.real, to=maxi.real, length.out=10), nrow=1)
      image(x=1, y=1:10, z=iml, axes=FALSE, ylab="", xlab="", col=heatcol, ...)
      if (trim.heatmap != 1){
	axis(side=4, labels=c(paste("<",round(mini.real,1)),round((mini.real + maxi.real)/2,1),paste(">",round(maxi.real,1))), at=c(1,5.5,10), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE,line=-0.7, ...)
      }else{
	axis(side=4, labels=c(round(mini.real,1),round((mini.real + maxi.real)/2,1),round(maxi.real,1)), at=c(1,5.5,10), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE,line=-0.7, ...)
      }
      par(op1)
      par(op2)
    }
  }
  
  ##One tree
  else{
    
    bmar=floor(max(nchar(tree$order.lab))/3.5)
    par(mar=c(bmar,4,3,3))
    if(!missing(lab)){
      pl=0.1+0.1*ceiling(ncol(lab)/3)
      layout(matrix(c(1:2), 2, 1, byrow=TRUE), heights=c(1-pl,pl))#, respect=TRUE)
    }
    else{
	layout(matrix(c(1:2), 2, 1, byrow=TRUE), heights=c(0.9,0.1))#, respect=TRUE)
    }
    
    plot(as.dendrogram(as.hclust(tree)), edge.root=TRUE, nodePar=list(col=c('blue','black')), edgePar=list(col=c('blue', 'black')), leaflab="none", dLeaf=NULL,horiz=FALSE, frame.plot=FALSE, ylab=paste(tree$method,"linkage"), cex.lab=0.8, main=title, cex.main=0.8, ...)
    parmar=par("mar")
    paru=par("usr")

  if (legend && !missing(lab)){
      leg <- lab
      if (is.element(NA,leg)){
          leg[which(is.na(leg))] <- "NA"
      }
      
      legend(legend.pos,legend=unique(as.vector(leg)), fill=unique(as.colors(as.vector(lab), palette=palette)), bty="n", cex=0.7, text.col="gray50")
    }

    if (names){
      mtext(text=tree$order.lab, at=1:length(tree$order), side=1, line=1, col="black", las=2, cex=0.65, font=1, ...)
    }

    mtext(paste("Hierarchical Clustering :",tree$method,attr(tree$diss,"Metric"), sep=" "), at=length(tree$order)/2, side=3, las=1, col="gray50", cex=0.8)
    
    if (!missing(lab)){
      lab.num <- lab
      ind <- 0
      for (i in unique(as.vector(lab))){
        if (is.na(i)){
            lab.num[which(is.na(lab))]<-ind
        }
        else{
            lab.num[which(lab==i)]<-ind
        }
        ind <- ind+1
      }
      im<-matrix(as.numeric(lab.num),ncol=ncol(lab), nrow=nrow(lab))
      im <- as.matrix(im[tree$order,])
      
      par(mar=c(2,parmar[2],parmar[3],parmar[4]))
      image(x=1:length(tree$order), y=1:ncol(im),xlim=c(paru[1],paru[2]), z=im,axes=FALSE, ylab="", xlab="", col=unique(as.colors(as.vector(lab),palette=palette)),...)
  
      if (ncol(lab)>1){
        axis(side=4, labels=colnames(lab), at=1:ncol(im), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE, ...)
      }
      else{
        axis(side=4, labels=FALSE, at=1:ncol(im), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE, ...)
      }

    }
  }
   
  if (!missing(pdfname)){
    dev.off()
  }


 }
