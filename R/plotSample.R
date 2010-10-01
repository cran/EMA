plotSample <-
    function(acp, axes=c(1,2), new.plot=FALSE, lab="quality", palette="rainbow", lim.cos2.sample=0, text=TRUE, lab.title=NULL, ellipse=FALSE, ...) {
    
    nind <- dim(acp$ind$coord)[1]
    col.PC <- rep("black", nind)
    qual <- qualitySample(acp, axes=axes)[,3]
    
    if (!is.null(lab)) {
        if ((length(lab) == 1)&&(lab == "quality")) {
            col.PC <- rep("grey", nind)
            col.PC[qual < 0.3] <- "red"
            col.PC[qual > 0.7] <- "green"
        }
        else{
            col.PC <- as.colors(lab, palette=palette)
            corresColLab <- unique(col.PC)
            names(corresColLab) <- unique(lab)
        }
    }
    
    col.PC[qual < lim.cos2.sample] <- "white" 
    if (text){
        text="all"
    }else{
        text=NULL
    }
    
    
    
    plot.PCA(acp, axes=c(axes[1],axes[2]), col.ind=col.PC, title=paste("Sample representation \n Axes", axes[1], "and", axes[2]), new.plot=new.plot, label=text, cex=0.8, ...)
    
    ## ellipses have to be plotted by hand
    if ((ellipse) && !is.null(lab) && (lab != "quality")){
        aux <- cbind.data.frame(lab, acp$ind$coord[, axes])
        coord.ell <- coord.ellipse(aux, bary=TRUE)$res
        nbre.ellipse <- nlevels(coord.ell[, 1])
        for (e in 1:nbre.ellipse) {
            data.elli <- coord.ell[coord.ell[, 1] == levels(coord.ell[, 1])[e], -1]
            lines(x = data.elli[, 1], y = data.elli[, 2], col = corresColLab[levels(coord.ell[, 1])[e]])
        }
    }
    
    
    if (!is.null(lab)) {
        if ((length(lab) == 1)&&(lab == "quality")) {
            legend("topright", paste("Quality index", c("<0.3",">=0.3 and <=0.7",">0.7")), col=c("red","grey","green"), pch=20, title="Quality", bty="n", cex=0.7, text.col="gray50")
        } else {
            legend("topright", legend=unique(lab), col=unique(col.PC), pch=20, title=lab.title, bty="n", cex=0.7, text.col="gray50")
        }
    }
}
