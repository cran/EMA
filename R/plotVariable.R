plotVariable <-
    function(acp, axes=c(1,2), new.plot=FALSE, lab=NULL, lim.cos2.var=0, palette="rainbow", ...) {

        if (!is.null(lab)){
            col.lab <- as.colors(lab, palette=palette)
        }
        else {
            col.lab <- rep(1,dim(acp$var$coord)[1])
        }
        
        plot.PCA(acp, choix="var", axes=axes, title=paste("Variable Representation on axes", axes[1], "and", axes[2]) , new.plot=new.plot, col.var=col.lab, lim.cos2.var=lim.cos2.var, cex=0.8, ...)
        
        if (!missing(lab)) {
            legend("topright", legend = unique(lab), col=unique(col.lab), pch=20, bty="n", cex=0.7, text.col="gray50")
        }
    }
