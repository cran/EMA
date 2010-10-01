plot3dSample <-
    function(acp, lab=NULL, palette="rainbow", ...) {  
        
        v <- acp$eig[1:3, 2]
        cp1 <- round(v[1], digit = 2)
        cp2 <- round(v[2], digit = 2)
        cp3 <- round(v[3], digit = 2)
        lab.x <- paste("Dimension ",1," (",cp1,"%)",sep="")
        lab.y <- paste("Dimension ",2," (",cp2,"%)",sep="")
        lab.z <- paste("Dimension ",3," (",cp3,"%)",sep="")
        
        
        if (is.null(lab)) {
            plot3d(acp$ind$coord[, 1], acp$ind$coord[, 2], acp$ind$coord[, 3], xlab=lab.x, ylab=lab.y, zlab=lab.z, col="red", size=4, main="Sample Representation", ...)
            
        } else {
            col <- as.colors(lab, palette=palette)
            plot3d(acp$ind$coord[, 1], acp$ind$coord[, 2], acp$ind$coord[, 3], xlab=lab.x, ylab=lab.y, zlab=lab.z, col=col, size=4, main="Sample Representation", ...)
            
        }
    }
