probePlots <-
    function(abatch, path, pbsList, labAxisProbes=TRUE, labAxisArrays=TRUE, legendArrays=TRUE, legendProbes=TRUE, cex.axis=0.9, cex.legend=0.8, pdfName)
{

        ## check if arguments are ok
        if(missing(abatch) && missing(path))
            stop("** ERROR : 'abatch' or 'path' arguments are empty")
        
        if(missing(pbsList))
            stop("** ERROR : 'pbsList' argument is missing")
        if(class(pbsList)!="character")
            stop("** ERROR : 'pbsList' argument is not a character class")
        
        if(!missing(abatch)){
            if(class(abatch) ==  "AffyBatch")
                abatch <- abatch
            else
                stop("** ERROR : 'abatch' argument is not an AffyBatch class")
        }
        
        if(!missing(path)){
            if(class(path) ==  "character")
                abatch <- ReadAffy(celfile.path=path)
            else
                stop("** ERROR : 'path' argument is not a character class")
        }

        if(!missing(pdfName)){
          if(class(pdfName) ==  "character")
            pdf(pdfName, width=9, height=8)
          else
            stop("** ERROR : 'pdfName' argument is not a character")
        }


        ## select the probesets in the affybatch
        pbs <- probeset(abatch, pbsList)
        
        ## prepare the plot
        par(mfrow=c(2,3), bg = "gray80", mai=c(5, 4, 4, 2), mar=c(5, 4, 4, 2), oma=c(0,0,2,0))
        
        for(i in 1:length(pbs)){
            pmi <- pm(pbs[[i]])
            mmi <- mm(pbs[[i]])
            
            ## INTER CHIPS PLOT
            plot(pbs[[i]], type="l", main=paste("Probeset ", names(pbs)[i], "\nInter arrays\nxaxis = probes", sep=""), axes=FALSE, ylab="Intensity", xlab="", col=rainbow(nrow(pmi)), lty=1:5)
            if(!is.null(rownames(pmi)) && labAxisProbes){
              xlabels <- rownames(pmi)
              axis(1, 1:nrow(pmi), xlabels, las=2, cex.axis=cex.axis)
            }
            if(is.null(rownames(pmi)) && labAxisProbes){
              xlabels <- as.character(1:nrow(pmi))
              axis(1, 1:nrow(pmi), xlabels, las=2, cex.axis=cex.axis)
            }
            if(!labAxisProbes)
              axis(1, 1:nrow(pmi), labels=FALSE, las=2, cex.axis=cex.axis)
            axis(2)
            box()
            if(!is.null(colnames(pmi)) && legendArrays)
              legend("topright", colnames(pmi), lty=1:5, col=rainbow(nrow(pmi)), cex=cex.legend)
            
            ## INTER PROBES PLOT - PM
            matplot(t(pmi), type="l", main=paste("Probeset ", names(pbs)[i], "\nInter probes - Perfect Match\nxaxis = arrays", sep=""), axes=FALSE, ylab="Intensity", xlab="", col=rainbow(nrow(pmi)), lty=1:5)
            if(!is.null(colnames(pmi)) && labAxisArrays){
              xlabels <- colnames(pmi)
              axis(1, 1:ncol(pmi), xlabels, las=2, cex.axis=cex.axis)
            }
            if(is.null(colnames(pmi)) && labAxisArrays){
              xlabels <- as.character(1:ncol(pmi))
              axis(1, 1:ncol(pmi), xlabels, las=2, cex.axis=cex.axis)
            }
            if(!labAxisArrays)
              axis(1, 1:ncol(pmi), labels=FALSE, las=2, cex.axis=cex.axis)
            axis(2)
            box()
            if(!is.null(rownames(pmi)) && legendProbes)
              legend("topright", rownames(pmi), lty=1:5, col=rainbow(nrow(pmi)), cex=cex.legend)
            if(is.null(rownames(pmi)) && legendProbes)
              legend("topright", as.character(1:nrow(pmi)), lty=1:5, col=rainbow(nrow(pmi)), cex=cex.legend)
            
            ## INTER PROBES PLOT - MM
            matplot(t(mmi), type="l", main=paste("Probeset ", names(pbs)[i], "\nInter probes - Mis Match\nxaxis = arrays", sep=""), axes=FALSE, ylab="Intensity", xlab="", col=rainbow(nrow(mmi)), lty=1:5)
            if(!is.null(colnames(mmi)) && labAxisArrays){
              xlabels <- colnames(mmi)
              axis(1, 1:ncol(mmi), xlabels, las=2, cex.axis=cex.axis)
            }
            if(is.null(colnames(mmi)) && labAxisArrays){
              xlabels <- as.character(1:ncol(mmi))
              axis(1, 1:ncol(mmi), xlabels, las=2, cex.axis=cex.axis)
            }
            if(!labAxisArrays)
              axis(1, 1:ncol(mmi), labels=FALSE, las=2, cex.axis=cex.axis)
            axis(2)
            box()
            if(!is.null(rownames(mmi)) && legendProbes)
              legend("topright", rownames(mmi), lty=1:5, col=rainbow(nrow(mmi)), cex=cex.legend)
            if(is.null(rownames(mmi)) && legendProbes)
              legend("topright", as.character(1:nrow(mmi)), lty=1:5, col=rainbow(nrow(mmi)), cex=cex.legend)
            
        }

        if(!missing(pdfName)){
          if(class(pdfName) ==  "character")
            dev.off()
        }

      }
