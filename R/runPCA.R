runPCA <-
    function(X, ncp=5, scale=TRUE, ind.sup=NULL, quanti.sup=NULL, quali.sup=NULL, sample.qual=TRUE, variable.qual=FALSE, sample.cont=TRUE, variable.cont=FALSE, plotSample=TRUE, plotVariable=FALSE, plotInertia = TRUE, plotBiplot=FALSE, plot3dSample=FALSE, lab.sample="quality", lab.var=NULL, palette="rainbow", lim.cos2.sample=0, lim.cos2.var=0, pdf=FALSE, pdfname= NULL, verbose=TRUE, ...) {
        
        ## lot of options here ! try with ... ?
        ##  if (class(X) == "ExpressionSet") X <- exprs(X)
        if (!is.null(pdfname)) {
            pdf <- TRUE
        }
        
        if (pdf){
            if (!is.null(pdfname))
                pdf(pdfname, width=8, height=8, pointsize=8)
            else {
                pdf(width=8, height=8, pointsize=8)
            }
            op <- par(mfrow=c(2,1))
        }
        
        acp <- PCA(X, ncp=min(ncp,dim(X)[2]), scale.unit=scale, ind.sup=ind.sup, quanti.sup=quanti.sup, quali.sup=quali.sup, graph=FALSE)
        
        if (verbose) { 
            print("---------------- Sample coordinates on axes -----------------")
            print(acp$ind$coord)
        }
        
        if ((sample.qual) & (verbose)) {
            rep.sample <- qualitySample(acp, axes=c(1:min(ncp, dim(X)[2])))
            print("---------------- Sample Quality -----------------")
            print(rep.sample)
        }
        
        if ((variable.qual) & (verbose)) {
            print("---------------- Variable Quality -----------------")
            print(acp$var$cos2)
        }
        
        if ((sample.cont) & (verbose)) {
            print("---------------- Sample Contribution (%) -----------------")
            print(acp$ind$contrib)
        }
        
        if ((variable.cont) & (verbose)) {
            print("---------------- Variable Contribution (%) -----------------")
            print(acp$var$contrib)
        }
        
        if (plotInertia) {
            if (!pdf) {
                x11()
            }
            plotInertia(acp, ncp=ncp, ...)
        }
        
        if (plotBiplot) {
            plotBiplot(acp, ...)
        }
        
        if (plotSample) {
            plotSample(acp, axes=c(1,2), new.plot=!pdf, lab=lab.sample, lim.cos2.sample=lim.cos2.sample, palette=palette, ...)
            plotSample(acp, axes=c(1,3), new.plot=!pdf, lab=lab.sample, lim.cos2.sample=lim.cos2.sample, palette=palette, ...)
            plotSample(acp, axes=c(2,3), new.plot=!pdf, lab=lab.sample, lim.cos2.sample=lim.cos2.sample, palette=palette, ...)
        }
        
        if (plotVariable) {
            plotVariable(acp, axes=c(1,2), new.plot=!pdf, lab=lab.var, lim.cos2.var=lim.cos2.var, palette=palette, ...)
            plotVariable(acp, axes=c(1,3), new.plot=!pdf, lab=lab.var, lim.cos2.var=lim.cos2.var, palette=palette, ...)
            plotVariable(acp, axes=c(2,3), new.plot=!pdf, lab=lab.var, lim.cos2.var=lim.cos2.var, palette=palette, ...)
        }
        
        if (plot3dSample) {
            if (!pdf) {
                x11()
            }
            plot3dSample(acp, lab=lab.sample, palette=palette, ...)
        }
  
  if (pdf) {
      par(op)
      dev.off()
  }
  
  res <- acp
  return(res)
}
