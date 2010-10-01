expFilter <-
function(data, threshold=3.5, graph=TRUE){
  vecMax <- apply(data, 1, max)
  if (graph){
    hist(data, breaks=200, main="Data distribution", xlab="Expression Level")
    abline(v=threshold, col=2)
    mtext(threshold, at=threshold, side=1, las=1, col=2, cex=0.6)

  }
  data <- data[which(vecMax>threshold),]
}

