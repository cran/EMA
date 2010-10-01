plotBiplot <-
function(acp, ...) {

## Goal
## Sample and variable representation on a same graph

## Input
## acp : result from PCA function

## Output
## Plot of samples and variables on a same graph

  v<-acp$eig[1:2,2]
  cp1 <- round(v[1], digit = 2)
  cp2 <- round(v[2], digit = 2)
  lab.x <- paste("Dimension ",1," (",cp1,"%)",sep="")
  lab.y <- paste("Dimension ",2," (",cp2,"%)",sep="")
  nind<-dim(acp$ind$coord)[1]
 
  indcoord<-t(t(acp$ind$coord[,1:2])/sqrt(acp$eig)[1:2,1])/sqrt(nind)
  x <- indcoord
  y <- acp$var$coord[,1:2]*(sqrt(nind))
  biplot(x,y, ...)

  unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
  rangx1 <- unsigned.range(x[, 1])
  rangx2 <- unsigned.range(x[, 2])
  rangx1 <- rangx2 <- range(rangx1, rangx2)
  rangy1 <- unsigned.range(y[, 1])
  rangy2 <- unsigned.range(y[, 2])
  ratio <- max(rangy1/rangx1, rangy2/rangx2)
  points(indcoord*ratio,pch=3)
  abline(h=0,lty=2)
  abline(v=0,lty=2)
}

