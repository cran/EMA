`test.nested.model` <-
function(X,X0,Y) {
  r <- ncol(X)
  r0 <- ncol(X0)
  n <- nrow(X)
  theta <- solve(t(X)%*%X)%*%t(X)%*%Y
  Y.pred <- X%*%theta
  Y0.pred <- X0%*%solve(t(X0)%*%X0)%*%t(X0)%*%Y

### calcul du numerateur du test de fisher
  NUM <- Y.pred-Y0.pred
  NUM <- NUM^2
  NUM <- colSums(NUM)/(r-r0)

### calcul du denominateur du test de fisher
  DEN <- Y-Y.pred
  DEN <- DEN^2
  DEN <- colSums(DEN)/(n-r)

### calcul du F
  F <- NUM/DEN

### calcul de la p-value
  pvalue <- 1-pf(F,r-r0,n-r)

  return(list(theta=theta, F=F, pvalue=pvalue, residu=Y-Y.pred, sigma2=DEN, X=X))
}

