inverse <- function(X)
{
  res.svd <- svd(X)
  U <- res.svd$u
  V <- res.svd$v
  p <- length(res.svd$d)
  D <- matrix(0,p,p)
  diag(D) <- 1/res.svd$d
  return(V%*%D%*%t(U))
}


`test.LC` <-
function(C, X, Y, global=FALSE, cor.multtest=TRUE,typeFDR="FDR-BH") {
### test sur plusieurs combinaisons lineaires
## C:      vecteur ou matrice des combinaisons lineaires a tester (matrice: une ligne = une CL)
## X:      matrice de design, obtenue a partir du modele lineaire
## Y:      vecteur ou matrice de reponse (matrice: une colonne = un gene)
## global: booléen, si TRUE test de F global de l'hypothèse H_0 = 0 pour l'ensemble des contrastes donnés, si FALSE test des k contrastes un à un, H_0i = 0 pour i dans [1..k] 
## cor.multtest: booléen, si global=TRUE et cor.multtest=TRUE une correction de test multiple est appliquée pour le calcul des pvalues en appelant la fonction multiple.correction()
## typeFDR: argument typeFDR pour multiple.correction(). 
    
    if(!is.logical(global))
      stop("Error : logical value expected for option 'global' ")
    if(!is.logical(cor.multtest))
      stop("Error : logical value expected for option 'cor.multtest' ")
    r <- ncol(X)
    q <- nrow(C)   
    n <- nrow(X)
    if (class(C)!="matrix") C=t(as.matrix(C))

    invX <- solve(t(X) %*% X)
    theta <- invX %*% t(X) %*% Y
    Y.pred <- X %*% theta
    DEN <- (Y - Y.pred)^2
    sigma2 <- colSums(DEN)/(n - r)
    phi <- C %*% theta

    if(global){
      pvalue=vector(mode="numeric",length=length(sigma2))
      if(is.matrix(Y) & !is.null(colnames(Y))) names(pvalue)=colnames(Y)
      for(i in 1:length(sigma2)){
            var.theta <- invX*sigma2[i]
            Fvalue <- t(phi[,i])%*%inverse(C%*%var.theta%*%t(C))%*%phi[,i]/q
            pvalue[i] <- 1-pf(Fvalue,q,n-r)
      }
    }

    if(!global){
      Fvalue <- phi^2/(as.matrix(rowSums((C %*% invX) * C)) %*% as.matrix(t(sigma2)))
      pvalue <- 1 - pf(Fvalue, 1, n - r)
      # Here there are 2 possible configurations : 
      # 1. Y is a vector, p-value is a matrix with 1 column and p lines (one per contrast), multiple correction applies to the p different contrasts. This is the situation where we test many (more than 10) linear combinations of parameters of a linear model.
      # 2. Y is a matrix (gene expression typical setting), p-value is a matrix with n column (one per gene) and p lines (one per contrast), multiple correction applies to the n genes. This is the situation where we test a few combinations of parameters of a linear model on many (more than 30) variables or genes.
      # We treat the 2 situations the same way, that's why if we are in the first one we first transform the one column matrix in a one row matrix
      if(ncol(pvalue)==1) pvalue=t(pvalue)
      if(cor.multtest==TRUE){
            for(i in 1:nrow(pvalue)){
                  pvalue[i,]=multiple.correction(pvalue[i,],typeFDR)
            }
      }
    }

    return(list(Estimate = phi, Fvalue = Fvalue, pvalue = pvalue, Y.pred = Y.pred, resid = Y-Y.pred, sigma2 = sigma2, theta = theta))
}



