PLS <-
function(E,F,n=1, scale=TRUE)
  {

    library(MASS)
    n.ind <- nrow(E)
    if(scale)
      {
        E <- scale(E)*sqrt(n.ind/(n.ind-1))
        F <- scale(F)*sqrt(n.ind/(n.ind-1))

        ## récupération des paramètres de centrage et de scaling
        param.scale <- attr(E,"scaled:scale")/sqrt(n.ind/(n.ind-1))
        param.center <- attr(E,"scaled:center")    
        
      }

    if(!scale)
      {
        E <- scale(E,scale=FALSE)
        F <- scale(F,scale=FALSE)

        ## récupération des paramètres de centrage et de scaling
        param.scale <- NULL
        param.center <- attr(E,"scaled:center")            
      }
    
    norme2 <- function(x)
      {
        if(dim(x)[2]!=1)
          {
            stop("x is not a vector")
          }

        return(sum(x*x))
      }

    
    E0 <- E
    F0 <- F
    W <- matrix(0,ncol(E),n)
    T <- matrix(0,nrow(E),n)
    P <- matrix(0,ncol(E),n)
    C <- matrix(0,n,2)

    for(i in 1:n)
      {
        print(paste("composante",i))
#        print(F0)
        tEF <- t(E0) %*% F0

        ### w1 est vecteur propre de tEFtFE associée à la plus grande valeur propre
        w1 <- t(E0)%*%F0
        w1 <- w1/sqrt(norme2(w1))
        W[,i] <- w1


        ### calcul de t1
        t1 <- E0%*%w1
        T[,i] <- t1
        C[i,1] <- (cov(F0,t1)*(n.ind-1)/n.ind)^2
###        C[i,2] <- calc.cov.theo(E0,F0)
        

        ### on effectue deux régressions
        p1 <- t(E0)%*%t1/norme2(t1)
        P[,i] <- p1
        r1 <- t(F0)%*%t1/norme2(t1)

        ### mise a jour des matrices
        E0 <- E0 - t1%*%t(p1)
        F0 <- F0 - t1%*%t(r1)
                
      }

    Wstar <- W%*%ginv(t(P)%*%W)
    return(list(W=W, Wstar=Wstar, T=T, P=P, C=C, E0=E0, F0=F0, param.center=param.center, param.scale=param.scale))
  }

