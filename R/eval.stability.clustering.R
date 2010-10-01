eval.stability.clustering<-
  function(X,nb=c(2:4),f=0.8,nsub=10,s0=0.98, list_DIS=c("euclidean","pearson"), list_ALGO=c("average","complete","ward"), pdfname = NULL, verbose = TRUE){

  
  
  Transform.vector.to.list <-function (v) {
    if (is.integer(v) == FALSE)
      stop("Transform.vector.to.list: the elements of the input vector must be integers", call.=FALSE);
    ## relabeling of the classes in order to avoid (that is the labels of the classes will be consecutive integers)
    n.examples <- length(v);
    new.v <- integer(n.examples); ## vector with relabeled elements
    max.class <- max(v);
    v.index <- integer(max.class); ## vector of the class indices for label translation
    label.class <- 0;
    
    for (i in 1:n.examples) {
      old.label <- v[i]; # old label of the ith example
      if (v.index[old.label] == 0) {
        label.class <- label.class + 1;
        v.index[old.label] <- label.class;
      }
      new.v[i] <- v.index[old.label];	
    }
    ## building of the list of clusters
    cl =list();
    for (i in 1:n.examples) {
      if (length(cl) < new.v[i]) 
        cl[[new.v[i]]] <- i
      else  
        cl[[new.v[i]]][length(cl[[new.v[i]]]) + 1] <- i;
    }
    return(cl);
  }
  
  Do.boolean.membership.matrix <-function(cl, dim.M, examplelabels) {
    M <- matrix(integer(dim.M*dim.M), nrow=dim.M);
    colnames(M) <- rownames(M) <- examplelabels;
    singletons <- integer(dim.M);  
    c <- length(cl); # number of clusters 
    for (j in 1:c) {
      n.ex <- length(cl[[j]]);
      if (n.ex == 1)
        singletons[cl[[j]][1]] <- 1
      else {
        for (x1 in 1:(n.ex-1)) {
          for (x2 in (x1+1):n.ex) {
            x <- cl[[j]][x1];
            y <- cl[[j]][x2];
            M[x,y] <- 1;
          }
        }
      }
    }
    for (x1 in 1:(dim.M-1)) 
      for (x2 in (x1+1):dim.M) 
        M[x2,x1] <- M[x1,x2];
    for (x in 1:(dim.M)) 
      M[x,x] <- singletons[x];
    return(M);
  }
  
  
  Compute.Chi.sq<-function (M, s0){
    n <- ncol(M)
    K <- nrow(M)
    x <- numeric(K)
    for (k in 1:K) {
      for (j in 1:n) {
        if (is.na(M[k, j])==FALSE & M[k, j] > s0) x[k] <- x[k] + 1
      }
    }
    theta <- sum(x)/(n * K)
    
    if (theta == 0) {
      p.value = NA
    } else if (theta == 1) {
      p.value = 1
    } else {
      chi.statistic <- sum((x - n * theta)^2)/(n * theta * (1 - theta))
      p.value <- 1 - pchisq(chi.statistic, K - 1)
    }
    return(p.value)
  }
  
  
  teste.stab<-function (sim.matrix, s0 ){
    n.clusterings <- nrow(sim.matrix)
    n.measures <- ncol(sim.matrix)
    ordered.clusterings <- integer(n.clusterings)
    p.value <- numeric(n.clusterings)
    means <- numeric(n.clusterings)
    variance <- numeric(n.clusterings)
    
    means <- mean(as.data.frame(t(sim.matrix)),na.rm=TRUE)
    for (i in 1:n.clusterings) variance[i] <- var(sim.matrix[i,],na.rm=TRUE)
    
    sorted.means <- sort(means, decreasing = TRUE)
    sorted.indices <- order(means, decreasing = TRUE)
    means <- sorted.means
    variance <- variance[sorted.indices]
    ordered.sim.matrix <- sim.matrix[sorted.indices, ]
    
    ordered.clusterings <- sorted.indices 
    p.value[1] <- 1
    for (k in n.clusterings:2)
      p.value[k] <- Compute.Chi.sq(ordered.sim.matrix[1:k,], s0)
    
    d <- data.frame(ordered.clusterings = ordered.clusterings, p.value = p.value, means = means, variance = variance)
    rownames(d) <- 1:n.clusterings
    
    return(d)
  }
  
  t0=Sys.time()
    
  VALID.STAB<-matrix(0,0,4)
  
  for (j in list_DIS){
    for (q in list_ALGO){
      
      ## STABILITY
      ## with  partition   
      n <- ncol(X)
      n.sub.ex <- ceiling(n * f)
      
      for (t in 1:nsub) {
        Jsim.vector <- rep(NA,length(nb))         
        sub1 <- sample(n, n.sub.ex)
        sub2 <- sample(n, n.sub.ex)
        Xsub1 <- X[, sub1]
        colnames(Xsub1) <- sub1
        Xsub2 <- X[, sub2]
        colnames(Xsub2) <- sub2
        
        if (q=="kcentroids"){
          S1 <- clustering(Xsub1,j,q, nb=nb)
          S2 <- clustering(Xsub2,j,q, nb=nb)
        }
        else{
          S1 <- clustering(Xsub1,j,q)
          S2 <- clustering(Xsub2,j,q)
        }
     
        ##Works on 2 partitions 80%
        for (l in nb){
          if (q!="kcentroids") {
            t1 <- cutree(as.hclust(S1), k=l)
            cl1<-Transform.vector.to.list(t1)
            t2 <- cutree(as.hclust(S2), k=l)
            cl2<-Transform.vector.to.list(t2)
            
            M1<-Do.boolean.membership.matrix(cl1, n.sub.ex, sub1)
            M2<-Do.boolean.membership.matrix(cl2, n.sub.ex, sub2)
          }  
          else{
            cl1<-Transform.vector.to.list(S1[[l-1]])   
            M1<-Do.boolean.membership.matrix(cl1, n.sub.ex, sub1)
            cl2<-Transform.vector.to.list(S2[[l-1]])   
            M2<-Do.boolean.membership.matrix(cl2, n.sub.ex, sub2)
          }
          
          sub.common <- intersect(sub1, sub2)
          label.examples <- as.character(sub.common)
          M1 <- M1[label.examples, label.examples]
          M2 <- M2[label.examples, label.examples]
          
          ## calculate Jaccard coefficient
          Jsim.vector[l-1] <-sum(M1 * M2)/(sum(M1 * M1) + sum(M2 * M2) - sum(M1 * M2))
        }    
        
        STAB<-data.frame(dist=rep(j,length(nb)), algo=rep(q,length(nb)), nb, res=rep(t,length(nb)), Jsim.vector=Jsim.vector)
        VALID.STAB<-rbind(VALID.STAB,STAB)
      } # t resampling
    } # q algo
  } # j linkage
  
  t1=Sys.time()
  tdiff=t1-t0
print(tdiff)
  ##================================================
  ##                CHI-2 TEST - Which one is the most stable ?
  ##================================================
  
  NAMES   <-list()
  RESULTS <-list()
  NBSET   <-c()
  
  for (i in nb) {
    VALID.STABsub<-VALID.STAB[which(VALID.STAB$nb==i),]
    
    sim.matrix   <- matrix(VALID.STABsub$Jsim.vector, ncol = nsub, byrow = TRUE )
    names.matrix <- matrix(paste(paste(VALID.STABsub$nom,VALID.STABsub$dist),VALID.STABsub$algo),ncol = nsub, byrow = TRUE )
    
    D<-teste.stab(sim.matrix, s0 ) 
    NAMES   <-c(NAMES,list(names.matrix[,1]))
    RESULTS <-c(RESULTS,list(D))
    NBSET   <-c(NBSET,paste("nb.clust=",i))
  }
  
  for( i in 1:length(RESULTS)) {
    corr.p.value<-p.adjust(RESULTS[[i]]$p.value,method ="holm")
    RESULTS[[i]]<-cbind(RESULTS[[i]], corr.p.value)
  }
  
  RES <-list()
  NAM <-list()
  for( i in 1:length(RESULTS)) {
    res<- RESULTS[[i]]
    nam<- NAMES[[i]]
    res1<-res[which(res$corr.p.value>0.05),]
    RES <-c(RES,list(res1))
    nam.ord<-nam[res1$ordered.clusterings]
    NAM <-c(NAM,list(nam.ord))
  }
  
  stab.methods<-list()
  for (i in 1:(max(nb)-1))  stab.methods<-c(stab.methods, c(nclass=i+1,data.frame(methods=NAMES[[i]][RES[[i]][,1]],  p.value=RES[[i]][,2] ))) 
  if ( verbose == TRUE ) print(stab.methods) 
  
  ##================================================
  ##                GRAPHICS
  ##================================================
  
  spp<-NULL
  for (w in 1:length(NAM)) spp<-c(spp,NAM[[w]])
  PERC<-sort(   table(spp)*100/(max(nb)-1)  ,decreasing=TRUE)
  
  if (!is.null(pdfname)) {
    pdf(paste(pdfname,".pdf",sep=""))
  }
  par(mfrow=c(1,1))
  par(mar=c(12,4,4,2)+0.1)
  
  barplot(PERC, las=2, ylab ="%",ylim=c(0,100), cex.names = 0.75, axes=FALSE)   
  axis(2)
  mtext("Frequency of stable methods")
  
  if (!is.null(pdfname)) {
    dev.off()
  }
  

  t2=Sys.time()
  tdiff=t2-t1
  print(tdiff)
  
  stab.methods
}


