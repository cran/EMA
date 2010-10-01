clust.dist <-
function(mat, meth.dis="euclidean") { 
  
  MEASURE <- c("euclidean", "manhattan", "pearson", "pearsonabs", "spearman", "spearmanabs", "jaccard", "dice")
  mea <- pmatch(meth.dis, MEASURE)
  
  if (is.na(mea)) stop("Error :Unknown Metric.")
  if (mea==1) DIS <- as.matrix(daisy(t(mat), metric="euclidean"))      
  if (mea==2) DIS <- as.matrix(daisy(t(mat), metric="manhattan"))    
  if (mea==3) DIS <- (1-cor(mat, method="pearson"))/2
  if (mea==4) DIS <- 1-abs(cor(mat, method="pearson"))
  if (mea==5) DIS <- (1-cor(mat, method="spearman"))/2
  if (mea==6) DIS <- 1-abs(cor(mat, method="spearman"))
  if (mea==7) DIS <- jaccard(mat)                                    
  if (mea==8) DIS <- dice(mat)                                        
  attr(DIS,"Metric")<-meth.dis

  return(DIS)
}
