clustering <-
  function(data, metric="euclidean", method="ward", nb) {

  METHOD <- c("average", "single", "complete", "ward", "weighted","diana","kcentroids")
  mea <- pmatch(method, METHOD)
  
  if (is.na(mea)) stop("Error : Unknown Linkage.")
  if (!is.na(mea) && mea != 7)  DIS <- clust.dist(data, metric)

  if (mea<=5) HIERA <- agnes(DIS, method=method, diss=TRUE, keep.diss=TRUE)
  else if (mea==6) HIERA <- diana(DIS,diss=TRUE)
  else if (mea==7) {
    if (missing(nb)) stop("Number of clusters 'k' have to be selected for kcentroids algorithm")
    HIERA<-list()
    if (metric=="euclidean") {
      for (i in nb) HIERA <- c(HIERA, km=kmeans(t(data), centers=i)$cluster )
    }
    else{
      DIS <- clust.dist(data, metric)
      for (i in nb) HIERA <- c(HIERA, pm=pam(DIS$DIS, k=i, diss=TRUE)$cluster )
    }
  }
 
  attr(HIERA$diss,"Metric")<-attr(DIS,"Metric")
  
  return(HIERA)             
}
