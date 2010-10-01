clustering.kmeans <-
function(data, N=100, iter.max=20, title="Kmeans - Hierarchical Clustering", dist.s="pearson", dist.g="pearsonabs", method="ward") {

  print ("Running kmeans ...")
  km <- kmeans(data, centers=N, iter.max=iter.max)
  data2cluster<-km$centers

  print ("Running clustering ...")
  print (paste(dist.s,method))
  c.sample <- clustering(data2cluster, metric=dist.s, method=method)
  print (paste(dist.g, method))
  c.gene <- clustering(t(data2cluster), metric=dist.g, method=method)
  clustering.plot(tree=c.sample, tree.sup=c.gene, data=data2cluster, title=title)
  
  return(list(c.km=km, c.sample=c.sample, c.kcenters=c.gene))
}

