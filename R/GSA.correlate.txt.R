GSA.correlate.txt <-
function (GSA.genesets.obj, genenames) {
  
  nsets <- length(GSA.genesets.obj$genesets)
  ngenes <- unlist(lapply(GSA.genesets.obj$genesets, length))
  allgenes <- unlist(GSA.genesets.obj$genesets)
  sets.in.exp <- match(unique(allgenes), genenames)
  exp.in.sets <- match(genenames, allgenes)
  
  print(paste("Number of gene-sets :", nsets, sep=" "))
  
  print(paste("Total number of unique genes in gene-set collection :", length(unique(allgenes)),sep=" "))
    
  print(paste("Total number of unique genes in genenames list :", length(unique(genenames)),sep=" "))
    
  print(paste("Number of unique genes in both collections :", sum(!is.na(sets.in.exp)),sep=" "))
  nn <- rep(NA, nsets)
  for (i in 1:nsets) {
    nn[i] <- sum(!is.na(match(GSA.genesets.obj$genesets[[i]],genenames)))
  }
  print("Quantiles of fraction coverage of gene-sets")
  print(quantile(nn/ngenes, seq(0, 1, by = 0.1)))
  
  return()
}

