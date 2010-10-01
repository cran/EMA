runGSA <-
function(nData, labels, gmtfile, chip="hgu133plus2" ,np=1000, minsize=10, maxsize=800, resp.type="Two class unpaired", fdr=0.25){

  #pkg <- chip
  pkg <- paste(chip, "db", sep=".")
  require(pkg, character.only=TRUE)

  if(resp.type != "Two class paired" && is.element(c(0,1), unique(labels))[1] == TRUE && is.element(c(0,1), unique(labels))[2] == TRUE)
    labels <- labels+1
  
  geneset.obj <- GSA.read.gmt(gmtfile)
  genenames <- mget(rownames(nData), eval(as.name(paste(chip,"SYMBOL", sep=""))))
  genenames <- unlist(genenames)
  minsize<-as.numeric(minsize)
  maxsize<-as.numeric(maxsize)

  print("GSA Analysis is running...")
  GSA.obj <- GSA(nData,labels, genenames=genenames, genesets=geneset.obj$genesets, resp.type=resp.type, method="maxmean", nperms=np, minsize=minsize ,maxsize=maxsize)

  print("GSA parameters :")
  print(paste("Number of permutations :", np, sep=" "))
  print(paste("Type of analysis :", resp.type, sep=" "))
  print(paste("Minimun number of genes in a gene set :", minsize, sep=" "))
  print(paste("Maximun number of genes in a gene set :", maxsize, sep=" "))
  print(paste("FDR chosen :", fdr, sep=" "))
  GSA.correlate.txt(geneset.obj, genenames)

  results <- GSA.listsets(GSA.obj,geneset.names=geneset.obj$geneset.names, FDRcut=fdr)

  ## genes scores for negative set
  gscore.negative=list()
  if(dim(results$negative)[[1]]!=0){
      for (i in 1:dim(results$negative)[[1]]){
          gscore.negative[[results$negative[,"Gene_set_name"][i]]]=GSA.genescores(as.numeric(results$negative[,"Gene_set"][i]),geneset.obj$genesets,GSA.obj,as.vector(genenames),negfirst=TRUE)
      }
  }
  ## genes scores for positive sets
  gscore.positive=list()
  if(dim(results$positive)[[1]]!=0){
      for (i in 1:dim(results$positive)[[1]]){
          gscore.positive[[results$positive[,"Gene_set_name"][i]]]=GSA.genescores(as.numeric(results$positive[,"Gene_set"][i]),geneset.obj$genesets,GSA.obj,as.vector(genenames))
      }
  }

  results$genescore.positive=gscore.positive
  results$genescore.negative=gscore.negative

  par(mfrow=c(nr=2, nc=2))
  hist(GSA.obj$pvalues.lo, breaks=50)
  hist(GSA.obj$pvalues.hi, breaks=50)
  hist(GSA.obj$fdr.lo, breaks=50)
  hist(GSA.obj$fdr.hi, breaks=50)

  print("GSA Analysis OK !")
  return(results)
}

