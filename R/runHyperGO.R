`runHyperGO` <-
function(list, pack.annot, categorySize = 1, verbose = TRUE, name = "hyperGO", htmlreport = TRUE, txtreport = TRUE, tabResult = FALSE, pvalue = 0.05) {
  require(pack.annot, character.only=TRUE)

  filehtmlGO <- paste(name, ".GO.html", sep="")
  filetxtCC <- paste(name, ".GO.CC.txt", sep="")
  filetxtBP <- paste(name, ".GO.BP.txt", sep="")
  filetxtMF <- paste(name, ".GO.MF.txt", sep="")

  pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep = "")))
  pack.annot.ACC <- eval(as.name(paste(gsub(".db", "", pack.annot), "ACCNUM", sep = "")))
  pack.annot.GO <- eval(as.name(paste(gsub(".db", "", pack.annot), "GO", sep = "")))
  
  listALL <- list    

  ## ENTREZ
  designALL <- names(unlist(as.list(pack.annot.ACC)))
  entrezIds <- mget(designALL, envir = pack.annot.EID, ifnotfound=NA)
  haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
  print ("ENTREZ ID done.")
  
  ## GO
  haveGo <- sapply(mget(haveEntrezId, pack.annot.GO, ifnotfound=NA),
                   function(x) {
                     if (length(x) == 1 && is.na(x)) {
                       FALSE
                     }else {
                       TRUE
                     }
                   })
  numNoGO <- sum(!haveGo)
  if (verbose)
   print(paste(numNoGO, "pbsets in Universe have no GO ids"))
  
  designGO <- haveEntrezId[haveGo]

  ## UNIQUE UNIVERSE
  entrezUniverse <- na.omit(unique(unlist(mget(designALL, pack.annot.EID, ifnotfound=NA))))
  if (any(duplicated(entrezUniverse))){
    stop("error in gene universe: can't have duplicate Ent")
  }

  ## SELECTED GENE LIST
  listGO <- intersect(listALL, designGO)
  numNoGO <- length(listALL) - length(listGO)
  if (verbose)
   print (paste(numNoGO, " pbsets in the list have no GO ids"))
 
  selectedEntrezIds <- na.omit(unique(unlist(mget(listGO, pack.annot.EID, ifnotfound=NA))))

  ## GO ANALYSIS

  htmlheader(paste(date(), "<br>Gene Ontology HyperGeometric Analysis<br>","<br>List length:", length(listALL),"<br><br>"), filename = filehtmlGO)
  
  params.BP.over <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = entrezUniverse, annotation = pack.annot, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
  params.MF.over <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = entrezUniverse, annotation = pack.annot, ontology = "MF", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
  params.CC.over <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = entrezUniverse, annotation = pack.annot, ontology = "CC", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
  
  hgOver.BP <- hyperGTest(params.BP.over)
  hgOver.MF <- hyperGTest(params.MF.over)
  hgOver.CC <- hyperGTest(params.CC.over)

  if (htmlreport){
    htmlresult(hgOver.BP, filename = filehtmlGO, app = TRUE, categorySize = categorySize, pvalue = pvalue)
    htmlresult(hgOver.MF, filename = filehtmlGO, app = TRUE, categorySize = categorySize, pvalue = pvalue)
    htmlresult(hgOver.CC, filename = filehtmlGO, app = TRUE, categorySize = categorySize, pvalue = pvalue)
  }

  if (txtreport){
    goReport(hgOver.BP, fileout = filetxtBP, type = "BP", pack.annot = pack.annot, categorySize = categorySize, pvalue = pvalue)
    goReport(hgOver.MF, fileout = filetxtMF, type = "MF", pack.annot = pack.annot, categorySize = categorySize, pvalue = pvalue)
    goReport(hgOver.CC, fileout = filetxtCC, type = "CC", pack.annot = pack.annot, categorySize = categorySize, pvalue = pvalue)
  }

  if (tabResult){
    res1 <- summary(hgOver.BP, pvalue = pvalue, categorySize = categorySize)
    res2 <- summary(hgOver.MF, pvalue = pvalue, categorySize = categorySize)
    res3 <- summary(hgOver.CC, pvalue = pvalue, categorySize = categorySize)
    res <- list(BP = res1, MF = res2, CC = res3)
    return(res)
  }else{
      return(list(BP = hgOver.BP, MF = hgOver.MF, CC = hgOver.CC))
  }
}

