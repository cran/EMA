`runHyperKEGG` <-
function(list, pack.annot, categorySize = 1, name = "hyperKEGG", htmlreport = TRUE, txtreport = TRUE, tabResult = FALSE, pvalue = 0.05) {

  require(pack.annot, character.only=TRUE)
  
  pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep="")))
  pack.annot.ACC <- eval(as.name(paste(gsub(".db", "", pack.annot), "ACCNUM", sep="")))
  pack.annot.KEGG <- gsub(".db", "",  pack.annot)

  listALL <- list

  filehtmlKEGG <- paste(name, ".KEGG.html", sep="")
  filetxtKEGG <- paste(name, ".KEGG.txt", sep="") 

  ## ENTREZ
  designALL <- names(unlist(as.list(pack.annot.ACC)))
  entrezIds <- mget(designALL, envir = pack.annot.EID, ifnotfound=NA)
  haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
  designKEGG<-haveEntrezId
  
  ## UNIQUE UNIVERSE
  entrezUniverse <- unique(unlist(mget(designKEGG, pack.annot.EID, ifnotfound=NA)))
  
  listKEGG <- intersect(listALL, designKEGG)
  selectedEntrezIds <- unlist(mget(listKEGG, pack.annot.EID, ifnotfound=NA))

  ## UNIQUE LISTE
  selectedEntrezIds <- na.omit(unique(unlist(selectedEntrezIds)))
  params.KEGG.over <-new("KEGGHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse, annotation = pack.annot.KEGG, pvalueCutoff=0.05, testDirection="over")
  
  hgOver.KEGG <- hyperGTest(params.KEGG.over)

  if (htmlreport){
    htmlheader(paste(date(), "<br>KEGG HyperGeometric Analysis<br>", "<br>List length:", length(listALL),"<br><br>"), filename = filehtmlKEGG)
    htmlresult(hgOver.KEGG, filename = filehtmlKEGG, app = TRUE, categorySize, pvalue = pvalue)
  }

  if (txtreport)
    keggReport(hgOver.KEGG, fileout = filetxtKEGG, pack.annot = pack.annot, pvalue = pvalue)

  if (tabResult) {
    res <- summary(hgOver.KEGG, pvalue = pvalue)
    return(res)
  }
  
}

