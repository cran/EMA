`keggReport` <-
function(hgOver, fileout = "report.txt", pack.annot, pvalue = 0.05) {
  pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep="")))
  pack.annot.SYMBOL <- eval(as.name(paste(gsub(".db", "", pack.annot), "SYMBOL", sep="")))
  
  write(file = fileout, paste("KEGG HyperGTest report\n"), append = FALSE, sep = ",")
  entrez <- unlist(as.list(pack.annot.EID))

  a <- geneIdsByCategory(hgOver)
  b <- geneIdUniverse(hgOver, cond=conditional(hgOver))
  
  a <- a[sigCategories(hgOver, pvalue)]
  b <- b[sigCategories(hgOver, pvalue)]
  
  for (i in as.vector(unlist(attributes(a)))) {
    if (length(b[[i]])>10) {
      write(file = fileout, paste("\n",i,":", length(a[[i]]), " geneIds | ", length(b[[i]]), " universIds"), append = TRUE, sep = ",", ncol = length(a[[i]]))
      
      pbset <- unique(names(entrez[which(is.element(entrez, a[[i]]))]))
      pbset <- intersect(pbset, unique(names(entrez[is.element(entrez,geneIds(hgOver))])))
      write(file = fileout, paste(length(pbset),"ProbeSets :"), append = TRUE, sep = ",")
      write(file = fileout, pbset, append = TRUE, sep=",", ncol = length(pbset))
      
      gs <- unique(unlist(mget(pbset, pack.annot.SYMBOL)))
      write(file = fileout, paste(length(gs), "GeneSymbols :"), append=TRUE, sep=",")
      write(file = fileout, gs, append = TRUE, sep=",", ncol = length(gs))     
    }
  }
}

