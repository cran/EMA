`goReport` <-
function(hgOver, fileout = "report.txt", type = c("CC", "MF", "BP"), pack.annot, pvalue = 0.05,  categorySize = 1) {

  pack.annot.ALLPROBES <- eval(as.name(paste(gsub(".db", "", pack.annot), "GO2ALLPROBES", sep = "")))
  pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep = "")))
  pack.annot.SYMBOL <- eval(as.name(paste(gsub(".db", "", pack.annot), "SYMBOL", sep = "")))

  if (type == "CC") {
    write(file = fileout, paste("Component Cellular - GO HyperGTest report\n"), append = FALSE, sep = ",")
    probes <- mget('GO:0005575', pack.annot.ALLPROBES)
  } else if (type == "BP") {
    write(file = fileout, paste("Biological Process - GO HyperGTest report\n"), append = FALSE, sep = ",")
    probes <- mget('GO:0008150', pack.annot.ALLPROBES)
  } else if (type == "MF") {
    write(file = fileout, paste("Molecular Function - GO HyperGTest report\n"), append = FALSE, sep = ",")
    probes <- mget('GO:0003674', pack.annot.ALLPROBES)
  }
  entrez <- mget(unique(unlist(probes)), pack.annot.EID)

  a <- geneIdsByCategory(hgOver)
  b <- geneIdUniverse(hgOver, cond=conditional(hgOver))
  
  a <- a[sigCategories(hgOver, pvalue)]
  b <- b[sigCategories(hgOver, pvalue)]
  
  for (i in as.vector(unlist(attributes(a)))) {
    if (length(b[[i]]) > categorySize) {
      write(file = fileout, paste("\n",i,":", length(a[[i]]), " geneIds | ", length(b[[i]]), " universIds"), append = TRUE, sep=",", ncol = length(a[[i]]))
      
      pbset <- unique(names(unlist(entrez[which(entrez%in%a[[i]])])))
      pbset <- intersect(pbset, unique(names(entrez[is.element(entrez,geneIds(hgOver))])))
      write(file = fileout, paste(length(pbset),"ProbeSets :"), append = TRUE, sep = ",")
      write(file = fileout, pbset, append = TRUE, sep = ",",ncol = length(pbset))
      
      gs <- unique(unlist(mget(pbset, pack.annot.SYMBOL)))
      write(file = fileout, paste(length(gs), "GeneSymbols :"), append = TRUE, sep = ",")
      write(file = fileout, gs, append = TRUE, sep=",", ncol = length(gs))     
    }
  }
}
