bioMartAnnot <-  
function (data, inputTypeId, outputTypeId = c("entrezgene", "hgnc_symbol", "ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position", "band", "strand"), dataset= c("hsapiens_gene_ensembl"), database = "ensembl",  sort.by = NULL, outfile = NA) {

    if (missing(data) || missing(inputTypeId)){
      stop("Error : Data or inputTypeId not found")  
    }
    if (missing(database) || missing(dataset)){
      stop("Error : You have to define the database and the dataset you want to use")  
    }
    else{
      mart <- useMart(database, dataset = dataset)
    }
    
    ##Input - vector or matrix
    if (is.vector(data)) {	
      namesGenes <- data
      data <- data.frame(row.names = namesGenes)
      data[[inputTypeId]] <- namesGenes
    }
    else {
      namesGenes <- rownames(data)
      data[[inputTypeId]] <- rownames(data)
    }

    d<-getBM(attributes = unique(c(inputTypeId, outputTypeId, if (!("hgnc_symbol" %in% c(outputTypeId, inputTypeId))) "hgnc_symbol", if (!( "ensembl_gene_id" %in% c(outputTypeId, inputTypeId)))  "ensembl_gene_id")), filters = inputTypeId, values = namesGenes, mart = mart) 
  
    ## Unannotated Gene
    noAnnot <- namesGenes[!(namesGenes %in% d[[1]])]
    
    if (length(noAnnot) > 0) {
      nr <- nrow(d)
      d[(nr+1):(nr+length(noAnnot)),inputTypeId] = noAnnot
    }

    d <- merge(data, d, by = inputTypeId)

    ## Sort data.frame
    if (!is.null(sort.by)) {
      d <- d[order(d[[sort.by]], decreasing = TRUE),]
    }
    else{
      d <- d[unlist(sapply(namesGenes, function(x,col){which(col==x)},col=d[,inputTypeId])),]
    }

    if (!is.na(outfile)) {
      dhtml = d
      ##Add 2 columns :  link_geneCards = gene cards web address for the gene
      ##link_proteinAtlas = protein atlas web address for the gene
      dhtml$link_geneCards = paste("<p><a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=", dhtml$hgnc_symbol, "\">", dhtml$hgnc_symbol,"</a></p>",sep = "")
      dhtml$link_proteinAtlas = paste("<p><a href='http://www.proteinatlas.org/gene_info.php?ensembl_gene_id=", dhtml$ensembl_gene_id, "'>", dhtml$ensembl_gene_id,"</a>/</p/>",sep = "") 
      dhtml = dhtml[c(inputTypeId, setdiff(names(data),c(inputTypeId,outputTypeId)),"link_geneCards", "link_proteinAtlas", if (!is.null(outputTypeId)) outputTypeId[!(outputTypeId %in% inputTypeId)])]
      d = d[c(inputTypeId, setdiff(names(data),c(inputTypeId,outputTypeId)), if (!is.null(outputTypeId)) outputTypeId[!(outputTypeId %in% inputTypeId)])]
      
      write.table(d, file = paste(outfile,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE) 
      print(xtable(dhtml), type = "html", file = paste(outfile,sep=""), sanitize.text.function = force,include.rownames=FALSE)
    }
    return(d)
}

