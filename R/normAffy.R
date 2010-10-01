normAffy <-
function(filenames, celfile.path, method=c("GCRMA","RMA","MAS5"), cdfname = NULL, rmaffx=TRUE){
    ## check if arguments are ok
    if( missing(filenames) & missing(celfile.path) )
        stop("** FAILURE : 'filenames' or 'celfile.path' arguments are empty")
    if(class(celfile.path)!="character")
        stop("** FAILURE : 'celfile.path' argument is not a character class")
   
    celfile.path <- file.path(celfile.path)
    if (missing(filenames)){
        filenames <- list.celfiles(path=celfile.path, full.names=FALSE)
    }
    else{
        filenames<-as.vector(filenames)
    }
    print (filenames)
    method <- match.arg(method)
    
    if(method == "RMA"){
        print("--- RMA normalization ---")
        normData <- justRMA(filenames=filenames, celfile.path=celfile.path, cdfname=cdfname)
        normData <- exprs(normData)
    }
    if (method == "MAS5"){
        rawData <- ReadAffy(filenames=filenames, celfile.path=celfile.path, cdfname=cdfname)
        print("--- MAS5 normalization ---")
        normData <- mas5(rawData)
        normData <- log2(exprs(normData))
    }
    if (method == "GCRMA"){
        require(gcrma)
        print("--- justGCRMA normalization ---")
        normData <- justGCRMA(filenames=filenames, celfile.path=celfile.path, type="affinities", fast=TRUE, cdfname=cdfname)
        normData <- exprs(normData)
    }
    
    if(rmaffx){
        normData=normData[-grep("AFFX", rownames(normData)),]
    }
    
    normData
}

