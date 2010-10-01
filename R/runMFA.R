`runMFA` <-
function(Data, group=NULL, ncp=5, name.group=NULL, type=NULL, ind.sup = NULL, num.group.sup = NULL, graph = TRUE, report.file=NULL, report.pdf=NULL){
    
    ## donn<U+00E9>es sous formes de liste
    if (!is.data.frame(Data) & is.list(Data)){
        nbGrp <- length(Data)
        if (nbGrp < 2) stop("Provide at least 2 groups")
        group <- c()
        Data2MFA <- NULL
        nbInd <- nrow(Data[[1]])
        for (i in 1:nbGrp){
            if (nrow(Data[[i]]) != nbInd) stop("All data sets must have the same number of rows")
            Data2MFA <- cbind(Data2MFA, Data[[i]])
            group <- c(group, ncol(Data[[i]]))
        }
        if (is.null(name.group) & !is.null(names(Data))) name.group <- names(Data)
    }
    ## donn<U+00E9>es en data.frame
    if (is.data.frame(Data) | is.matrix(Data)){   
        if (is.null(group)) stop("No groups provided")
        Data2MFA <- Data
        nbGrp <- length(group)
    }
    
    if (is.null(type)) type <- rep("s", nbGrp)
    
    resMFA <- MFA(Data2MFA, group=group, type=type, ind.sup=ind.sup, ncp=ncp, name.group=name.group, num.group.sup=num.group.sup, graph=graph, axes = c(1,2))
    
    if (!is.null(report.file)|| !is.null(report.pdf)){
        MFAreport(resMFA, file.txt=report.file, file.pdf=report.pdf)
    }
    
    return(resMFA)
}

