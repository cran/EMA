`MFAreport` <-
function(resMFA, file.txt=NULL, file.pdf=NULL){

    if (!is.null(file.txt)){
        cat("GROUPS\n", file=file.txt)
        for (i in 1:length(resMFA$group)){
            cat(paste("\n", names(resMFA$group)[i], "\n"), file=file.txt, append=TRUE)
            write.infile(resMFA$group[[i]], file=file.txt, append=TRUE, sep="\t")
        }
        cat("\n\n\n", file=file.txt, append=TRUE)
        cat(rep("#", 25), file=file.txt, append=TRUE)
        cat("\nPartial axes\n", file=file.txt, append=TRUE)
        for (i in 1:(length(resMFA$partial.axes)-1)){
            cat(paste("\n", names(resMFA$partial.axes)[i], "\n"), file=file.txt, append=TRUE)
            write.infile(resMFA$partial.axes[[i]], file=file.txt, append=TRUE, sep="\t")
        }
        cat("\n cor entre les facteurs partiels\n", file=file.txt, append=TRUE)
        write.infile(resMFA$partial.axes[[4]], file=file.txt, append=TRUE, sep="\t")
        cat(rep("#", 25), file=file.txt, append=TRUE)
        cat("\nIndividuals\n", file=file.txt, append=TRUE)
        for (i in 1:length(resMFA$ind)){
            cat(paste("\n", names(resMFA$ind[i]), "\n"), file=file.txt, append=TRUE)
            write.infile(resMFA$ind[[i]], file=file.txt, append=TRUE, sep="\t")
        }
    }

    if (!is.null(file.pdf)){
        pdf(file=file.pdf, height=8.26, width=11.69)
        plot(resMFA, axes=c(1,2), choix="ind", new.plot=FALSE, title="global")
        plot(resMFA, axes=c(3,4), choix="ind", new.plot=FALSE, title="global")
        plot(resMFA, axes=c(1,2), choix="ind", new.plot=FALSE, title="global with partial ind", partial="all")
        plot(resMFA, axes=c(3,4), choix="ind", new.plot=FALSE, title="global with partial ind", partial="all")
        for (i in 1:length(resMFA$separate.analyses)){
            tt <- try(plot(resMFA$separate.analyses[[i]], choix="ind", title=paste("partial -", names(resMFA$separate.analyses)[i]), new.plot=FALSE, axes=c(1,2)))
            tt <- try(plot(resMFA$separate.analyses[[i]], choix="ind", title=paste("partial -", names(resMFA$separate.analyses)[i]), new.plot=FALSE, axes=c(3,4)))
        }
        plot(resMFA, choix="group", new.plot=FALSE, axes=c(1,2))
        plot(resMFA, choix="group", new.plot=FALSE, axes=c(3,4))
        dev.off()
    }
}

