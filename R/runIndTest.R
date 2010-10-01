runIndTest<-function (data, labels, gene.names = NULL, plot = TRUE, dirname= NULL, grp.name=c("Group1","Group2"))
{
    if (is.null(gene.names)) {
        if (is.null(rownames(data))) {
            gene.names = 1:nrow(data)
        }
        else {
            gene.names = rownames(data)
        }
    }
    if (is.vector(data)) {
        grp1 <- matrix(data[which(labels == 0)], nrow = 1)
        grp2 <- matrix(data[which(labels == 1)], nrow = 1)
        cpt <- 1
    }
    
    else {
        grp1 <- data[, which(labels == 0)]
        grp2 <- data[, which(labels == 1)]
        cpt <- nrow(data)
    }
    out <- as.data.frame(matrix(NA, ncol = 3, nrow = cpt))
   colnames(out) <- c("Probeset", "StatisticalTest", "P-values")
    if (cpt>25){
        plot=FALSE
    }
 
    for (i in 1:cpt) {
        if(sd(round(grp1[i, ]),3)!= 0 && sd(round(grp2[i, ]),3)!= 0){
            testnorm.grp1 <- shapiro.test(grp1[i, ])
            testnorm.grp2 <- shapiro.test(grp2[i, ])
            
            ##loi non normale
            if (testnorm.grp1$p.value < 0.05 && testnorm.grp2$p.value < 0.05) {
                test <- "Wilcoxon"
                testmoyenne <- wilcox.test(grp1[i, ], grp2[i, ])
            }
            else {
                testvariance <- var.test(grp1[i, ], grp2[i, ])
                ##var non egale
                if (testvariance$p.value < 0.05) {
                    test <- "Welch"
                    testmoyenne <- t.test(grp1[i, ], grp2[i, ], var.equal = FALSE)
                }
                ##var egale
                else {
                    test <- "Student"
                    testmoyenne <- t.test(grp1[i, ], grp2[i, ], var.equal = TRUE)
                }
            }
            pval <- round(testmoyenne$p.value, 3)
            if (plot){
                if (!is.null(dirname)){
                    dirnamepng=paste(dirname,"/",rownames(data)[i],".png",sep="")
                    bitmap(file=dirnamepng,type="png16m",taa=4, gaa=4, height = 6, width = 6, res=150)
                }else{x11()}
                
                hg1 <- hist(grp1[i, ], plot = FALSE)
                hg2 <- hist(grp2[i, ], plot = FALSE)
                par(mfrow = c(1, 3))
                ym <- max(c(hg1$counts, hg2$counts), na.rm=TRUE)+1
                plot(hg1, col="light blue", main=grp.name[1], ylim=c(0,ym), xlab="Expression Level")
                points(density(grp1[i, ], na.rm = TRUE), type = "l", lwd = 1, col = "red")
                boxplot(grp1[i, ], grp2[i, ], main = paste(gene.names[i],"\n",test,"-",pval),names = grp.name, col = c("light blue", "#2896C8"))
                plot(hg2, main=grp.name[2], col="#2896C8", xlab="Expression Level", ylim=c(0,ym))
                points(density(grp2[i, ], na.rm = TRUE), type = "l",lwd = 1, col = "red")
                if (!is.null(dirname)){
                    dev.off()
                }
            }
            out[i, ] <- c(gene.names[i], test, pval)
       }else {out[i,]<-c(gene.names[i], "NA", "NA")}
    }
    return(out)
}
