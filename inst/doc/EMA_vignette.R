### R code from vignette source 'EMA_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: B0 (eval = FALSE)
###################################################
## ## ILLUMINA CHIPS
## ## From lumi1.14.0 help page
## ## load example data
## require(lumi)
## data(example.lumi)
## ## Do all the default preprocessing in one step
## lumi.N <- lumiExpresso(example.lumi)
## ## Convert data into matrix for EMA
## lumi.mat <- exprs(lumi.N)


###################################################
### code chunk number 2: B1 (eval = FALSE)
###################################################
## ## cDNA CHIPS
## ## From vsn3.16.0 help page
## ## VSN normalisation
## require(vsn)
## data(lymphoma)
## lym <- justvsn(lymphoma)
## ## Convert data into matrix for EMA
## lym.mat <- exprs(lym)


###################################################
### code chunk number 3: B2
###################################################
##Graphical parameters
par(cex.axis=0.6, cex.main=0.7, cex.lab=0.7)
##Load package
require(EMA)


###################################################
### code chunk number 4: B3 (eval = FALSE)
###################################################
## ##Load EMA package
## require(EMA)
## 
## ## GCRMA Normalisation
## ## Not run because cel files are not available from this package
## cel.path=paste(getwd(),"/Data/E-GEOD-13787", sep="")
## marty<-normAffy(celfile.path=cel.path, method="GCRMA")


###################################################
### code chunk number 5: C
###################################################
##Load GCRMA Normalised data
##marty.type.cl=Her2+ corresponds to Her2+ breast cancer
##marty.type.cl=Basal corresponds to Basal-Like carcinoma
data(marty)

##And discard probesets with a maximum log2 expression value below 3.5
marty.f<-expFilter(marty)
dim(marty.f)


###################################################
### code chunk number 6: D
###################################################
acp<-runPCA(t(marty.f), scale=FALSE,lab.sample=marty.type.cl, 
            plotSample=FALSE, plotInertia=FALSE)
plotInertia(acp)


###################################################
### code chunk number 7: D1
###################################################
## Individual map (axe 1 and 2)
plotSample(acp,axes=c(1,2),lab=marty.type.cl)


###################################################
### code chunk number 8: D2 (eval = FALSE)
###################################################
## ## Or create a pdf report with selected plots
## runPCA(t(marty.f), scale=FALSE, pdfname="PCA.pdf",lab.sample=marty.type.cl)


###################################################
### code chunk number 9: E
###################################################
## PCA after normalisation and without filtering, but with scaling
## done on 1000 most variant variables
mvgenes<-genes.selection(marty.f, thres.num=1000)
acp<-runPCA(t(marty[mvgenes,]), scale=TRUE, lab.sample=marty.type.cl, 
            plotSample=FALSE, plotInertia=FALSE)

## Individual map (axe 1 and 2)
plotSample(acp,axes=c(1,2),lab=marty.type.cl)


###################################################
### code chunk number 10: F
###################################################
##Gene representation (only genes the most correlated to the two first components)
level.cl<-ifelse(apply(marty,1,max)>3.5,"High","Low")
plotVariable(acp,axes=c(1,2),lim.cos2.var=0.8,lab=level.cl,label="none")


###################################################
### code chunk number 11: G
###################################################
## Sample Hierarchical Clustering (Pearson's correlation coefficient and Ward method)
c.sample<-clustering(data=marty.f, metric="pearson", method="ward")
clustering.plot(tree=c.sample,  lab=marty.type.cl, 
                title="GCRMA Data - filtered")


###################################################
### code chunk number 12: G1
###################################################
## Heatmap performed on the 100 probesets with the highest IQR values
mvgenes<-genes.selection(marty.f, thres.num=100)
c.sample<-clustering(data=marty.f[mvgenes,], metric="pearson", method="ward")
c.gene<-clustering(data=t(marty.f[mvgenes,]), metric="pearsonabs", method="ward")
clustering.plot(tree=c.sample, tree.sup=c.gene, data=marty.f[mvgenes,],
                names.sup=FALSE, lab= marty.type.cl, trim.heatmap=0.99)


###################################################
### code chunk number 13: H (eval = FALSE)
###################################################
## 
## ### MFA
## 
## ### BiClustering - BiCare package


###################################################
### code chunk number 14: I
###################################################
### Student test with BH correction
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)
rt<-runTtest(marty.f, labels=marty.type.num, algo="t.equalvar", q=0.05, plot=FALSE)
head(rt)


###################################################
### code chunk number 15: J (eval = FALSE)
###################################################
## rs<-runSAM(marty.f, labels=marty.type.num)
## head(rs)


###################################################
### code chunk number 16: J4 (eval = FALSE)
###################################################
## ##Look at the annotation of the DEG
## rt.sign <- rt[order(rt$AdjpValue),]
## rt.annot<-bioMartAnnot(rt.sign, inputTypeId ="affy_hg_u133_plus_2",
##                        outputTypeId=c("entrezgene","hgnc_symbol"), 
##                        dataset=c("hsapiens_gene_ensembl"),database = "ensembl")


###################################################
### code chunk number 17: J5 (eval = FALSE)
###################################################
## ## Not run because cel files are not available from this package
## filenames <- list.files("Data/E-GEOD-13787", pattern=".CEL", ignore.case=TRUE)
## rawdata <- ReadAffy(filenames=filenames, celfile.path="Data/E-GEOD-13787", cdfname=NULL)
## probePlots(rawdata, pbsList=rt.annot$affy_hg_u133_plus2[1:10])


###################################################
### code chunk number 18: J3 (eval = FALSE)
###################################################
## ## http://www.broad.mit.edu/gsea/msigdb/msigdb_index.html
## ## You have to register first and then download the gmt file from their site
## gsaOUT <- runGSA(marty, marty.type.num , 
## 	gmtfile="c2.kegg.v2.5.symbols.gmt", 
## 	chip="hgu133plus2")


###################################################
### code chunk number 19: K (eval = FALSE)
###################################################
## ## GO and KEGG analysis on the DEG by the SAM procedure
## runHyperGO(list=rownames(rt.sign), pack.annot="hgu133plus2.db", name="HyperGO_type")
## runHyperKEGG(list=rownames(rt.sign), pack.annot="hgu133plus2.db", name="HyperKEGG_type")


###################################################
### code chunk number 20: L
###################################################
set.seed(5000)
gene<-rnorm(100)
gene[51:100]<-gene[51:100]+2
group<-ifelse(gene<=median(gene),"Low gene expression","High gene expression")
time<-abs(rnorm(100))
time[51:100]<-time[51:100]+2
status<-sample(c(0,1),size=100,replace=TRUE)

res<-km(time,status,group,main="Kaplan Meier curve")
res$fit.km
res$lr
res$p.lr


###################################################
### code chunk number 21: sessionInfo
###################################################
toLatex(sessionInfo())


