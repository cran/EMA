\name{runRankprod}
\alias{runRankprod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rank product for small samples size}
\description{
  Rank product method to identify  differentially expressed genes. This
  method is useful for small samples size. 
}
\usage{
runRankprod(data, labels, q = 0.05, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix, a data frame, or an ExpressionSet object. Each row of 'data' (or 'exprs(data)', respectively) must correspond to a gene, and each column to a sample.}
  \item{labels}{A vector of integers corresponding to observation (column) class labels. For 2 classes, the labels must be 0 and 1.}
  \item{q}{A numeric value specifying the pvalue threshold.}
  \item{plot}{A logical value specifying if drawing plots or not.}
}
\value{
  A list of two dataframes, the identification of  up-regulated  and
  down-regulated genes in class 2 compared to class 1, respectively.
  RP/Rsum : The rank product.
  AdjpValue : The adjusted pvalues.
  RawpValue : The raw pvalues.
  FC(class1/class2) : The fold change calculation.
}
\references{
Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, P.(2004) Rank Products: A simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments, FEBS Letter, 57383-92
 }
\author{EMA group}
\seealso{\code{\link{RPadvance}}}
\examples{
## load data
data(marty)

\dontrun{
## filtering data
marty <- expFilter(marty, threshold=3.5, graph=FALSE)
}

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

## run differential analysis on example set
example.subset <- marty[1:100,]
out <- runRankprod(example.subset, labels=marty.type.num, q=0.05, plot=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
