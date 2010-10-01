\name{genes.selection}
\alias{genes.selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genes selection}
\description{
Selects the most variant genes of a expression dataset. 
}
\usage{
genes.selection(data, thres.IQR, thres.num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A expression matrix, genes on rows and samples on columns}
  \item{thres.IQR}{IQR threshold - Genes are selected based on the IQR threshold}
  \item{thres.num}{Genes number threshold - Selects the N genes with the highest IQR}
}
\details{
The IQR (Inter Quartile Range) is performed on each genes. If the thres.IQR option is chosen, the most variant genes are selected according to an IQR threshold. If the thres.num option is chosen, the genes are ordered according to their IQR level, and the N first genes are selected.
}
\value{
 The name of the selected genes.
}

\author{EMA group}

\seealso{ \code{\link{IQR}},}
\examples{
data(marty)
data.f<-expFilter(marty, graph=FALSE)

##Select the most variant genes
sel<-genes.selection(data.f, thres.num=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

