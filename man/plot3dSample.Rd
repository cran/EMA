\name{plot3dSample}
\alias{plot3dSample}

\title{Sample representation in 3 dimensions for PCA}
\description{
  Sample representation in 3 dimensions for Principal Component Analysis
}
\usage{
plot3dSample(acp, lab = NULL, palette="rainbow", ...)
}

\arguments{
  \item{acp}{result from PCA or do.pca function}
  \item{lab}{colored by lab}
  \item{palette}{character, name of color palette, by default = "rainbow"}
  \item{...}{Arguments to be passed to methods, such as graphical parameters (see 'par').}

}

\value{
  Sample representation in 3 dimensions 
}

\author{EMA group}

\seealso{\code{\link{runPCA}},\code{\link{PCA}}}
\examples{
\dontrun{
data(marty)

## PCA on sample
example.subset <- marty[1:100,1:100]
pca <- runPCA(t(example.subset), plotSample = FALSE, plotInertia = FALSE, verbose = FALSE)

## Sample plot in 3 dimensions colored by tumour type
plot3dSample(pca, lab = marty.type.cl)
}
}

\keyword{plot}

