\name{groupGO}
\alias{groupGO}
\title{Functional Profile of a gene set at specific GO level.
Given a vector of genes, this function will return the GO profile at
specific level.}
\usage{
  groupGO(gene, organism = "human", ont = "CC", level = 2,
    readable = FALSE)
}
\arguments{
  \item{gene}{a vector of entrez gene id.}

  \item{organism}{Currently, only "human" and "mouse"
  supported.}

  \item{ont}{One of "MF", "BP", and "CC" subontologies.}

  \item{level}{Specific GO Level.}

  \item{readable}{if readable is TRUE, the gene IDs will
  mapping to gene symbols.}
}
\value{
  A \code{groupGOResult} instance.
}
\description{
  Functional Profile of a gene set at specific GO level.
  Given a vector of genes, this function will return the GO
  profile at specific level.
}
\examples{
data(gcSample)
	yy <- groupGO(gcSample[[1]], organism="human", ont="BP", level=2)
	head(summary(yy))
	#plot(yy)
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\seealso{
  \code{\link{groupGOResult-class}},
  \code{\link{compareCluster}}
}
\keyword{manip}

