\name{enrichGO}
\alias{enrichGO}
\title{GO Enrichment Analysis of a gene set.
Given a vector of genes, this function will return the enrichment GO
categories with FDR control.}
\usage{
  enrichGO(gene, organism = "human", ont = "MF",
    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
    readable = FALSE)
}
\arguments{
  \item{gene}{a vector of entrez gene id.}

  \item{organism}{Currently, only "human", "mouse" and
  "yeast" supported.}

  \item{ont}{One of "MF", "BP", and "CC" subontologies.}

  \item{pvalueCutoff}{Cutoff value of pvalue.}

  \item{qvalueCutoff}{Cutoff value of qvalue.}

  \item{readable}{if readable is TRUE, the gene IDs will
  mapping to gene symbols.}
}
\value{
  A \code{enrichResult} instance.
}
\description{
  GO Enrichment Analysis of a gene set. Given a vector of
  genes, this function will return the enrichment GO
  categories with FDR control.
}
\examples{
#data(gcSample)
	#yy <- enrichGO(gcSample[[1]], organism="human", ont="BP", pvalueCutoff=0.01)
	#head(summary(yy))
	#plot(yy)
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}
\seealso{
  \code{\link{enrichResult-class}},
  \code{\link{compareCluster}}
}
\keyword{manip}

