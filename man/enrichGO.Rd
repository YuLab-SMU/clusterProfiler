\name{enrichGO}
\alias{enrichGO}
\title{GO Enrichment Analysis of a gene set.}
\description{
  Given a vector of genes, this function will return the enrichment GO categories with FDR control.
}
\usage{
	enrichGO(gene, organism="human", ont="MF", pvalueCutoff = 0.01, testDirection="over")
}

\arguments{
  \item{gene}{ a vector of entrez gene id. }
  \item{organism}{ Currently, only "human", "mouse" and "yeast" supported. }
  \item{ont}{ One of "MF", "BP", and "CC" subontologies. }
  \item{pvalueCutoff}{ Cutoff value of pvalue. }
  \item{testDirection}{ one of "over" or "under" for over- or under-representation test. }
}

\value{
  A \code{enrichGOResult} instance.
}

\seealso{
	\code{\link{enrichGOResult-class}},
	\code{\link{compareCluster}}
}

\examples{
	#data(gcSample)
	#yy <- enrichGO(gcSample[[1]], organism="human", ont="BP", pvalueCutoff=0.01, testDirection="over")
	#head(summary(yy))
	#plot(yy)
}

\keyword{ manip }
