\name{enrichKEGG}
\alias{enrichKEGG}
\title{KEGG Enrichment Analysis of a gene set.}
\description{
  Given a vector of genes, this function will return the enrichment KEGG categories with FDR control.
}
\usage{
	enrichKEGG(gene, organism="human", pvalueCutoff = 0.05, readable=FALSE)
}

\arguments{
  \item{gene}{ a vector of entrez gene id. }
  \item{organism}{ Currently, only "human" and "mouse" supported. }
  \item{pvalueCutoff}{ Cutoff value of pvalue. }
  \item{readable}{ if readable is TRUE, the gene IDs will mapping to gene symbols. }
}

\value{
  A \code{enrichKEGGResult} instance.
}

\seealso{
	\code{\link{enrichKEGGResult-class}},
	\code{\link{compareCluster}}
}

\examples{
	data(gcSample)
	yy = enrichKEGG(gcSample[[5]], pvalueCutoff=0.01)
	head(summary(yy))
	#plot(yy)
}
\keyword{ manip }
