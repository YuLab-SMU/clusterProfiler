\name{enrichKEGG}
\alias{enrichKEGG}
\title{KEGG Enrichment Analysis of a gene set.
Given a vector of genes, this function will return the enrichment KEGG
categories with FDR control.}
\usage{
  enrichKEGG(gene, organism = "human", pvalueCutoff = 0.05,
    pAdjustMethod = "BH", universe, minGSSize = 5,
    qvalueCutoff = 0.2, readable = FALSE)
}
\arguments{
  \item{gene}{a vector of entrez gene id.}

  \item{organism}{Currently, only "human" and "mouse"
  supported.}

  \item{pvalueCutoff}{Cutoff value of pvalue.}

  \item{pAdjustMethod}{one of "holm", "hochberg", "hommel",
  "bonferroni", "BH", "BY", "fdr", "none"}

  \item{universe}{background genes}

  \item{qvalueCutoff}{qvalue cutoff}

  \item{minGSSize}{minimal size of genes annotated by
  Ontology term for testing.}

  \item{readable}{whether mapping gene ID to gene Name}
}
\value{
  A \code{enrichResult} instance.
}
\description{
  KEGG Enrichment Analysis of a gene set. Given a vector of
  genes, this function will return the enrichment KEGG
  categories with FDR control.
}
\examples{
data(gcSample)
	yy = enrichKEGG(gcSample[[5]], pvalueCutoff=0.01)
	head(summary(yy))
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

