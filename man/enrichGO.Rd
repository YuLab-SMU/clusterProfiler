\name{enrichGO}
\alias{enrichGO}
\title{GO Enrichment Analysis of a gene set.}
\usage{enrichGO(gene, organism="human", ont="MF", pvalueCutoff=0.01,
    readable=FALSE)}
\description{GO Enrichment Analysis of a gene set.
Given a vector of genes, this function will return the enrichment GO
categories with FDR control.}
\value{A \code{enrichGOResult} instance.}
\seealso{\code{\link{enrichGOResult-class}}, \code{\link{compareCluster}}}
\keyword{manip}
\arguments{\item{gene}{a vector of entrez gene id.}
\item{organism}{Currently, only "human", "mouse" and "yeast" supported.}
\item{ont}{One of "MF", "BP", and "CC" subontologies.}
\item{pvalueCutoff}{Cutoff value of pvalue.}
\item{readable}{if readable is TRUE, the gene IDs will mapping to gene
symbols.}}
\examples{#data(gcSample)
#yy <- enrichGO(gcSample[[1]], organism="human", ont="BP", pvalueCutoff=0.01)
#head(summary(yy))
#plot(yy)}
