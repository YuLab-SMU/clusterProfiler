\name{compareCluster}
\alias{compareCluster}
\title{Compare gene clusters functional profile...}
\usage{compareCluster(geneClusters, fun=enrichGO, ...)}
\description{Compare gene clusters functional profile
Given a list of gene set, this function will compute profiles of each gene
cluster.}
\value{A \code{clusterProfResult} instance.}
\author{Guangchuang Yu \url{http://ygc.name}}
\seealso{\code{\link{compareClusterResult-class}}, \code{\link{groupGO}}
\code{\link{enrichGO}}}
\keyword{manip}
\arguments{\item{geneClusters}{a list of entrez gene id.}
\item{fun}{One of groupGO and enrichGO.}
\item{...}{Other arguments.}}
\examples{data(gcSample)
xx <- compareCluster(gcSample, fun=enrichKEGG, organism="human", pvalueCutoff=0.05)
#summary(xx)
#plot(xx, type="dot", caption="KEGG Enrichment Comparison")}
