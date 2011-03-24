\name{compareCluster}
\alias{compareCluster}
\title{Compare gene clusters functional profile}
\description{
  Given a list of gene set, this function will compute profiles of each gene cluster.
}
\usage{
	compareCluster(geneClusters, fun=enrichGO, ...)
}

\arguments{
  \item{geneClusters}{ a list of entrez gene id. }
  \item{fun}{ One of groupGO and enrichGO. }
  \item{...}{ Other arguments. }
}

\value{
  A \code{clusterProfResult} instance.
}

\seealso{
	\code{\link{compareClusterResult-class}},
	\code{\link{groupGO}}
	\code{\link{enrichGO}}
}

\examples{
	data(gcSample)
	xx <- compareCluster(gcSample, fun=enrichKEGG, organism="human", pvalueCutoff=0.05)
	#summary(xx)
	#plot(xx, type="dot", caption="KEGG Enrichment Comparison")
}
\keyword{ manip }
