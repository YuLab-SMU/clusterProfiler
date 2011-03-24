\name{compareClusterResult-class}
\docType{class}
\alias{compareClusterResult-class}
\alias{show,compareClusterResult-method}
\alias{summary,compareClusterResult-method}
\alias{plot,compareClusterResult-method}
\alias{show}
\alias{summary}
\alias{plot}


\title{Class "compareClusterResult"}
\description{
  This class represents the comparison result of gene clusters by GO categories at specific level or GO enrichment analysis.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("compareClusterResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{compareClusterResult}:}{containing a dataframe of functional profiles.}
	\item{\code{geneClusters}:}{ input gene clusters. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="compareClusterResult")}}
	\item{summary}{\code{signature(object="compareClusterResult")}}
	\item{plot}{\code{signature(x="compareClusterResult")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{groupGOResult}}
   \code{\linkS4class{enrichGOResult}}
  \code{\link{compareCluster}}
}

\keyword{classes}
