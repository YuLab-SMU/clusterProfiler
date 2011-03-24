\name{enrichKEGGResult-class}
\docType{class}
\alias{enrichKEGGResult-class}
\alias{show,enrichKEGGResult-method}
\alias{summary,enrichKEGGResult-method}
\alias{plot,enrichKEGGResult-method}

\title{Class "enrichKEGGResult"}
\description{
  This class represents the result of KEGG enrichment analysis.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("enrichKEGGResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{enrichKEGGResult}:}{containing a dataframe of KEGG enrichment analysis result.}
	\item{\code{pvalueCutoff}:}{ Cutoff value of pvalue. }
	\item{\code{Organism}:}{ Organism. }
	\item{\code{Gene}:}{ input gene vector. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="enrichKEGGResult")}}
	\item{summary}{\code{signature(object="enrichKEGGResult")}}
	\item{plot}{\code{signature(x="enrichKEGGResult")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{compareClusterResult}}
  \code{\link{compareCluster}}
  \code{\link{enrichKEGG}}
}

\keyword{classes}
