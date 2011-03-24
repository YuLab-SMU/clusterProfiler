\name{enrichGOResult-class}
\docType{class}
\alias{enrichGOResult-class}
\alias{show,enrichGOResult-method}
\alias{summary,enrichGOResult-method}
\alias{plot,enrichGOResult-method}
\alias{plot,GOHyperGResult-method}

\title{Class "enrichGOResult"}
\description{
  This class represents the result of GO enrichment analysis with FDR control.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("enrichGOResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{enrichGOResult}:}{containing a dataframe of GO enrichment analysis result.}
	\item{\code{pvalueCutoff}:}{ Cutoff value of pvalue. }
	\item{\code{testDirection}:}{ test direction for over- or under-representation. }
	\item{\code{Ont}:}{ sub-ontology for enrichment analysis }
	\item{\code{Organism}:}{ Organism. }
	\item{\code{Gene}:}{ input gene vector. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="enrichGOResult")}}
	\item{summary}{\code{signature(object="enrichGOResult")}}
	\item{plot}{\code{signature(x="enrichGOResult")}}
	\item{plot}{\code{signature(x="GOHyperGResult")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{compareClusterResult}}
  \code{\link{compareCluster}}
  \code{\link{enrichGO}}
}

\keyword{classes}
