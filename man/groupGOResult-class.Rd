\name{groupGOResult-class}
\docType{class}
\alias{groupGOResult-class}
\alias{show,groupGOResult-method}
\alias{summary,groupGOResult-method}
\alias{plot,groupGOResult-method}

\title{Class "groupGOResult"}
\description{
  This class represents the result of functional Profiles of a set of gene at specific GO level.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("groupGOResult", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{groupGOResult}:}{containing a dataframe of functional profiles.}
	\item{\code{Ont}:}{Ontology used in analysis. }
	\item{\code{Level}:}{GO Level. }
	\item{\code{Organism}:}{Organism. }
	\item{\code{Gene}:}{input Gene vector. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object="groupGOResult")}}
	\item{summary}{\code{signature(object="groupGOResult")}}
	\item{plot}{\code{signature(x="groupGOResult")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{compareClusterResult}}
  \code{\link{compareCluster}}
  \code{\link{groupGO}}
}

\keyword{classes}
