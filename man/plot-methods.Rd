\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,compareClusterResult,ANY-method}
\title{plot method}
\usage{
\S4method{plot}{compareClusterResult,ANY}(x, type = "dot", title = "",
  font.size = 12, showCategory = 5, by = "geneRatio",
  colorBy = "p.adjust")
}
\arguments{
  \item{...}{Additional argument list}

  \item{x}{compareClusterResult object}

  \item{type}{one of bar or dot}

  \item{title}{figure title}

  \item{font.size}{font size}

  \item{showCategory}{category numbers}

  \item{by}{one of geneRatio, Percentage or count}

  \item{colorBy}{one of pvalue or p.adjust}
}
\value{
plot
}
\description{
plot method generics
}
\author{
Guangchuang Yu \url{http://ygc.name}
}

