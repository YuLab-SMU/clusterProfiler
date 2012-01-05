\name{plot.categoryNet}
\alias{plot.categoryNet}
\title{plot gene net by categories}
\usage{
  plot.categoryNet(inputList, categorySize = "geneNum",
    showCategory = 5, pvalue = NULL, foldChange = NULL,
    output = "fixed")
}
\arguments{
  \item{inputList}{a list of gene IDs}

  \item{categorySize}{setting category size}

  \item{showCategory}{number of categories to plot}

  \item{pvalue}{pvalue}

  \item{foldChange}{foldChange}

  \item{output}{output type}
}
\value{
  plotted igraph object.
}
\description{
  plot function of gene Concept Net.
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

