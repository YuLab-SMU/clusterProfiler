\name{viewKEGG}
\alias{viewKEGG}
\title{viewKEGG function is for visualize KEGG pathways
works with enrichResult object to visualize enriched KEGG pathway}
\usage{
  viewKEGG(obj, pathwayID, foldChange, color.low = "green",
    color.high = "red", kegg.native = TRUE,
    out.suffix = "clusterProfiler")
}
\arguments{
  \item{obj}{enrichResult object}

  \item{pathwayID}{pathway ID or index}

  \item{foldChange}{fold change values}

  \item{color.low}{color of low foldChange genes}

  \item{color.high}{color of high foldChange genes}

  \item{kegg.native}{logical}

  \item{out.suffix}{suffix of output file}
}
\description{
  viewKEGG function is for visualize KEGG pathways works
  with enrichResult object to visualize enriched KEGG
  pathway
}

