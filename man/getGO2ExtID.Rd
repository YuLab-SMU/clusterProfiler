\name{getGO2ExtID}
\alias{getGO2ExtID}
\title{query genes annotated by given GOIDs}
\usage{
  getGO2ExtID(GOID, organism)
}
\arguments{
  \item{GOID}{the query GO IDs}

  \item{organism}{one of human, mouse and yeast.}
}
\value{
  a list of gene IDs, the names of the list is the GOIDs
}
\description{
  provide a vector of GOIDs, and organism, this function
  will return the species specific gene list annotated by
  the given GOIDs.
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

