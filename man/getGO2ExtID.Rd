\name{getGO2ExtID}
\alias{getGO2ExtID}
\title{query genes annotated by given GOIDs}
\usage{getGO2ExtID(GOID, organism)}
\description{provide a vector of GOIDs, and organism, this function will return the species specific gene list annotated by the given GOIDs.}
\value{a list of gene IDs, the names of the list is the GOIDs}
\author{Guangchuang Yu}
\arguments{\item{GOID}{the query GO IDs}
\item{organism}{one of human, mouse and yeast.}}
