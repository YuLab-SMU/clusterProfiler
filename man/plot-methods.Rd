\docType{methods}
\name{plot}
\alias{plot}
\title{plot method}
\arguments{
  \item{x}{A \code{compareClusterResult} instance.}

  \item{type}{one of "dot" and "bar".}

  \item{title}{graph title}

  \item{font.size}{graph font size}

  \item{showCategory}{numeric parameter, restrict the top
  categories for plotting.}

  \item{by}{one of "percentage" and "count"}

  \item{x}{A \code{enrichGOResult} instance.}

  \item{title}{graph title}

  \item{font.size}{graph font size}

  \item{showCategory}{number of KEGG categories to show.}

  \item{x}{A \code{enrichKEGGResult} instance.}

  \item{title}{graph title}

  \item{font.size}{graph font size}

  \item{showCategory}{number of KEGG categories to show.}

  \item{x}{A \code{groupGOResult} instance}

  \item{order}{logical parameter, order the result by
  *Count*.}

  \item{title}{graph title}

  \item{font.size}{graph font size}

  \item{showCategory}{number of GO categories to show.}

  \item{drop}{logical parameter, drop void category.}
}
\value{
  ggplot object

  ggplot object

  ggplot object

  ggplot object
}
\description{
  plot method for \code{compareClusterResult} instance

  plot method for \code{enrichGOResult} instance

  plot method for \code{enrichKEGGResult} instance

  plot method for \code{groupGOResult} instance
}
\author{
  Guangchuang Yu \url{http://ygc.name}

  Guangchuang Yu \url{http://ygc.name}

  Guangchuang Yu \url{http://ygc.name}

  Guangchuang Yu \url{http://ygc.name}
}

