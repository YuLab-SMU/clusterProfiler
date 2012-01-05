\name{plotting.clusterProfile}
\alias{plotting.clusterProfile}
\title{plotting-clusterProfile}
\usage{
  plotting.clusterProfile(clProf.reshape.df, type = "dot",
    by = "percentage", title = "", font.size = 12)
}
\arguments{
  \item{clProf.reshape.df}{data frame of compareCluster
  result}

  \item{type}{one of dot and bar}

  \item{by}{one of percentage and count}

  \item{title}{graph title}

  \item{font.size}{graph font size}
}
\value{
  ggplot object
}
\description{
  Internal plot function for plotting compareClusterResult
}
\author{
  Guangchuang Yu \url{http://ygc.name}
}

