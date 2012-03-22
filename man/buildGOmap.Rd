\name{buildGOmap}
\alias{buildGOmap}
\title{buildGOmap}
\usage{
  buildGOmap(gomap)
}
\arguments{
  \item{gomap}{data.frame with two columns names
  "entrezgene", and "go_accession"}
}
\value{
  files save in the the working directory
}
\description{
  building GO mapping files
}
\details{
  provided by a data.frame of gene and GO directly
  annotation file this function will building gene to GO
  and GO to gene mapping, with directly and undirectly
  annotation.
}
\author{
  Yu Guangchuang
}

