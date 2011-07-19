

#' statistical analysis and visulization of functional profiles for genes and
#' gene clusters
#' The package implements methods to analyze and visualize functional profiles
#' of gene and gene clusters.
#' 
#' This package is designed to compare gene clusters functional profiles.
#' 
#' \tabular{ll}{ Package: \tab clusterProfiler\cr Type: \tab Package\cr
#' Version: \tab 0.99.0\cr Date: \tab 03-15-2011\cr biocViews:\tab GO,
#' Clustering, Visulization\cr Depends:\tab AnnotationDbi, GO.db, org.Hs.eg.db,
#' GOstats, ggplot2, plyr, methods\cr Suggests:\tab GOSemSim\cr License: \tab
#' Artistic-2.0\cr }
#' 
#' @name clusterProfiler-package
#' @aliases clusterProfiler-package clusterProfiler
#' @docType package
#' @author Guangchuang Yu
#' 
#' Maintainer: Guangchuang Yu <guangchuangyu@@gmail.com>
#' @seealso \linkS4class{compareClusterResult}, \linkS4class{groupGOResult}
#'   \linkS4class{enrichGOResult}
#' @keywords package
NULL





#' Class "compareClusterResult"
#' This class represents the comparison result of gene clusters by GO
#' categories at specific level or GO enrichment analysis.
#' 
#' 
#' @name compareClusterResult-class
#' @aliases compareClusterResult-class show,compareClusterResult-method
#'   summary,compareClusterResult-method plot,compareClusterResult-method show
#'   summary plot
#' @docType class
#' @author Guangchuang Yu <guangchuangyu@@gmail.com>
#' @seealso \code{\linkS4class{groupGOResult}}
#'   \code{\linkS4class{enrichGOResult}} \code{\link{compareCluster}}
#' @keywords classes
NULL





#' Datasets
#' gcSample contains a sample of gene clusters.
#' 
#' 
#' @name DataSet
#' @aliases gcSample
#' @docType data
#' @keywords datasets
NULL





#' Class "enrichGOResult"
#' This class represents the result of GO enrichment analysis with FDR control.
#' 
#' 
#' @name enrichGOResult-class
#' @aliases enrichGOResult-class show,enrichGOResult-method
#'   summary,enrichGOResult-method plot,enrichGOResult-method
#'   plot,GOHyperGResult-method
#' @docType class
#' @author Guangchuang Yu <guangchuangyu@@gmail.com>
#' @seealso \code{\linkS4class{compareClusterResult}}
#'   \code{\link{compareCluster}} \code{\link{enrichGO}}
#' @keywords classes
NULL





#' Class "enrichKEGGResult"
#' This class represents the result of KEGG enrichment analysis.
#' 
#' 
#' @name enrichKEGGResult-class
#' @aliases enrichKEGGResult-class show,enrichKEGGResult-method
#'   summary,enrichKEGGResult-method plot,enrichKEGGResult-method
#' @docType class
#' @author Guangchuang Yu <guangchuangyu@@gmail.com>
#' @seealso \code{\linkS4class{compareClusterResult}}
#'   \code{\link{compareCluster}} \code{\link{enrichKEGG}}
#' @keywords classes
NULL





#' Class "groupGOResult"
#' This class represents the result of functional Profiles of a set of gene at
#' specific GO level.
#' 
#' 
#' @name groupGOResult-class
#' @aliases groupGOResult-class show,groupGOResult-method
#'   summary,groupGOResult-method plot,groupGOResult-method
#' @docType class
#' @author Guangchuang Yu <guangchuangyu@@gmail.com>
#' @seealso \code{\linkS4class{compareClusterResult}}
#'   \code{\link{compareCluster}} \code{\link{groupGO}}
#' @keywords classes
NULL



