#' A list with 6980 K genes of human gut genome
#'
#' All 6980 K genes were reported on
#' "An integrated catalog of reference genes in the
#'  human gut microbiome.[J]. Nature Biotechnology, 2014."
#'
#'  data(hgmlist)
#'  kk <- enrichKEGG(gene = hgmlist[1:20], organism = "ko", universe = hgmlist )
#'  dotplot(kk)
#'
#' @format A list with 6980 K genes
"hgmlist"
