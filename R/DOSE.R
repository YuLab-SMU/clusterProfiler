
##' enrichment map
##'
##' enrichMap
##' @title enrichMap
##' @param x gseaResult or enrichResult object
##' @param n maximum number of category to shown
##' @param fixed if set to FALSE, will invoke tkplot
##' @param vertex.label.font font size of vertex label
##' @param ... additional parameter
##' @return figure
##' @export
##' @author ygc
enrichMap <- DOSE::enrichMap



##' category-gene-net plot
##'
##' category gene association
##' @title cnetplot
##' @param x enrichResult object
##' @param showCategory number of category plotted
##' @param categorySize one of geneNum or pvalue
##' @param foldChange fold change of expression value
##' @param fixed logical
##' @param ... additional parameter
##' @return plot
##' @export
##' @author ygc
cnetplot <- DOSE::cnetplot
##cnetplot <- DOSE:::cnetplot.enrichResult
