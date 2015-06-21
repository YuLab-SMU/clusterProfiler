
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

##' dot plot of enrichResult
##'
##' 
##' @title dotplot
##' @param object an instance of enrichResult
##' @param x variable for x axis
##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param showCategory number of category
##' @param font.size font size
##' @param title plot title
##' @return ggplot object
##' @export
##' @author Guangchuang Yu
dotplot <- DOSE::dotplot

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
cnetplot <- DOSE:::cnetplot.enrichResult
