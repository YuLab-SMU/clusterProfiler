#  clusterProfiler

[![platform](http://www.bioconductor.org/shields/availability/devel/clusterProfiler.svg)](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#archives)
[![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/clusterProfiler.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/clusterProfiler/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/clusterProfiler.svg)](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#since)
[![post](http://www.bioconductor.org/shields/posts/clusterProfiler.svg)](https://support.bioconductor.org/t/clusterProfiler/)
[![commit](http://www.bioconductor.org/shields/commits/bioc/clusterProfiler.svg)](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#svn_source)
[![download](http://www.bioconductor.org/shields/downloads/clusterProfiler.svg)](http://bioconductor.org/packages/stats/bioc/clusterProfiler.html)
[![gitter](https://img.shields.io/badge/GITTER-join%20chat-green.svg)](https://gitter.im/GuangchuangYu/Bioinformatics)


This package implements methods to analyze and visualize functional profiles of genomic coordinates (supported by [ChIPseeker](http://www.bioconductor.org/packages/ChIPseeker)), gene and gene clusters.

It supports both *hypergeometric test* and *Gene Set Enrichment Analysis* for many ontologies/pathways, including:

+ Disease Ontology (via [DOSE](http://www.bioconductor.org/packages/DOSE))
+ [Network of Cancer Gene](http://ncg.kcl.ac.uk/) (via [DOSE](http://www.bioconductor.org/packages/DOSE))
+ Gene Ontology (supports many species with GO annotation query online via [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/))
+ KEGG Pathway and Module with latest online data (supports more than 2000 species listed in <http://www.genome.jp/kegg/catalog/org_list.html>)
+ Reactome Pathway (via [ReactomePA](http://www.bioconductor.org/packages/ReactomePA))
+ [Molecular Signatures Database](http://software.broadinstitute.org/gsea/msigdb)
  * hallmark gene sets
  * positional gene sets
  * curated gene sets
  * motif gene sets
  * computational gene sets
  * GO gene sets
  * oncogenic signatures
  * immunologic signatures
+ DAVID (via [RDAVIDWebService](http://www.bioconductor.org/packages/RDAVIDWebService))
+ Other Annotations
  * from other sources (e.g. [DisGeNET](http://www.disgenet.org/web/DisGeNET/menu/home) as [an example](http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/))
  * user's annotation
  * customized ontology
  * and many others
  
*clusterProfiler* also provides several visualization methods to help interpreting enriched results, including:

+ plotGOgraph (via [topGO](http://www.bioconductor.org/packages/topGO) package)
+ barplot
+ dotplot
+ cnetplot
+ enrichMap
+ gseaplot

and several useful utilities:

+ gofilter (to restrict the result at specific GO level)
+ dropGO (to screen out GO term of specific level or specific term)
+ simplify (to remove redundant GO terms)

A discussion forum can be found at <https://groups.google.com/forum/#!forum/clusterprofiler>.

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong <http://guangchuangyu.github.io>

## Citation ##

Please cite the following article when using `clusterProfiler`:

```
Yu G, Wang L, Han Y and He Q. 
clusterProfiler: an R package for comparing biological themes among gene clusters.
OMICS: A Journal of Integrative Biology, 2012, 16(5):284-287. 
```

URL: [http://online.liebertpub.com/doi/abs/10.1089/omi.2011.0118](http://online.liebertpub.com/doi/abs/10.1089/omi.2011.0118)

## License ##

All source code is copyright, under the Artistic-2.0 License.
For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Installation ##

To install:
 * the latest released version:
   `biocLite("clusterProfiler")`
 * the latest development version:
 `devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))`

## Documentation ##

To view the vignette of `clusterProfiler` installed in your system, start `R` and enter:
```r
vignette("clusterProfiler", package="clusterProfiler")
```

More documents can be found in <http://www.bioconductor.org/packages/DOSE>, <http://www.bioconductor.org/packages/ReactomePA> and <http://guangchuangyu.github.io/tags/clusterprofiler>.



## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/clusterProfiler/issues). Thx!

