#  clusterProfiler

[![platform](http://www.bioconductor.org/shields/availability/devel/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#archives)
[![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/clusterProfiler.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/clusterProfiler/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#since)
[![post](http://www.bioconductor.org/shields/posts/clusterProfiler.svg)](https://support.bioconductor.org/t/clusterProfiler/)
[![commit](http://www.bioconductor.org/shields/commits/bioc/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#svn_source)
[![download](http://www.bioconductor.org/shields/downloads/clusterProfiler.svg)](https://bioconductor.org/packages/stats/bioc/clusterProfiler)
[![gitter](https://img.shields.io/badge/GITTER-join%20chat-green.svg)](https://gitter.im/GuangchuangYu/Bioinformatics)


This package implements methods to analyze and visualize functional profiles of genomic coordinates (supported by [ChIPseeker](http://www.bioconductor.org/packages/ChIPseeker)), gene and gene clusters.

It supports both *hypergeometric test* and *Gene Set Enrichment Analysis* for many ontologies/pathways, including:

+ Disease Ontology (via [DOSE](http://www.bioconductor.org/packages/DOSE))
+ [Network of Cancer Gene](http://ncg.kcl.ac.uk/) (via [DOSE](http://www.bioconductor.org/packages/DOSE))
+ Gene Ontology (supports many species with GO annotation query online via [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/))
+ KEGG Pathway and Module with latest online data (supports more than 2000 species listed in <http://www.genome.jp/kegg/catalog/org_list.html>)
+ Reactome Pathway (via [ReactomePA](http://www.bioconductor.org/packages/ReactomePA))
+ DAVID (via [RDAVIDWebService](http://www.bioconductor.org/packages/RDAVIDWebService))
+ [Molecular Signatures Database](http://software.broadinstitute.org/gsea/msigdb)
  * hallmark gene sets
  * positional gene sets
  * curated gene sets
  * motif gene sets
  * computational gene sets
  * GO gene sets
  * oncogenic signatures
  * immunologic signatures
+ Other Annotations
  * from other sources (e.g. [DisGeNET](http://www.disgenet.org/web/DisGeNET/menu/home) as [an example](http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/))
  * user's annotation
  * customized ontology
  * and many others
  
*clusterProfiler* also provides several visualization methods to help interpreting enriched results, including:

+ barplot
+ cnetplot
+ dotplot
+ enrichMap
+ gseaplot
+ plotGOgraph (via [topGO](http://www.bioconductor.org/packages/topGO) package)
+ upsetplot
  
and several useful utilities:

+ bitr (Biological Id TranslatoR)
+ bitr_kegg (bitr using KEGG source)
+ compareCluster (biological theme comparison)
+ dropGO (screen out GO term of specific level or specific term)
+ go2ont (convert GO ID to Ontology)
+ go2term (convert GO ID to descriptive term)
+ gofilter (restrict result at specific GO level)
+ gsfilter (restrict result by gene set size)
+ search_kegg_organism (search kegg supported organism)
+ setReadable (convert IDs stored `enrichResult` object to gene symbol)
+ simplify (remove redundant GO terms, supported via [GOSemSim](http://www.bioconductor.org/packages/GOSemSim))

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

More documents can be found in <http://www.bioconductor.org/packages/DOSE>, <http://www.bioconductor.org/packages/clusterProfiler> and <http://guangchuangyu.github.io/clusterProfiler>.



## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/clusterProfiler/issues). Thx!

