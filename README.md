#  clusterProfiler

[![platform](http://www.bioconductor.org/shields/availability/devel/clusterProfiler.svg)](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#archives)
[![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/clusterProfiler.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/clusterProfiler/)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/clusterProfiler.svg)](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#since)
[![post](http://www.bioconductor.org/shields/posts/clusterProfiler.svg)](https://support.bioconductor.org/t/clusterProfiler/)
[![commit](http://www.bioconductor.org/shields/commits/bioc/clusterProfiler.svg)](http://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#svn_source)
[![download](http://www.bioconductor.org/shields/downloads/clusterProfiler.svg)](http://bioconductor.org/packages/stats/bioc/clusterProfiler.html)

This package implements methods to analyze and visualize functional profiles (GO and KEGG) of gene and gene clusters.

## Authors ##

Guangchuang YU, School of Public Health, The University of Hong Kong [http://ygc.name](http://ygc.name)

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
 `install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))`

## Documentation ##

+ [Why clusterProfiler fails](http://ygc.name/2014/08/07/why-clusterprofiler-fails/)
+ [use simplify to remove redundancy of enriched GO terms](http://ygc.name/2015/10/21/use-simplify-to-remove-redundancy-of-enriched-go-terms/)
+ [KEGG enrichment analysis with latest online data using clusterProfiler](http://ygc.name/2015/02/01/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/)
+ [DAVID functional analysis with clusterProfiler](http://ygc.name/2015/03/16/david-functional-analysis-with-clusterprofiler/)
+ [use clusterProfiler as an universal enrichment analysis tool](http://ygc.name/2015/05/11/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/)
+ [a formula interface for GeneOntology analysis <- bioinfoblog.it](http://bioinfoblog.it/2015/02/a-formula-interface-for-geneontology-analysis/)
+ [Enrichment map](http://ygc.name/2014/08/03/enrichment-map/)
+ [dotplot for enrichment result](http://ygc.name/2015/06/23/dotplot-for-enrichment-result/)
+ [functional enrichment for GTEx paper](http://ygc.name/2015/08/13/functional-enrichment-for-gtex-paper/)
+ [functional enrichment analysis with NGS data](http://ygc.name/2015/08/21/functional-enrichment-analysis-with-ngs-data/)


To view the vignette of `clusterProfiler` installed in your system, start `R` and enter:
```r
vignette("clusterProfiler", package="clusterProfiler")
```


## Bugs/Feature requests ##

 - If you have any, [let me know](https://github.com/GuangchuangYu/clusterProfiler/issues). Thx!

