---
html_preview: False
output:
  md_document:
    variant: markdown
---

clusterProfiler: statistical analysis and visualization of functional profiles for genes and gene clusters
==========================================================================================================

<!-- AddToAny BEGIN -->
<div class="a2a_kit a2a_kit_size_32 a2a_default_style">

<a class="a2a_dd" href="//www.addtoany.com/share"></a>
<a class="a2a_button_facebook"></a> <a class="a2a_button_twitter"></a>
<a class="a2a_button_google_plus"></a>
<a class="a2a_button_pinterest"></a> <a class="a2a_button_reddit"></a>
<a class="a2a_button_sina_weibo"></a> <a class="a2a_button_wechat"></a>
<a class="a2a_button_douban"></a>

</div>

<script async src="//static.addtoany.com/menu/page.js"></script>
<!-- AddToAny END -->
<link rel="stylesheet" href="https://guangchuangyu.github.io/css/font-awesome.min.css">
<link rel="stylesheet" href="https://guangchuangyu.github.io/css/academicons.min.css">

[![releaseVersion](https://img.shields.io/badge/release%20version-3.2.11-blue.svg?style=flat)](https://bioconductor.org/packages/clusterProfiler)
[![develVersion](https://img.shields.io/badge/devel%20version-3.3.6-blue.svg?style=flat)](https://github.com/guangchuangyu/clusterProfiler)
[![total](https://img.shields.io/badge/downloads-34261/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/clusterProfiler)
[![month](https://img.shields.io/badge/downloads-1401/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/clusterProfiler)

The `clusterProfiler` package implements methods to analyze and
visualize functional profiles of genomic coordinates (supported by
ChIPseeker), gene and gene clusters.

`clusterProfiler` is released within the
[Bioconductor](https://bioconductor.org/packages/clusterProfiler)
project and the source code is hosted on
<a href="https://github.com/GuangchuangYu/clusterProfiler"><i class="fa fa-github fa-lg"></i>
GitHub</a>.

<i class="fa fa-user"></i> Author
---------------------------------

Guangchuang Yu, School of Public Health, The University of Hong Kong.

<a href="https://twitter.com/guangchuangyu"><i class="fa fa-twitter fa-3x"></i></a>
<a href="https://guangchuangyu.github.io/blog_images/biobabble.jpg"><i class="fa fa-wechat fa-3x"></i></a>
<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=Guangchuang+Yu[Author+-+Full]"><i class="ai ai-pubmed ai-3x"></i></a>
<a href="https://scholar.google.com.hk/citations?user=DO5oG40AAAAJ&hl=en"><i class="ai ai-google-scholar ai-3x"></i></a>
<a href="https://orcid.org/0000-0002-6485-8781"><i class="ai ai-orcid ai-3x"></i></a>
<a href="https://impactstory.org/u/0000-0002-6485-8781"><i class="ai ai-impactstory ai-3x"></i></a>

<i class="fa fa-book"></i> Citation
-----------------------------------

Please cite the following article when using `clusterProfiler`:

[![doi](https://img.shields.io/badge/doi-10.1089/omi.2011.0118-blue.svg?style=flat)](http://dx.doi.org/10.1089/omi.2011.0118)
[![citation](https://img.shields.io/badge/cited%20by-160-blue.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=2349076811020942117)
[![Altmetric](https://img.shields.io/badge/Altmetric-23-blue.svg?style=flat)](https://www.altmetric.com/details/681089)

**Yu G**, Wang L, Han Y and He Q<sup>\*</sup>. clusterProfiler: an R
package for comparing biological themes among gene clusters. ***OMICS: A
Journal of Integrative Biology***. 2012, 16(5):284-287.

<i class="fa fa-pencil"></i> Featured Articles
----------------------------------------------

<img src="https://guangchuangyu.github.io/featured_img/clusterProfiler/elife-02077-fig5-v1.jpg" width="650">

<i class="fa fa-hand-o-right"></i> Find out more on
<i class="fa fa-pencil"></i> [Featured
Articles](https://guangchuangyu.github.io/clusterProfiler/featuredArticles/).

<i class="fa fa-download"></i> Installation
-------------------------------------------

Install `clusterProfiler` is easy, follow the guide on the [Bioconductor
page](https://bioconductor.org/packages/clusterProfiler/):

``` {.r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("clusterProfiler")
```

<i class="fa fa-cogs"></i> Overview
-----------------------------------

#### <i class="fa fa-angle-double-right"></i> Supported Analyses

-   Over-Representation Analysis
-   Gene Set Enrichment Analysis
-   Biological theme comparison

#### <i class="fa fa-angle-double-right"></i> Supported ontologies/pathways

-   Disease Ontology (via
    [DOSE](https://www.bioconductor.org/packages/DOSE))
-   [Network of Cancer Gene](http://ncg.kcl.ac.uk/) (via
    [DOSE](https://www.bioconductor.org/packages/DOSE))
-   [DisGeNET](http://www.disgenet.org/web/DisGeNET/menu/home) (via
    [DOSE](https://www.bioconductor.org/packages/DOSE))
-   Gene Ontology (supports many species with GO annotation query online
    via
    [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/))
-   KEGG Pathway and Module with latest online data (supports more than
    4000 species listed in
    <http://www.genome.jp/kegg/catalog/org_list.html>)
-   MeSH enrichment analysis (via
    [meshes](https://www.bioconductor.org/packages/meshes))
-   Reactome Pathway (via
    [ReactomePA](https://www.bioconductor.org/packages/ReactomePA))
-   DAVID (via
    [RDAVIDWebService](https://www.bioconductor.org/packages/RDAVIDWebService))
-   [Molecular Signatures
    Database](http://software.broadinstitute.org/gsea/msigdb)
    -   hallmark gene sets
    -   positional gene sets
    -   curated gene sets
    -   motif gene sets
    -   computational gene sets
    -   GO gene sets
    -   oncogenic signatures
    -   immunologic signatures
-   Other Annotations
    -   from other sources (e.g.
        [DisGeNET](http://www.disgenet.org/web/DisGeNET/menu/home) as
        [an
        example](https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/))
    -   user's annotation
    -   customized ontology
    -   and many others

#### <i class="fa fa-angle-double-right"></i> Visualization

-   barplot
-   cnetplot
-   dotplot
-   enrichMap
-   gseaplot
-   plotGOgraph (via
    [topGO](https://www.bioconductor.org/packages/topGO) package)
-   upsetplot

#### <i class="fa fa-angle-double-right"></i> Useful utilities:

-   bitr (Biological Id TranslatoR)
-   compareCluster (biological theme comparison)
-   dropGO (screen out GO term of specific level or specific term)
-   go2ont (convert GO ID to Ontology)
-   go2term (convert GO ID to descriptive term)
-   gofilter (restrict result at specific GO level)
-   gsfilter (restrict result by gene set size)
-   simplify (remove redundant GO terms, supported via
    [GOSemSim](https://www.bioconductor.org/packages/GOSemSim))

<i class="fa fa-hand-o-right"></i> Find out details and examples on
<i class="fa fa-book"></i>
[Documentation](https://guangchuangyu.github.io/clusterProfiler/documentation/).

<i class="fa fa-wrench"></i> Related Tools
------------------------------------------

<ul class="fa-ul">
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/DOSE">DOSE</a> for Disease Ontology Semantic and Enrichment analyses</li>
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/GOSemSim">GOSemSim</a> for GO semantic similarity measurement</li>
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/meshes">meshes</a> for MeSH Enrichment and Semantic analysis</li>
    <li><i class="fa-li fa fa-angle-double-right"></i><a href="https://guangchuangyu.github.io/ReactomePA">ReactomePA</a> for Reactome pathway analysis</li>

</ul>
<i class="fa fa-hand-o-right"></i> Find out more on
[projects](https://guangchuangyu.github.io/#projects).

<i class="fa fa-code-fork"></i> Projects that depend on *clusterProfiler*
-------------------------------------------------------------------------

#### <i class="fa fa-angle-double-right"></i> Bioconductor packages

-   [bioCancer](https://www.bioconductor.org/packages/bioCancer):
    Interactive Multi-Omics Cancers Data Visualization and Analysis
-   [debrowser](https://www.bioconductor.org/packages/debrowser):
    Interactive Differential Expresion Analysis Browser
-   [eegc](https://www.bioconductor.org/packages/eegc): Engineering
    Evaluation by Gene Categorization (eegc)
-   [LINC](https://www.bioconductor.org/packages/LINC): co-expression of
    lincRNAs and protein-coding genes
-   [MoonlightR](https://www.bioconductor.org/packages/MoonlightR):
    Identify oncogenes and tumor suppressor genes from omics data
-   [TCGAbiolinks](https://www.bioconductor.org/packages/TCGAbiolinks):
    An R/Bioconductor package for integrative analysis with TCGA data
-   [TCGAbiolinksGUI](https://www.bioconductor.org/packages/TCGAbiolinksGUI):
    "TCGAbiolinksGUI: A Graphical User Interface to analyze cancer
    molecular and clinical data"

#### <i class="fa fa-angle-double-right"></i> Other applications

-   [APOSTL](https://github.com/bornea/APOSTL): An Interactive Galaxy
    Pipeline for Reproducible Analysis of Affinity Proteomics Data

<i class="fa fa-comment"></i> Feedback
--------------------------------------

<ul class="fa-ul">
    <li><i class="fa-li fa fa-hand-o-right"></i> Please make sure you have followed <a href="https://guangchuangyu.github.io/2016/07/how-to-bug-author/"><strong>the important guide</strong></a> before posting any issue/question</li>
    <li><i class="fa-li fa fa-bug"></i> For bugs or feature requests, please post to <i class="fa fa-github-alt"></i> [github issue](https://github.com/GuangchuangYu/clusterProfiler/issues)</li>
    <li><i class="fa-li fa fa-support"></i>  For user questions, please post to [Bioconductor support site](https://support.bioconductor.org/) and [Biostars](https://www.biostars.org/). We are following every post tagged with **clusterProfiler**</li>
    <li><i class="fa-li fa fa-commenting"></i> Join the group chat on <a href="https://twitter.com/hashtag/clusterProfiler"><i class="fa fa-twitter fa-lg"></i></a> and <a href="http://huati.weibo.com/k/clusterProfiler"><i class="fa fa-weibo fa-lg"></i></a></li>

</ul>
