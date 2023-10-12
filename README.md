# clusterProfiler

<img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/clusterProfiler/clusterProfiler.png" height="200" align="right" />

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://img.shields.io/badge/release%20version-4.0.2-green.svg)](https://www.bioconductor.org/packages/clusterProfiler)
[![](https://img.shields.io/badge/devel%20version-4.1.3-green.svg)](https://github.com/guangchuangyu/clusterProfiler)
[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#since)

[![platform](http://www.bioconductor.org/shields/availability/devel/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#archives)
[![Build
Status](http://www.bioconductor.org/shields/build/devel/bioc/clusterProfiler.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/clusterProfiler/)
[![Linux/Mac Travis Build
Status](https://img.shields.io/travis/GuangchuangYu/clusterProfiler/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/clusterProfiler)
[![AppVeyor Build
Status](https://img.shields.io/appveyor/ci/Guangchuangyu/clusterProfiler/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/clusterProfiler)
[![codecov](https://codecov.io/gh/GuangchuangYu/clusterProfiler/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/clusterProfiler/)

<!--
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--08--16-green.svg)](https://github.com/GuangchuangYu/clusterProfiler/commits/master)
-->

-   [clusterProfiler](http://bioconductor.org/packages/clusterProfiler)
    supports exploring functional characteristics of both coding and
    non-coding genomics data for thousands of species with up-to-date
    gene annotation.
-   It provides a universal interface for gene functional annotation
    from a variety of sources and thus can be applied in diverse
    scenarios.
-   It provides a tidy interface to access, manipulate, and visualize
    enrichment results to help users achieve efficient data
    interpretation
-   Datasets obtained from multiple treatments and time points can be
    analyzed and compared in a single run, easily revealing functional
    consensus and differences among distinct conditions

For details, please visit
<https://yulab-smu.top/biomedical-knowledge-mining-book/>.

<img src="graphic-abstract-The-Innovation-2021.jpg" width="890"/>

## :writing_hand: Authors

Guangchuang YU <https://yulab-smu.top>

School of Basic Medical Sciences, Southern Medical University

[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=clusterProfiler&url=http://online.liebertpub.com/doi/abs/10.1089/omi.2011.0118&screen_name=guangchuangyu)
[![saythanks](https://img.shields.io/badge/say-thanks-ff69b4.svg)](https://saythanks.io/to/GuangchuangYu)
[![](https://img.shields.io/badge/follow%20me%20on-WeChat-green.svg)](https://guangchuangyu.github.io/blog_images/biobabble.jpg)

------------------------------------------------------------------------

If you use
[clusterProfiler](http://bioconductor.org/packages/clusterProfiler) in
published research, please cite the most appropriate paper(s) from this
list:

1.  T Wu<sup>#</sup>, E Hu<sup>#</sup>, S Xu, M Chen, P Guo, Z Dai, T
    Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo<sup>\*</sup>, **G
    Yu**<sup>\*</sup>. clusterProfiler 4.0: A universal enrichment tool
    for interpreting omics data. ***The Innovation***. 2021,
    2(3):100141. doi:
    [10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)
2.  **G Yu**<sup>\*</sup>. Gene Ontology Semantic Similarity Analysis
    Using GOSemSim. In: Kidder B. (eds) Stem Cell Transcriptional
    Networks. ***Methods in Molecular Biology***. 2020, 2117:207-215.
    Humana, New York, NY. doi:
    [10.1007/978-1-0716-0301-7_11](https://doi.org/10.1007/978-1-0716-0301-7_11)
3.  **G Yu**<sup>\*</sup>. Using meshes for MeSH term enrichment and
    semantic analyses. ***Bioinformatics***. 2018, 34(21):3766–3767.
    doi:
    [10.1093/bioinformatics/bty410](https://doi.org/10.1093/bioinformatics/bty410)
4.  **G Yu**, QY He<sup>\*</sup>. ReactomePA: an R/Bioconductor package
    for reactome pathway analysis and visualization. ***Molecular
    BioSystems***. 2016, 12(2):477-479. doi:
    [10.1039/C5MB00663E](https://doi.org/10.1039/C5MB00663E)
5.  **G Yu**<sup>\*</sup>, LG Wang, and QY He<sup>\*</sup>. ChIPseeker:
    an R/Bioconductor package for ChIP peak annotation, comparison and
    visualization. ***Bioinformatics***. 2015, 31(14):2382-2383. doi:
    [10.1093/bioinformatics/btv145](https://doi.org/10.1093/bioinformatics/btv145)
6.  **G Yu**<sup>\*</sup>, LG Wang, GR Yan, QY He<sup>\*</sup>. DOSE: an
    R/Bioconductor package for Disease Ontology Semantic and Enrichment
    analysis. ***Bioinformatics***. 2015, 31(4):608-609. doi:
    [10.1093/bioinformatics/btu684](https://doi.org/10.1093/bioinformatics/btu684)
7.  **G Yu**, LG Wang, Y Han and QY He<sup>\*</sup>. clusterProfiler: an
    R package for comparing biological themes among gene clusters.
    ***OMICS: A Journal of Integrative Biology***. 2012, 16(5):284-287.
    doi: [10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118)
8.  **G Yu**, F Li, Y Qin, X Bo<sup>\*</sup>, Y Wu, S Wang<sup>\*</sup>.
    GOSemSim: an R package for measuring semantic similarity among GO
    terms and gene products. ***Bioinformatics***. 2010, 26(7):976-978.
    doi:
    [10.1093/bioinformatics/btq064](https://doi.org/10.1093/bioinformatics/btq064)

<!--



 r badge_custom("1st most cited paper", "in OMICS", "green",
  "http://online.liebertpub.com/action/showMostCitedArticles?journalCode=omi")`
 r badge_custom("ESI", "Highly Cited Paper", "green")`
 r badge_doi("10.1089/omi.2011.0118", "green")`


------------------------------------------------------------------------

### Citation




<img src="https://guangchuangyu.github.io/software/citation_trend/clusterProfiler.png" width="890"/>


### Download stats

r badge_download_bioc("clusterProfiler")
r badge_bioc_download("clusterProfiler", "total", "blue")
r badge_bioc_download("clusterProfiler", "month", "blue")


<img src="https://guangchuangyu.github.io/software/dlstats/clusterProfiler.png" width="890"/>

-->
