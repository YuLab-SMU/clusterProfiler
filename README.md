clusterProfiler
===============

[![](https://img.shields.io/badge/release%20version-3.4.1-green.svg?style=flat)](https://bioconductor.org/packages/clusterProfiler) [![](https://img.shields.io/badge/devel%20version-3.5.2-green.svg?style=flat)](https://github.com/guangchuangyu/clusterProfiler) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#since) [![](https://img.shields.io/badge/download-39609/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/clusterProfiler) [![](https://img.shields.io/badge/download-2029/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/clusterProfiler)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/GuangchuangYu/clusterProfiler/branch/master/graph/badge.svg)](https://codecov.io/gh/GuangchuangYu/clusterProfiler/) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--05--16-green.svg)](https://github.com/GuangchuangYu/clusterProfiler/commits/master) [![commit](http://www.bioconductor.org/shields/commits/bioc/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#svn_source) [![GitHub forks](https://img.shields.io/github/forks/GuangchuangYu/clusterProfiler.svg)](https://github.com/GuangchuangYu/clusterProfiler/network) [![GitHub stars](https://img.shields.io/github/stars/GuangchuangYu/clusterProfiler.svg)](https://github.com/GuangchuangYu/clusterProfiler/stargazers)

[![platform](http://www.bioconductor.org/shields/availability/devel/clusterProfiler.svg)](https://www.bioconductor.org/packages/devel/bioc/html/clusterProfiler.html#archives) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/clusterProfiler.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/clusterProfiler/) [![Linux/Mac Travis Build Status](https://img.shields.io/travis/GuangchuangYu/clusterProfiler/master.svg?label=Mac%20OSX%20%26%20Linux)](https://travis-ci.org/GuangchuangYu/clusterProfiler) [![AppVeyor Build Status](https://img.shields.io/appveyor/ci/Guangchuangyu/clusterProfiler/master.svg?label=Windows)](https://ci.appveyor.com/project/GuangchuangYu/clusterProfiler)

This package implements methods to analyze and visualize functional profiles of genomic coordinates (supported by [ChIPseeker](http://www.bioconductor.org/packages/ChIPseeker)), gene and gene clusters.

For details, please visit our project website, <https://guangchuangyu.github.io/clusterProfiler>.

-   [Documentation](https://guangchuangyu.github.io/clusterProfiler/documentation/)
-   [Featured Articles](https://guangchuangyu.github.io/clusterProfiler/featuredArticles/)
-   [Feedback](https://guangchuangyu.github.io/clusterProfiler/#feedback)

[![Twitter](https://img.shields.io/twitter/url/https/github.com/GuangchuangYu/clusterProfiler.svg?style=social)](https://twitter.com/intent/tweet?hashtags=clusterProfiler&url=http://online.liebertpub.com/doi/abs/10.1089/omi.2011.0118&screen_name=guangchuangyu)

------------------------------------------------------------------------

Please cite the following article when using `clusterProfiler`:

***Yu G***, Wang L, Han Y and He Q<sup>\*</sup>. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

[![](https://img.shields.io/badge/doi-10.1089/omi.2011.0118-green.svg?style=flat)](http://dx.doi.org/10.1089/omi.2011.0118) [![](https://img.shields.io/badge/Altmetric-22-green.svg?style=flat)](https://www.altmetric.com/details/681089)

------------------------------------------------------------------------

### Citation

[![](https://img.shields.io/badge/2nd%20most%20cited%20paper-in%20OMICS-green.svg?style=flat)](http://online.liebertpub.com/action/showMostCitedArticles?journalCode=omi) [![citation](https://img.shields.io/badge/cited%20by-198-green.svg?style=flat)](https://scholar.google.com.hk/scholar?oi=bibs&hl=en&cites=2349076811020942117) [![](https://img.shields.io/badge/cited%20in%20Web%20of%20Science%20Core%20Collection--green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000303653300007) [![](https://img.shields.io/badge/ESI-Highly%20Cited%20Paper-green.svg?style=flat)](http://apps.webofknowledge.com/InboundService.do?mode=FullRecord&customersID=RID&IsProductCode=Yes&product=WOS&Init=Yes&Func=Frame&DestFail=http%3A%2F%2Fwww.webofknowledge.com&action=retrieve&SrcApp=RID&SrcAuth=RID&SID=Y2CXu6nry8nDQZcUy1w&UT=WOS%3A000303653300007)

       +-+---------+---------+---------+---------+---------+---+
    80 +                                         *             +
       |                                                       |
       |                                                       |
    60 +                                                       +
       |                                                       |
       |                                                   *   |
    40 +                                                       +
       |                                                       |
       |                               *                       |
    20 +                     *                                 +
       |           *                                           |
       | *                                                     |
       +-+---------+---------+---------+---------+---------+---+
       2012      2013      2014      2015      2016      2017   

### Download stats

[![download](http://www.bioconductor.org/shields/downloads/clusterProfiler.svg)](https://bioconductor.org/packages/stats/bioc/clusterProfiler) [![](https://img.shields.io/badge/download-39609/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/clusterProfiler) [![](https://img.shields.io/badge/download-2029/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/clusterProfiler)

         ++-------------+-------------+------------+-------------+-------------+-------------+---------+
         |                                                                                      *      |
         |                                                                                             |
    4000 +                                                                                             +
         |                                                                                             |
         |                                                                                        *    |
         |                                                                      *                      |
         |                                                                         *                   |
    3000 +                                                                    ** *                     +
         |                                                                          *       * *        |
         |                                                                                   * *       |
         |                                                                              * **           |
         |                                                                           *                 |
    2000 +                                                                            **               +
         |                                                                                             |
         |                                                                  **                         |
         |                                                                 *                           |
    1000 +                                                           *****                             +
         |                                                    *                                        |
         |                              *  *     ***    ***  * ** * *                                  |
         |                 **    * **  * ** * ***    * *   *     *                                     |
         |      *     *  **  * ** *   *               *                                                |
       0 +   **  *****  *                                                                              +
         ++-------------+-------------+------------+-------------+-------------+-------------+---------+
        2011          2012          2013         2014          2015          2016          2017
