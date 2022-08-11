
<!--

TODO:

+ Uniprot to KO, <https://www.biostars.org/p/415837/>
  - e.g. <https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+up:A0A059ZTB3>
+ aggregate multiple p values
  - <https://en.wikipedia.org/wiki/Fisher's_method>
  - e.g. independent test of cancer sample 1 vs control 1 and cancer sample 2 vs control 2, then combine results
  - e.g. independent test for different omics data and then combine results


## KEGG

+ <https://www.kegg.jp/kegg/rest/keggapi.html>
+ [x] get KEGG organism list using api
  - <http://rest.kegg.jp/list/organism>
  - Currently, the `get_kegg_species()` function parses html page





-->

# clusterProfiler 4.5.2

+ support passing a GSON object to `enricher(USER_DATA)` and `GSEA(USER_DATA)` (2022-8-01, Mon)
+ `gson_kegg_mapper()` allows building a gson object from outputs of KEGG Mapper service (2022-07-29, Fri, #492)
+ fix `show` method for `compareClusterResult` (2022-06-21, Tue, #473)
+ `gson_KEGG()` download latest KEGG and output a GSON object (2022-06-08, Wed)
+ support passing a GSON object to `gseKEGG(organism)` 
+ support passing a GSON object to `enrichKEGG(organism)` (2022-06-06, Mon)

# clusterProfiler 4.5.1

+ follow KEGG api upgrade that change from http to https (2022-06-06, Mon)
+ use 'wininet' to download KEGG data when `.Platform$OS.type = "windows"` (2022-06-03, Fri)
+ mv `read.gmt` and `read.gmt.wp` to the 'gson' package and reexport these two functions from 'gson' (2022-04-28, Thu)
+ fix `compareCluster` when fun = `enrichPathway`(2022-4-28, Thu)

# clusterProfiler 4.4.0

+ Bioconductor 3.15 release

# clusterProfiler 4.3.4

+ fix `enrichGO` , `gseGO` and `groupGO` when `keyType = 'SYMBOL'` && `readable=TRUE`(2022-4-9, Sat)

# clusterProfiler 4.3.3

+ parse GAF file to prepare GO annotation data (esp for proteomic study) (2022-03-08, Tue, #397, #418, #421, #442)
+ bug fixed in `compareCluster()` (2022-01-27, Thu, #424)

# clusterProfiler 4.3.2

+ bug fixed in `extract_params()` (2022-01-12, Wed, #392, @amcdavid)
+ make `simplify()` works for `gseGO()` in `compareCluster()`
+ support formula interface for GSEA methods in `compareCluster()` (2022-01-04, Tue, @altairwei, #416)

# clusterProfiler 4.3.1

+ `compareCluster()` supports GSEA algorithm (2021-12-11, Sat)
+ update error message of `download.KEGG.Path()` and `download.KEGG.Module()`(2021-11-21, Sun)
+ update `simplify()` function to support `ont = ALL` (2021-10-27, Wed)

# clusterProfiler 4.2.0

+ Bioconductor 3.14 release

# clusterProfiler 4.1.4

+ import `yulab.utils` (2021-08-20, Fri)

# clusterProfiler 4.1.3

+ Remove Human Gut Microbiome dataset as the functionalities are provided in <https://github.com/YuLab-SMU/MicrobiomeProfiler> (2021-08-15, Sun)

# clusterProfiler 4.1.2

+ update citation and DESCRIPTION (2021-08-15, Sun)
+ update kegg_species.rda and allow online download using KEGG api (2021-08-14, Sat)

# clusterProfiler 4.1.1

+ add citation (new paper published on The Innovation) (2021-07-04, Sun)

# clusterProfiler 4.0.0

+ Bioconductor 3.13 release

# clusterProfiler 3.99.1

+ Add new data set, `DE_GSE8057`, which contains DE genes obtained from GSE8057 (2020-03-08, Mon)

# clusterProfiler 3.99.0

+ Add KEGG enrichment analysis of Human Gut Microbiome data (2021-02-20, Sat)

# clusterProfiler 3.19.1

+ setting default timeout to 300 for downloads (2021-02-05, Fri)
+ fixed download method setting 
+ capable of setting KEGG download method via `options(clusterProfiler.download.method = METHOD)` (2020-12-31, Thu)

# clusterProfiler 3.18.0

+ Bioconductor 3.12 release (2020-10-28, Wed)

# clusterProfiler 3.17.5

+ update `[[.compareClusterResult` (2020-10-14, Wed)

# clusterProfiler 3.17.3

+ internal suports of enrichment analyses using WikiPathways (2020-09-09, Wed)
  - `enrichWP` for ORA analysis
  - `gseWP` for GSEA analysis
  - `get_wp_organisms` for listing supported organisms
  - `read.gmt.wp` for parsing gmt file downloaded from wikiPathways 

# clusterProfiler 3.17.2

+ use `libcurl` if capable (2020-09-08, Tue)
  - <https://github.com/YuLab-SMU/clusterProfiler/pull/290>

# clusterProfiler 3.17.1

+ bug fixed of `extract_params` (2020-08-18, Tue)
  - <https://github.com/YuLab-SMU/clusterProfiler/issues/282>

# clusterProfiler 3.16.0

+ Bioconductor 3.11 release

# clusterProfiler 3.15.3

+ incorporate clusterProfiler.dplyr (2020-03-12, Thu)
  + arrange, filter, group_by, mutate, rename, select, slice and summarize

# clusterProfiler 3.15.2

+ remove `Suggests` of `KEGG.db` as it will be deprecated in Bioconductor 3.11 (2020-01-14, Tue)
+ optimize `enrichGO` to use less memory (2019-12-13, Fri)
+ re-implement `read.gmt` without using GSEABase, and my own version is much more fasta :) 

# clusterProfiler 3.15.1

+ e.g. user can pass `fun=enrichGO` to `compareCluster` without quoting `enrichGO` (2019-12-02, Mon)
+ add `keytype` and `readable` info in `compareCluster` output
+ mv `compareClusterResult` class defintion to `DOSE` (2019-11-02, Sat)
+ mv `fortify`, `barplot` and  `dotplot` for `compareClusterResult` to `enrichplot`.

# clusterProfiler 3.14.0

+ Bioconductor 3.10 release

# clusterProfiler 3.12.0

+ Bioconductor 3.9 release

# clusterProfiler 3.11.1

+ `asis` parameter in `[.compareClusterResult` (2018-12-24, Mon)
  - <https://github.com/GuangchuangYu/enrichplot/issues/17>

# clusterProfiler 3.10.0

+ Bioconductor 3.8 release

# clusterProfiler 3.9.2

+ re-export `DOSE::gsfilter` and `DOSE::setReadable` (2018-05-25, Fri)

# clusterProfiler 3.9.1

+ change color scheme of dotplot of compareClusterResult back to red->purple
  (2018-05-17, Thu)
  - <https://support.bioconductor.org/p/108996/>

# clusterProfiler 3.8.0

+ Bioconductor 3.7 release

# clusterProfiler 3.7.1

+ uniprot_get function (2018-01-30, Tue)
+ import enrichplot (2018-01-29, Mon)
