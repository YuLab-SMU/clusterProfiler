
# clusterProfiler 4.21.4

+ **A protein name that maps to a gene is no longer treated as a family label.** TRANSFAC names a nuclear receptor by its protein -- `GR`, `PXR`, `FXR`, `HIF2A` -- and none of those are HGNC symbols, so `mechanism_is_family_label()` called each one a family. The cost fell on the perturbation channel, which skips family labels because a database indexes symbols: `GR` was never queried against GPSAdb and could not reach `supported` on a glucocorticoid experiment, where it is the one answer that is certain. On the airway dexamethasone list this had left `GR` at `single_source` level B while GPSAdb held two `NR3C1` datasets moving the list at padj 1e-22. With the classification corrected, `GR` is `supported` at level A and the list carries 24 regulators instead of 22 (2026-10-26, Mon)
+ Regulators are looked up in perturbation databases under every name they are known by. The collections say `GR`, the database says `NR3C1`; querying only the collection's spelling returned nothing, which is indistinguishable from "never perturbed". Both directions are now tried, and a name with no alias costs no extra request (2026-10-26, Mon)

+ **The perturbation sweep no longer runs concurrently.** The GPSAdb endpoint returns *different data* when queried in parallel and reports nothing unusual -- twenty concurrent lookups of one dataset gave the complete answer 16 times and a partial one 4 times (231 of 379 significant genes; another run 53 of 91), while six serial lookups agreed exactly. Nothing distinguishes the partial response: same request, same status, no warning. The consequence was visible in the verdicts, not in a message: the supported set flipped between runs (`LEF1` or `HOXA4`, never both), and a dataset could win a regulator on a truncated measurement. The sweep is I/O bound, so this costs about 45% more wall time (241 s against 173 s on the airway list) and buys a reproducible answer -- three consecutive runs now agree exactly. Set `options(clusterProfiler_perturb_workers = 4)` to take the trade knowingly (2026-10-26, Mon)
+ A p-value tie between two datasets for one regulator is now settled by the data rather than by row order. `padj` underflows to exactly 0 for well-powered datasets, so `SIX1`'s D25077 and D25078 tied and `order()` fell back on arrival order. The tie-break is now: system match, then p-value, then the number of genes moved, then the dataset id -- the last of which is unique, so ordering is total and no tie can survive to depend on how the rows happened to arrive (2026-10-26, Mon)

+ `gseKEGG()` reports the same `category` / `subcategory` columns as `enrichKEGG()`. The classification comes from the KEGG PATHWAY hierarchy that `enrichKEGG()` has always appended; `gseKEGG()` never called `append_kegg_category()`, so the same pathways carried the labels in one table and not in the other, and a reader that filtered on `category` silently behaved differently depending on which function produced the result (#787). `gseMKEGG()` used to call the same function on a result whose `setType` was already `"MKEGG"`, which only printed "not a KEGG enrichment result" on every run: `kegg_category_data()` carries pathway IDs (01100, ...) and no module IDs, so there is nothing to append for modules and the call is gone
+ `append_kegg_category()` documents that the classification is a cached copy of <https://www.genome.jp/kegg/pathway.html> rather than a live query: a pathway released after the cache was built comes back `NA` (2026-10-26, Mon, #787)
+ `mechanism_graph()` turns an `evidenceBundle` into a `mechgraph` object and `plot_mechanism()` draws it, so the verdicts can be seen rather than read out of a table. Regulators are coloured by the grade they earned -- A arrived with a measured perturbation, B with a recorded mechanism, C with an association -- and each edge names the collection that produced it. On the airway list that is 22 supported regulators linked to 67 of the input genes. A link means "this collection lists this gene", not "this factor regulates it": which members of a target repertoire are real targets is exactly what the analysis could not settle, and the caption says so (2026-10-26, Mon)

+ `simplify()` removes redundancy from **any** collection, not only GO: `measure = "JC"` compares the term gene sets by Jaccard index and `measure = "overlap"` by shared genes over the smaller set. Both read `@geneSets`, so they need no ontology and no `GOSemSimDATA`, which is what KEGG, Reactome and MSigDB results were missing -- `simplify()` used to stop with "only applied to output from gseGO and enrichGO" and the error now points at `measure = "JC"`. Reactome's hierarchy is the case that needed it: 66 significant terms on the `geneList` example compress to 41 at `cutoff = 0.7` and 31 at `cutoff = 0.5`, leaving Cell Cycle Checkpoints, Mitotic Prometaphase, APC/C-mediated degradation and Kinesins instead of every nested step beneath them. The grouping is the one `simplify_internal()` already used; ties go to the larger gene set, mirroring the GO version's preference for the more specific term (2026-10-26, Mon, #743)
+ **`evidence_similarity()` and `evidence_biotrigger()` are exported but are not part of the bundle**, after measuring what they actually contribute. On the airway list `biotrigger` raised `supported` from 21 to 25 and added no transcription factor: the additions were TREM2 (a microglial receptor), SORL1 and TMEM106B (neurodegeneration), AK2 and PI4KB (housekeeping enzymes) and C22orf46P (an unnamed locus) -- genes with no bearing on airway smooth muscle, and often reached through the same GPSAdb dataset that `perturbagen` had already used. A channel that raises a count without improving the mechanism is not an improvement, so it stays callable on its own but out of `upstream` (2026-10-26, Mon)
+ `consensus()` now requires the corroborating sources to cite **different datasets**, not merely different source names. `perturbagen` and `biotrigger` are two views of one GPSAdb experiment -- one says "knocking this down changes these genes", the other "perturbing this matches this list" -- so when both cited D25689 the agreement was the server agreeing with itself. Seventeen verdicts were exactly that (2026-10-26, Mon)
+ `options(clusterProfiler.progress =)` reports what a slow run is doing: `"major"` (the default) prints one line per evidence source, `"all"` also times each step, `"none"` is silent. A run with perturbation evidence takes minutes, mostly waiting, and the caller had no way to tell a slow run from a hung one (2026-10-26, Mon)

+ **A failed query is no longer reported as an empty answer.** Every source this package consults can fail, and the failures used to be caught with `error = function(e) NULL`, which made them indistinguishable from a successful query that found nothing. The two lead to opposite conclusions: "no knockdown dataset moves these genes" is a finding, "GPSAdb rejected the request" is not. A source that fails is now recorded as `failed` in `$diagnostics$failures` and shown by `print()`; a run in which **every** consulted source fails aborts rather than reporting `absent` for every entity; and a check that cannot run aborts rather than passing. This is the shape all four of this release's silent bugs shared (2026-10-25, Sun)
+ `$regulators` and `consensus()` now carry `evidence_kind` (`measured`, `predicted`, `observed` or `none`), and `supported` requires **at least one measured source**, corroborated by a second independent source. Counting sources is not counting evidence: a C3:TFT motif set and a miRDB repertoire are both predictions, so their agreement is one observation confirmed by the same reasoning twice, and a prediction plus a tissue signature is two things that were never perturbed. Independent measurements still corroborate when they share a kind -- two knockdown datasets are the strongest position available, and counting distinct kinds alone under-reported exactly that case. Measured on the airway list: `tf+mirna` gave 177 findings and **zero** `supported`, `cellstate+tissue` 379 and zero, `tf+cellstate` 492 and zero, while `tf+perturbagen` gave 22 (2026-10-25, Sun)
+ `mechanism_is_family_label()` no longer reads a hyphenated symbol as a family. TRANSFAC writes `NKX2_5` where HGNC writes `NKX2-5`; testing only the literal spelling called 7 such genes family labels, which also stripped them of their `measured` evidence kind since they arrive from GTRD (2026-10-25, Sun)
+ `mechanism_rule_for()` reports which xref rule a set identifier matches, walking the table in the same order as grading so the rule and the level can never disagree. `mechanism_kind_of()` maps that to the evidence kind (2026-10-25, Sun)

+ `collapse_findings()` groups findings whose gene sets largely overlap and keeps the strongest of each group, so a result reads as a few statements instead of every set that crossed a threshold. It is a separate step rather than something `consensus()` always does, because the need is not uniform: measured on the airway list, transcription-factor findings barely overlap (median pairwise Jaccard 0.000 over 155), cell-state findings likewise (0.008 over 337), and only tissue signatures are genuinely redundant (17% of pairs above 0.5; collapsing keeps 17 of 42). A collapse step that is not needed discards findings that say different things (2026-10-25, Sun)
+ `collapse_findings()` reports `collapsed` and `absorbed`, attributed at the moment a finding is dropped. Counting "how many were similar" afterwards double-counts, because one dropped finding can be similar to two that were kept (2026-10-25, Sun)

+ **The headline verdict could be produced by a symbol collision.** `consensus()` now refuses `supported` for a label that names different things in different sources, and says so in the sentence. Verified case: TRANSFAC's `NF1` is Nuclear Factor I (the `NF1_Q6` motif in MSigDB TFT_LEGACY), while the HGNC symbol `NF1` is neurofibromin -- the protein GPSAdb's D27128 knocks out. Matched on the string, the motif term and the knockout dataset "agreed" and made `NF1` a level-A `supported` finding on the airway list, for two molecules that have nothing to do with each other. `mechanism_symbol_collisions()` holds the known cases; roughly 610 of the 1114 C3:TFT terms are motif-style and may use a family name that collides with an unrelated symbol, so the table is a `shortcut:` until a TRANSFAC-to-HGNC map exists (2026-10-25, Sun)
+ Perturbation candidates are now the regulators whose target sets were **significant**, ordered by significance and capped (default 30, `options(clusterProfiler_perturb_candidates =)`), instead of every named regulator. Querying all 613 named regulators in the 900-gene airway list meant 613 requests and a run that never finished; the focused version answers in 35 s. How many were queried, how many were eligible and how many datasets were found is printed and carried on the returned table (2026-10-25, Sun)
+ GPSAdb has no data for most transcription-factor collection labels because those labels are family or complex names (`NFAT`, `AP1`, `MEF2`, `OCT1`), not gene symbols. Of the 30 strongest candidates for the airway list only 3 had a dataset. The `measured:` line in `print()` reports the coverage so a handful of hits is not mistaken for a handful of candidates (2026-10-25, Sun)
+ microRNA collections are no longer tested by default. They are predicted seed matches, and on the airway list 7 miRNAs were significant among the UP-regulated genes and 0 among the down ones -- the opposite of what a repressor predicts -- so the direction that would make the result interpretable is absent. `upstream = "mirna"` asks for them explicitly (2026-10-25, Sun)

+ `evidence(upstream = )` also accepts `"cellstate"` (C8) and `"tissue"` (C4:CGN), so a gene list can be compared with cell-state and tissue signatures as well as with regulators. These collections are read differently on purpose: their set identifiers already name the entity (`CAR_HPX`, `HAY_BONE_MARROW_MONOCYTE`), so they are passed through unchanged instead of through the transcription-factor rules, which would have turned a cardiac signature into the gene HPX and dropped the first word of a two-part state name. Only a TF label can be "unnamed", and both sources are graded C: `supported` would be a category error for an entity that means "this signature resembles that state" rather than "this factor may drive these genes" (2026-10-25, Sun)
+ `mechanism_entity()` collapses a set identifier according to the source that produced it, which is what makes the four collection families safe to mix in one bundle (2026-10-25, Sun)

+ `evidence_perturbation()` asks GPSAdb, over its MCP endpoint, which measured perturbation moves the input genes, and is the only source that can reach level A: a knockdown is an experiment, not an annotation. The question it answers is deliberately narrow -- whether knocking a factor down moves these genes and in which direction -- because `gpsadb_de()` returns the 500 most significant rows rather than the measured set, so an enrichment against them would be meaningless. `evidence(upstream = "perturbagen", design = list(treatment = "knockdown"))` attaches the table as `$perturbagen`; the result rows carry the same `regulator` label as the collections, which is what lets a factor named by both its target sets and a knockdown be one entity with two sources (2026-10-25, Sun)
+ `mechanism_alias_regulator()` brings a regulator label to one spelling before agreement is counted. C3:TFT names the hypoxia factor `HIF1` (the legacy symbol) while GPSAdb names it `HIF1A`; without the alias, a factor supported by both its target sets and a knockdown is two findings, and `consensus()` could never report the agreement a mechanism needs. With it, `consensus()` reaches `supported` for the first time on real data (2026-10-25, Sun)
+ `evidence()` accepts `upstream = "perturbagen"` alongside `"tf"` and `"mirna"`, and `design$background` names entities that belong to the experimental setting rather than the findings, so that `unsupported_claims` does not report them as unevidenced (2026-10-25, Sun)

+ `consensus()` combines an `evidenceBundle` into one verdict per regulator, separating `supported` (at least `min_independent` sources agree and one is graded B or above) from `single_source`, `absent` and `untested`. Agreement is counted only for an identical regulator label: `FOXO3` from a TF target collection and `MIR23C` from a miRNA collection that happen to share targets are two `single_source` findings, not one `supported` one. A source whose `coverage_tested` falls below `coverage_floor` reports `untested` for everything, because a source that never examined these genes cannot report an absence (2026-10-25, Sun)
+ `mechanism_writable()` turns a verdict and a level into the sentence the evidence permits, so a level C finding is offered wording that describes a correlation and never wording that asserts a mechanism. `evidence()` carries the table in `$writable` (2026-10-25, Sun)
+ `evidence_ppi()` grades STRING edges by the channel that supports them: `escore` (experiment) or `dscore` (database) at 0.4 or above gives B, anything else C, and a network edge never reaches A because no STRING channel measures a perturbation. `pscore` is deliberately not used -- it is 0 even for NFKB1-NFKBIA, which carries `escore = 0.999` -- so a `pscore` rule would grade a published physical interaction as an association. `evidence(network = TRUE)` attaches the table as `$ppi`; both the taxID default (human, since `getPPI()` can only infer a species from an enrichment result) and the channel thresholds are documented on the arguments (2026-10-25, Sun)

+ `evidence()` grades a gene list against *upstream participant* collections -- transcription factor targets, microRNA targets -- and returns an `evidenceBundle` whose `$regulators` table carries a per-regulator evidence level rather than a bare p-value. Each collection is corrected within itself and the column is named `padj_within_source`: C3:TFT holds about 1100 sets and C3:MIR about 2600 with different target definitions, so pooling them into one correction would let a large collection dilute a small one (2026-10-25, Sun)
+ `evidence()` reports `$diagnostics` with counts whose denominators differ and which must not be conflated: `n_input` is the number of genes handed in; `n_tested` is how many could enter the test at all -- inputs found in the universe and in some set, which is the `DESize` that `enrichit` divides `GeneRatio` by; `n_sig_covered` is how many a significant regulator covers. `$diagnostics$per_gene` additionally labels each gene `supported`, `absent` (annotated, but no significant regulator) or `untested` (not annotated at all), because an absence of evidence is not evidence of absence (2026-10-25, Sun)
+ `mechanism_sets()` returns an MSigDB collection as a `GSON` -- `"tf"` (C3:TFT), `"mirna"` (C3:MIR), `"perturbagen"` (C2:CGP), `"cellstate"` (C8) or `"tissue"` (C4:CGN) -- for `enricher(gson = )` or `GSEA(gson = )`. The GMT is fetched once into the session cache and reused only when its recorded release matches the `release` argument, so changing the release can never silently reuse an old collection; the file names are literal because they cannot be derived from the collection name (`c8.all`, `c4.cgn`, and `c2.cgn` does not exist). Every set is kept in symbol namespace, and the returned object must be passed through `gson=`, not `TERM2GENE=`: only the former carries species, keytype and gsname into the result (2026-10-25, Sun)
+ `mechanism_xref()` publishes the level rules behind `evidence()`: A for a measured perturbation, B for a recorded mechanism such as a ChIP-supported target set, C for an association such as a predicted target repertoire or a co-expression pattern. Unnamed motif sets such as `AACTTT_UNKNOWN` are reported as C regardless of the rule that matched, because `TATAAA_TATA_01` would otherwise inherit the B grade of the generic motif rule and a consensus sequence would be quotable as a mechanism (2026-10-25, Sun)
+ the session cache directory is no longer forced to `tempdir()` on attach. `.onAttach` overwrote a directory set in `.Rprofile`, and because `tempdir()` differs in every session it also made downloaded collections miss their cache; the option is now only filled when nobody has set it (2026-10-25, Sun)

+ `interpret()` gains `evidence_status()`, a deterministic grounding check: every gene symbol the interpretation names is classified as `observed` (a member gene of the enrichment result) or `out_of_evidence` (a valid symbol that is not), across `overview`, `key_mechanisms`, `regulatory_drivers`, `crosstalk`, `hypothesis`, `narrative`, `reasoning` and `network_evidence`. No model is called, so, unlike a second model, the check is exact and does not also mark grounded genes as uncertain; it labels a hypothesis as inferred rather than observed (2026-09-26, Sat)
+ `interpret()`, `interpret_agent()` and `interpret_hierarchical()` stamp their results with an execution record, read back with the new `provenance()` (printed by `print.interpret_provenance`; `provenance(x)$ok` says whether the call completed): the model actually addressed, the decoding parameters requested and whether the response confirmed them, SHA-256 digests of the system prompt, user prompt, raw response and input evidence, token usage (including Anthropic cached prompt tokens), `finish_reason`, per-call status and package versions. A remote model is not a versionable dependency, and without the record a call that failed at the transport layer arrives as an ordinary result object, indistinguishable from a model that had nothing to say (2026-09-03, Thu)
+ `interpret()` reads per-gene statistics through the new `gene_stats` argument, a differential-expression table (limma `topTable()`, DESeq2 `results()`, edgeR `topTags()` or Seurat `FindMarkers()`; the gene, log fold change and adjusted p-value columns are recognised by name or row names) and also supplies `gene_fold_change` when that is `NULL`. The genes of the enriched terms are then listed with log2 fold change and FDR, and a "Gene-Level Statistics" section adds the number of genes changed at FDR < 0.05 and the 50 most significant up- and down-regulated genes. A `gseaResult` without `gene_fold_change` gives its leading-edge genes the values of the ranked list it was run on (`x@geneList`). Without `gene_stats` the prompt for an `enrichResult` is unchanged (2026-10-18, Sun)
+ `interpret()`, `interpret_agent()` and `interpret_hierarchical()` no longer truncate the genes they send: every input gene is listed (most up-regulated first when `gene_fold_change` is supplied; previously the first 20), every fold change of the sent terms' genes is given (previously the 20 largest), the leading-edge genes of a `gseaResult` receive their fold changes, and an agent cluster without enriched terms lists every input gene (previously 50). The prompt headers drop "Top". `n_pathways` stays `NULL` (every significant term, one representative per theme) and its `n_pathways = 20` examples are removed: on five public lists a lower cap named fewer themes of the full result, so a short context should lower `reduce_threshold` before capping (2026-10-11, Sun)
+ `interpret()` now interprets an enrichment result and nothing else, so its `task` and `prior` arguments are removed; `interpret_hierarchical()` keeps its cell-type labels and loses its `task` argument. Combining enrichment evidence with the results of other tools (a SingleR or scGPT label, a TF-activity ranking, an enriched protein complex) moves to `interpret_agent(prior = )`, whose Synthesizer returns `conclusion` (the answer to the question the context poses) and, with `prior`, `prior_verdict` (Confirmed, Refined, Corrected or Not assessable) and `prior_rationale`. With `gene_fold_change` the top marker genes are also listed to the Detective and the Synthesizer (2026-09-30, Wed)
+ `reduce_redundancy()` and `reduce_go_redundancy()` are exported, and `interpret()` takes `reduce`/`reduce_threshold` (`reduce = TRUE` by default), so every significant theme reaches the model once: GO results are reduced by Wang semantic similarity on GOSemSim, KEGG and other gene-set collections by Jaccard overlap of member genes. Both are a no-op for data frames and lists. `reduce_go_redundancy()` keeps its rrvgo-style `rrvgo_*` attributes for backwards compatibility (2026-09-20, Sun)
+ `interpret()` and `interpret_agent()` take `reasoning`, passed through to the provider (e.g. `list(type = "enabled", budget_tokens = 4096)`), and record it in the execution stamp as `reasoning_requested`, alongside `reasoning_evidence`, which says whether the response actually carried a reasoning block. Left `NULL`, the request is byte-identical to one made before the argument existed (2026-09-20, Sun)
+ `interpret_agent()` gains an optional Reviewer agent (`review = TRUE`, `FALSE` by default) that audits the draft: the result gains `review` and `overview_original`, and `review` reports grounding separately from the draft (2026-09-03, Thu)
+ `interpret()` and `interpret_agent()` now warn when `add_ppi = TRUE` and the STRING query fails, instead of silently interpreting without PPI evidence (2026-09-30, Wed)
+ fix `interpret()` to read gene names from the object's symbols rather than from a mismatched identifier column; an unparsed model response is now reported instead of returning `NULL`, and the failure handler no longer destroys the failure it is handling (2026-09-26, Sat; 2026-09-03, Thu)
+ `digest` is now in `Suggests`, used for the prompt and response digests of `provenance()` (2026-10-23, Fri)

# clusterProfiler 4.21.3

+ fix `get_kegg_species()` and `kegg_taxa()` for KEGG's retirement of the `organism` special database name: `https://rest.kegg.jp/list/organism` now answers HTTP 400 while the rest of the REST API keeps working, so the weekly "Update KEGG Data" workflow and any user calling `search_kegg_organism(use_internal_data = FALSE)` failed with "cannot open URL ... HTTP status was '400 Bad Request'". Both now read `/list/genome`, which returns the same organisms as `<T number>\t<code>; <scientific name> (<common name>)`; the species list is parsed from the name and code of that single column. An organism without a common name keeps its scientific name in `common_name`, as it did before (2026-10-18, Sun)
+ `restore_hybrid_result()` and `restore_entrez_result_ids()` map a gene set's IDs back to the source `keyType` with one vectorized pass instead of one `labels[ids]` + `unname()` per set. On a BP result (~5,900 tested sets, 390k ID slots) that loop cost 7.5 s per call and `compareCluster()` calls it once per cluster; the vectorized form takes 0.2 s and is byte-identical, including list names, empty sets and `geneID` order (pinned by a test). `hybrid_ids_for_input()` likewise replaces a per-id `mapping$source == id` scan with a single split and named lookup (negligible wall time at cluster-sized inputs, but no longer quadratic in the mapping size). Measured on the PBMC 10k 13-cluster SYMBOL `compareCluster()` together with the enrichit fixes: 140 s -> about 70 s (run-to-run 69-78 s), with 7,489 rows, p-values (difference 0), Counts and geneID sets identical (2026-10-18, Sun)
+ `enrichGO()` caches the hybrid GSON it builds for a non-`ENTREZID` `keyType` when some input IDs cannot be mapped to Entrez. Without a `universe`, that build reads every key of the OrgDb and therefore depends only on (organism, ont, keyType), but `compareCluster()` calls `enrichGO()` once per cluster and rebuilt it every time (7.7 s per call on human GO BP, PBMC 10k, 13 clusters, 2026-10-18). The hybrid path itself is unchanged and still fires ONLY when an input ID has no Entrez mapping -- an all-mappable input keeps the existing fast path and never builds hybrid data -- while a universe-scoped or `evidence`-filtered call neither reads nor writes this cache (2026-10-18, Sun)
+ `read.gmt.pc()` now reads the fields of a Pathway Commons gene set description by their labels rather than by splitting the string at every `"; "`. A pathway name that holds a semicolon, such as `Phospholipase C-mediated cascade; FGFR2`, is no longer cut short and its `datasource`, `organism` and `idtype` stay in their own columns; previously the extra pieces were merged into the last field, so the fields were shifted one to the left and nothing said so. A description that declares fewer fields than expected returns `NA` for the missing one instead of failing the whole read, and `enrichPC()` and `gsePC()` take the species from the first gene set that declares an organism (#789, 2026-10-08, Thu)
+ fix `enrichGO()` with a non-`ENTREZID` `keyType` when some source IDs cannot be mapped to Entrez: the memory-saving Entrez annotation remains the fast path, while unmapped source IDs are added through direct native-key annotations and analyzed in one hybrid GSON; when `universe` is omitted, the effective annotated hybrid background is constructed and reported by message instead of silently dropping IDs (#829, 2026-10-03, Sat)
+ restore the original source `keyType` consistently after the Entrez optimization and readable conversion, including result metadata and projected gene-set memberships (#829, 2026-10-03, Sat)

# clusterProfiler 4.21.2

+ `enrichMKEGG()` and `gseMKEGG()` accept a `GSON` object for `organism`, as `enrichKEGG()`/`gseKEGG()` already did, so KEGG **module** enrichment can run against a locally built annotation (`gson_KEGG()`) when KEGG is unreachable; they also now stop with "organism should be a species name or a GSON object" instead of an unrelated error from the organism mapper. The `organism` documentation for the KEGG functions now mentions the GSON option at all (2026-09-22, Tue, #623)
+ `enrichGO()` takes an `evidence` argument to restrict the GO annotation by evidence code before the analysis, e.g. `evidence = c("IDA", "IPI", "IMP")` to keep manually supported annotations or `evidence = "IEA"` for electronic ones only; dropping `IEP` is the usual way to avoid circular reasoning when the gene sets come from co-expression clustering. Because the annotation (and hence the background) changes, the result is not simply a subset of the unfiltered run, and filtered runs are not cached (2026-09-22, Tue, #160)
+ `enrichGO()` no longer aborts with `missing value where TRUE/FALSE needed` when an OrgDb's organism cannot be identified — an AnnotationHub OrgDb, for example. The cached-GO-data check compared the cached organism against `get_organism(OrgDb)` in a bare `if()`, and that comparison is `NA` for such an OrgDb, so the second call in a session failed with an error that named none of this; an unidentifiable organism is now treated as "not cached" and the annotation is rebuilt (2026-09-22, Tue, #232)
+ `Gff2GeneTable()` now reads Ensembl/GTF-style GFF3 as well as NCBI/RefSeq files: it accepts a `gene_id=` attribute (Ensembl) in addition to `GeneID=` (NCBI) and takes the gene symbol from `Name=`/`gene_name=` when `gene=` is absent. Previously an Ensembl file left every key `NA`, and merging two all-`NA` key columns is a cartesian product — useful-looking but wrong on a small file and `negative length vectors are not allowed` on a real genome. A file with neither identifier style now stops with a message saying what is expected (2026-09-22, Tue, #193)
+ `add_go_level()` appends a `level` column giving the GO level of every enriched term (1 = the ontology root), so results can be kept to a band of levels with e.g. `subset(res, level >= 3 & level <= 6)`; the level is taken from the ontology the term itself belongs to, so `ont = "ALL"` and `compareCluster()` results are labelled term by term (2026-09-22, Tue, #793)
+ `simplify()` now infers the ontology from the GO IDs when the result object does not carry one, instead of refusing: this covers a `gseaResult` whose `@setType` holds a collection name rather than an ontology (MSigDB C5, a gson GO file, a `TERM2GENE` built from GO annotations) as well as an `enrichResult` from `enricher()`, which previously required assigning `x@ontology` by hand; results spanning several ontologies are simplified as `GOALL`, and non-GO gene sets are still rejected with a message that says so (2026-09-22, Tue, #753, #369)
+ `compareCluster()` now warns when a non-character `universe` is supplied: it wraps each per-cluster call in `suppressMessages()`, which used to swallow the underlying enrichment function's "`universe` is not in character and will be ignored" message, so a numeric universe appeared to be accepted while having no effect (2026-09-22, Tue, #654)
+ `gsePC()` and `gseWP()` now expose an `eps` argument and forward it to `GSEA()` (and on to `enrichit::gsea_gson()`), matching `GSEA()`, `gseGO()`, `gseKEGG()` and `gseMKEGG()` (2026-09-22, Tue)
+ fix `gson_cpd()` to use the `compound` KEGG REST endpoint; `enrichKEGG(organism = 'cpd')` had been failing since KEGG retired `/link/cpd/pathway`, which now returns HTTP 400 (2026-09-22, Tue, #828)
+ `kegg_rest()` now reports a failed download where it happens, telling an unreachable network apart from an HTTP error status and from an empty response; previously these surfaced much later as an obscure error (2026-09-22, Tue)
+ `download.KEGG.Path()` and `download.KEGG.Module()` now stop with an actionable message when the downloaded annotation is empty, or when the gene/pathway table shares no ID with the pathway-name table, instead of returning an empty annotation that downstream surfaced as a misleading `"No gene can be mapped"`; this replaces an `is.null()` guard that could never fire (2026-09-22, Tue, #146, #378, #561, #646)
+ align GSEA significance filtering with the historical behavior of clusterProfiler <= 4.18.x: `pvalueCutoff` now requires both the raw p-value and the adjusted p-value to pass the cutoff, so significant-pathway counts from `GSEA()`/`gseGO()`/`gseKEGG()`/`gseMKEGG()` are again comparable with the DOSE/fgsea backend (2026-08-14, Thu)
+ restore the `seed` argument on `GSEA()`, `gseGO()`, `gseMKEGG()` and `gseKEGG()` and forward it to `enrichit::gsea_gson()`, so permutation-based GSEA results can be reproduced with a fixed seed (or `set.seed()` before the call); the argument had been dropped when the enrichit engine was introduced (2026-08-14, Thu)
+ fix `get_data_from_KEGG_db()` to strip species suffix from pathway names when `use_internal_data=TRUE`, matching the online path behaviour (2026-08-04, Tue, #783)
+ fix `enrichGO()` to check whether the OrgDb supports `ENTREZID` before attempting the ENTREZID-remapping optimisation; custom OrgDb packages without `ENTREZID` keytype (e.g., non-model organisms) now fall back to the original `get_GO_data` path instead of throwing an error (2026-08-04, Tue, #823)

# clusterProfiler 4.21.1

+ support `UniProt accession -> KO` conversion in `bitr_kegg()` via KEGG REST and export `gson_KO()` for KEGG Orthology workflows (2026-06-24, Wed)
+ add `nseGO()`, `nseKEGG()`, `nseMKEGG()`, `nseWP()`, `mnseGO()`, `mnseKEGG()`, `mnseMKEGG()`, and `mnseWP()` as high-level topology-aware enrichment wrappers powered by `enrichit`, with metadata backfilling and lightweight wrapper tests (2026-06-24, Wed)
+ re-export `aggregate_omics()`, `aggregate_enrichment()`, `harmonize_ids()`, and `select_features_for_ora()` from enrichit for high-level multi-omics workflows (2026-06-24, Wed)
+ `get_ppi_network` to download and cache the full STRING PPI network for a specific species (2026-06-23, Tue)
  - the `getPPI` is now alias to `get_ppi`
+ fix GSEA wrappers to restore `eps`, forward extra arguments such as `sampleSize` to `enrichit::gsea_gson()`, and align the GSEA method documentation with the current enrichit interface (2026-06-22, Mon, #822)
+ re-export `bayes_enrich()` and `bayes_summary()` from enrichit (>= 0.1.5) for Bayesian term selection and posterior-ranked summaries of ORA enrichment results (2026-06-16, Tue)

# clusterProfiler 4.20.0

+ Bioconductor RELEASE_3_23 (2026-04-29, Wed)

# clusterProfiler 4.19.8

+ fix: map non-ENTREZID universe in enrichGO (2026-04-22, Wed)

# clusterProfiler 4.19.7

+ `interpret()`, `interpret_agent()`, and `interpret_hierarchical()` now use aisdk's global default model when `model = NULL`, so users can switch the package-wide default with `aisdk::set_model()` while still overriding per call with an explicit `model` argument (2026-03-31, Tue)

# clusterProfiler 4.19.6

+ update `ko2name()` to robustly parse KO names via KEGG REST, support vector input with deduplication, and return NA when NAME is missing (2026-02-25, Wed)
+ bug fixed for `plot.interpret` (2026-02-05, Thu)

# clusterProfiler 4.19.5

+ `interpret()` prompt optimized with 'Comparative Analysis' and 'Rule of Exclusion' to better distinguish cell types with shared functions (e.g. NK vs CD8+ T cells) using specific marker genes (2026-01-22, Thu)
+ fixed a bug in `interpret()` where empty interpretation results were returned due to incorrect list structure handling in `process_enrichment_input` (2026-01-22, Thu)
+ `interpret()` now considers specific marker genes in cell type annotation to avoid key markers being overshadowed by general pathways (2026-01-21, Wed)
+ optimize `enrichGO()` to avoid memory boom when keyType is not ENTREZID (2026-01-21, Wed, #805)
+ `gson_GO_local()` to support local GO annotation by adding ancestral terms (2026-01-21, Wed)
+ `interpret()` implements a gene-based fallback mode for clusters with no enriched pathways, ensuring comprehensive analysis (2026-01-20, Tue)
+ `plot()` method for interpretation object to visualize the LLM-inferred regulatory network using `ggtangle` (2026-01-20, Tue)
+ `interpret_agent()` supports multi-agent system (Deep Mode) for interpretation (2026-01-20, Tue)
  - Agent Cleaner: Filters noise and selects relevant pathways
  - Agent Detective: Identifies key regulators and functional modules using PPI/TF data
  - Agent Synthesizer: Synthesizes findings into a coherent narrative
+ `interpret()` supports 'Knowledge-Guided Interpretation' (2026-01-20, Tue)
  - `add_ppi` parameter to integrate PPI network and identify hub genes
  - `gene_fold_change` parameter to incorporate expression levels
  - Mixed-source enrichment analysis support (e.g. Pathways + TFs) for causal integration
  - LLM-guided network refinement to output core regulatory networks
+ `interpret_hierarchical()` for hierarchical interpretation (e.g. Major -> Minor clusters) (2026-01-20, Tue)
+ `interpret()` now supports `prior` parameter for reference-guided interpretation (e.g. from SingleR/scGPT) (2026-01-20, Tue)

# clusterProfiler 4.19.4

+ `interpret()` now supports `task` parameter to specify the task: 'interpretation', 'annotation' and 'phenotyping' (2025-01-18, Sat)
+ `interpret()` supports `enrichResult`, `gseaResult`, `compareClusterResult` and `list` of enrichment results (2025-01-18, Sat)

# clusterProfiler 4.19.3

+ instead of packing KEGG cache data in the package, we now download it from <https://yulab-smu.top/clusterProfiler> (2025-12-15, Mon)
+ add github action to automatically update KEGG cache data (2025-12-09, Tue)
+ use 'enrichit' as engine for enrichment analysis (2025-12-07, Sun)

# clusterProfiler 4.19.2

+ use 'quarto' as vignette engine (2025-11-20, Thu)
+ update KEGG cache data (2025-11-20, Thu, #792)
  - number of KEGG pathway with category information: 582
  - Number of species: 11344

# clusterProfiler 4.19.1

+ bug fixed in `enrichPC` (2025-11-01, Sat, #789) 

# clusterProfiler 4.18.0

+ Bioconductor RELEASE_3_22 (2025-10-30, Thu)

# clusterProfiler 4.16.0

+ Bioconductor RELEASE_3_21 (2025-04-17, Thu)

# clusterProfiler 4.15.2

+ more general regular pattern to remove species suffix in KEGG pathway name (2025-02-27, Thu)
+ remove input duplicated genes in `groupGO()` (2024-11-29, Fri, #741)
 
# clusterProfiler 4.15.1

+ `simplify()` keeps the most informative term if there exist multiple terms that meets the criteria (2024-11-29, Fri, #744)
+ add 'RichFactor', 'FoldEnrichment' and 'zScore' in `enrichDAVID()` result (2024-11-12, Tue)
+ update DAVID Web Service URL to make `enrichDAVID()` work properly (2024-11-09, Sat)
  - <https://davidbioinformatics.nih.gov/content.jsp?file=WS.html>
+ add new citation (2024-11-07, Thu)
  - <https://doi.org/10.1016/j.xinn.2024.100722>

# clusterProfiler 4.14.0

+ Bioconductor RELEASE_3_20 (2024-10-30, Wed)

# clusterProfiler 4.13.4

+ re-export `DOSE::enrichDO()` and `DOSE::gseDO()` (2024-10-01, Tue)

# clusterProfiler 4.13.3

+ fixed bug in `enrichPC()` (2024-08-26, Mon)

# clusterProfiler 4.13.2

+ fixed bug of `gson_KEGG()` (2024-08-19, Mon)

# clusterProfiler 4.13.1

+ update functions to access PathwayCommons data (2024-08-11, Sun, gson#9) 
+ use `yulab.utils::yulab_msg()` for startup message (2024-07-26, Fri)
+ update `kegg_category` information (7 categories and 572 subcategories) (2024-07-26, Fri)
  - Cellular Processes (36)                    
  - Drug Development (75)
  - Environmental Information Processing (41)
  - Genetic Information Processing (39)
  - Human Diseases (99)
  - Metabolism (190)
  - Organismal Systems (92)
  
# clusterProfiler 4.12.0

+ Bioconductor RELEASE_3_19 (2024-05-15, Wed)

# clusterProfiler 4.10.1

+ bug fixed in parsing KEGG category (2024-03-07, Thu, #664)
+ update citation (#656) and wikipedia data URL (2024-01-10, Wed, #633)

# clusterProfiler 4.10.0

+ Bioconductor RELEASE_3_18 (2023-10-25, Wed)

# clusterProfiler 4.9.5

+ fixed R check (2023-10-18, Wed)

# clusterProfiler 4.9.4

+ use `check_installed()` to check package dependency (2023-09-08, Fri, #621)
+ use `yread()` in WikiPathway utilities (2023-09-07, Thu)

# clusterProfiler 4.9.3

+ `enrichKEGG()` and `gseKEGG()` now supports `organism = 'cpd'` to accept KEGG Compound ID (2023-08-31, Thu)
  + `gson_cpd()` and `gson_ko()`
+ use `yulab.utils::yread()` to parse file (2023-08-15, Tue)
+ supports Pathways Common (2023-08-02, Wed, #613)

# clusterProfiler 4.9.2

+ `append_kegg_category()` function to add KEGG pathway category information to KEGG enrichment result and now it is the default behavior of `enrichKEGG()` and `gseKEGG()` (2023-07-12, Wed)
+ parse KEGG Pathway Category information (2023-07-11, Tue)
+ mv `parse_gff()` to `GOSemSim::read.gaf()` and re-export (2023-07-10, Mon)
+ mv `buildGOmap()` to `GOSemSim::buildGOmap() and re-export

# clusterProfiler 4.9.1

+ `getPPI()` to query PPI network from 'stringdb' (2023-05-15, Mon)
+ `getTaxID()` and `getTaxInfo()` functions to query taxonomy information (2023-05-14, Sun)

# clusterProfiler 4.8.0

+ Bioconductor RELEASE_3_17 (2023-05-03, Wed)

# clusterProfiler 4.7.2

+ change wikiPathways link. (2023-03-10, Fri)
+ update `get_data_from_KEGG_db()` for the KEGG api changes (2023-03-05, Sun)
+ removing species info at the end of KEGG pathway names (2023-03-05, Sun)
 
# clusterProfiler 4.7.1

+ update according to the KEGG api changes (2023-03-01, Wed)

# clusterProfiler 4.6.0

+ Bioconductor 3.16 release 

# clusterProfiler 4.5.3

+ `GSEA()` supports `GSONList` object (2022-09-21, Wed)
+ `enricher()` supports `GSONList` object (2022-09-06, Tue)

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
