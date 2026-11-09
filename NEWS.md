# clusterProfiler 4.21.5

## interpret()

+ `interpret()` takes `prior`: another tool's answer (a SingleR or scGPT label, a TF-activity ranking), sent after the evidence as a claim to check rather than as evidence. The record adds `conclusion`, `prior_verdict` and `prior_rationale`. On PBMC 10k with SingleR labels, 12 of 14 verdicts on misassigned subtypes were Corrected and 1 of 14 on correct ones (2026-11-02, Mon)
+ `reduce_redundancy()` defaults to 0.7, the default of the `simplify()` it calls, and its documentation now describes what `simplify()` does rather than the earlier greedy selection. Cisplatin GO ORA keeps 81 of 141 terms (0.5: 40), geneList GO ORA 76 of 195; KEGG lists are kept whole (2026-11-02, Mon)

## GO annotation cache

+ The built GO gene-set object is cached per (organism, ont, keytype) on top of the annotation cache: 1.4 s of a 2.8 s warm BP call, and `compareCluster()` pays it once per cluster (13 clusters: 39.5 s against 19.7 s, results identical). Up to three objects are kept, least recently used first out, and an `evidence`-filtered annotation is never cached (2026-10-29, Thu)

## pathway-level mechanism evidence

+ `pathway_evidence()` tests enriched pathways against measured perturbation data instead of hunting for transcription factors. Given candidate drivers it resolves each to its knockdown datasets and reports A/B/C with a `statement` column holding the claim a perturbation experiment could refute. Candidates are supplied by the caller on purpose: picking them automatically carried no usable signal (2026-10-30, Fri)
+ Level A requires replication across **two cell lines**, a direction that agrees with the pathway, and a pathway direction that does not come from the pathway's own members. Before this, a member-derived direction earned B on a tautology (the agreeing fraction is at least 50% by construction) and significance alone earned A even when the perturbation moved the pathway the way it was already running (2026-10-30, Fri)
+ `reason` records why a pathway stopped short of A -- `single-cell-line`, `direction-conflict`, `no-dataset`, `too-small`, `not-significant`, `no-test-completed`, `member-derived`, `low-support` and the input-side cases -- so "nothing found" is distinguishable from "nothing could be looked at". A failed GPSAdb call is `no-test-completed` with a warning, not weak evidence, and a timeout no longer reads as "no perturbation test identified a driver" (2026-10-31, Sat)
+ `mechanism_preflight()` decides before a run whether the data can carry a driver claim, and returns the numbers behind the decision. On eleven real datasets it blocked 7 of 7 that gave no usable answer and 0 of 4 that did, and every reason was visible in advance (2026-11-03, Mon)
+ `mechanism_report()` turns `reason` into "what to do next" and gives the candidate's full inventory, separating data, analysis and "no experiment would change this". It stops recommending experiments when the preflight says attribution is not supported, because adding a cell line cannot rescue a dataset with 23 differentially expressed genes, and it resolves aliases so a report on `GR` finds `NR3C1`'s dataset (2026-11-03, Mon)
+ `mechanism_report()` gains `network`, `min_overlap` and `global_effect_de`. `module` groups drivers whose regulator-target sets overlap, because paralogues, subunits of one complex and parallel pathway arms cannot be told apart by an expression experiment: on 17 benchmark datasets this covers the true driver 7/10 against 6/10 for the single top-1 and lifts module-level PPV from 19% to 54%, recovering HIF1A-EPAS1, STAT1-IRF1, XBP1-ATF6, TP53-TP73 and NR3C1-AR with no prior annotation. `global_effect_reason` marks a driver whose knockdown moved more than `global_effect_de` genes (default 3000), since such a perturbation drags any pathway of this size and a supported verdict on it says "large perturbation", not "specific driver". Neither field changes the A/B/C class (2026-11-08, Sat)
+ `mechanism_reproduction()` reports how many of a driver's independent perturbation datasets move the pathway the same way, and `pathway_evidence(reproduction = )` uses it to word a single-cell-line row. The majority sign is the mode of the per-dataset signs rather than `sign(sum(shift))`, so one dataset with a large effect cannot overturn two smaller ones that agree, and datasets whose members are all NA are dropped before `k` is counted (2026-11-07, Sat)
+ `mechanism_discover()` ranks candidate transcription factors for a profile, so `pathway_evidence()` gets a shortlist instead of searching 265 candidates at 1-2 s per query. On 17 experiments the true driver ranked first in 43%, the top five in 71% and the top ten in 93% (median 2) among the 14 drivers a typical network contains; the other three are absent from CollecTRI. The shortlist is not a verdict: verifying the top five let 70% of true drivers through but also 40% of the decoys, a positive predictive value near 19%, and the decoys that pass are the driver's neighbours -- ATF4 and ATF6 beside ER stress, EPAS1 beside HIF1A, IRF1 downstream of STAT1. The score also rises with the number of annotated targets, even among candidates that are not the driver; `n_targets` and `share` report this rather than correcting it (2026-11-07, Sat)
+ The attribution route explored on branch `mechanism-decompose` is not merged. Five statistics over one reference matrix landed on the same recall/precision trade-off, the best pointing at the true driver in 6 of 9 datasets against 6 of 9 for this implementation, with two failure modes no reweighting fixed. `dev_docs/mechanism-attribution-postmortem.md` records the five things not to try again (2026-11-03, Mon)
+ Fixes: a level-B row now says a single-cell-line limit is the data rather than the evidence, while a direction conflict keeps its own wording; a member-derived direction no longer reads "100% of members moved up, below the support floor of 50%"; `reason` reports `too-small` ahead of `low-support`; `coherent-refuted` is separated from `coherent-untested`; an empty result, an outage and an unreachable level A each warn instead of passing as findings; `print()` runs for `pathway_evidence` and a column subset no longer breaks it; empty and `NA` tokens are not counted as members; alias lookup runs in both directions, so `drivers = "GR"` finds `NR3C1` (2026-10-31, Sat)

## mechanism analysis: collect, adjudicate, write

+ `evidence()` grades a gene list against upstream collections and returns an `evidenceBundle` whose `$regulators` carries `evidence_kind` (`measured`, `predicted`, `observed` or `none`) and a level. `supported` requires one measured source corroborated by a second: counting sources is not counting evidence when two predictions agree (2026-10-25, Sun)
+ `adjudicate()`, formerly `consensus()`, decides which claims survive a change of evidence type and returns `verdicts`, the evidence behind them, `diagnostics`, `parameters` and a provenance fingerprint; `interpret()` requires an `adjudication`, so that decision cannot be skipped on the way to prose (2026-10-28, Wed)
+ `upstream = "perturbagen"` asks GPSAdb which measured perturbation moves the input genes -- the only source that reaches level A. `upstream = "tfunion"` adds ENCODE ChIP-seq, ChEA 2022, TF_Perturbations and TRRUST: 2869 sets naming 315 factors C3:TFT has no set for, including EPAS1, NFKB1, TP53, CTCF, ESR1, JUN and FOS (2026-10-27, Tue)
+ `mechanism_writable()`, `mechanism_xref()`, `collapse_findings()`, `mechanism_graph()` with `plot_mechanism()`, `mechanism_sets()`, `consensus_findings()`, `mechanism_significant()` and `mechanism_alias_regulator()` are exported, and `simplify()` removes redundancy from any collection rather than only GO (#743, 2026-10-26, Mon)
+ Fixes: `perturbation_effect` no longer calls an activator a repressor when the caller's own contrast knocks that factor down, and `measured_direction` now reads the same gene set it does; `unsupported_claims()` no longer drops real symbols such as `P53` and `P300`, no longer exempts a rejected regulator because its symbol was also an input gene, and no longer excuses one whose symbol prefixes an unrelated gene; `print()` reports coverage and a failed source instead of silence, and no longer counts aliases against candidate regulators (2026-10-29, Thu)

# clusterProfiler 4.21.4

## fixes

+ The headline verdict could be produced by a symbol collision. TRANSFAC's `NF1` is Nuclear Factor I while the HGNC symbol `NF1` is neurofibromin, and matched on the string the motif term and a knockout dataset "agreed". A collided label can no longer reach `supported` (2026-10-25, Sun)
+ A protein name that maps to a gene is no longer treated as a family label. TRANSFAC names nuclear receptors by their protein -- `GR`, `PXR`, `FXR`, `HIF2A` -- none of which are HGNC symbols, so the perturbation channel skipped them as families and `GR` could not reach `supported` on a glucocorticoid experiment, where it is the one answer that is certain (2026-10-26, Mon)
+ Regulators are looked up in perturbation databases under every name they are known by, in both directions: the collections say `GR` where the database says `NR3C1`, and `HIF1` where it says `HIF1A`. Querying one spelling returned nothing, which is indistinguishable from "never perturbed" (2026-10-27, Tue)
+ `mechanism_is_family_label()` no longer reads a hyphenated symbol as a family: TRANSFAC writes `NKX2_5` where HGNC writes `NKX2-5`, which also stripped those genes of their measured evidence kind (2026-10-25, Sun)
+ The perturbation sweep no longer runs concurrently: GPSAdb returns different data when queried in parallel and reports nothing unusual, so the supported set flipped between runs. `options(clusterProfiler_perturb_workers = 4)` takes the trade knowingly (2026-10-26, Mon)
+ A p-value tie between two datasets for one regulator is settled by system match, then p-value, then genes moved, then dataset id -- the last unique, so the ordering is total rather than dependent on arrival order (2026-10-26, Mon)
+ Corroboration now requires the sources to cite different datasets, not merely different source names: `perturbagen` and `biotrigger` are two views of one GPSAdb experiment, so when both cited the same dataset the server was agreeing with itself (2026-10-26, Mon)
+ `evidence_similarity()` and `evidence_biotrigger()` are exported but no longer part of the bundle, after measuring what they contribute: `biotrigger` raised the `supported` count without adding a transcription factor, and the additions were genes with no bearing on the system under study (2026-10-26, Mon)
+ microRNA collections are no longer tested by default. Predicted seed matches gave 7 significant miRNAs among the up-regulated genes and none among the down ones, the opposite of what a repressor predicts, so the direction that would make the result interpretable is absent (2026-10-25, Sun)
+ `upstream =` also accepts `"cellstate"` (C8) and `"tissue"` (C4:CGN), whose set identifiers already name the entity and are passed through unchanged rather than through the transcription-factor rules, which would turn a cardiac signature into the gene HPX (2026-10-25, Sun)
+ `network = TRUE` is documented as informational: `$ppi` carries the edges but nothing consumes them into a verdict, so the channel is reference material rather than corroboration (2026-10-28, Wed)
+ `evidence_ppi()` grades STRING edges by the channel behind them: `escore` or `dscore` at 0.4 or above gives B, anything else C, and no STRING channel measures a perturbation so an edge never reaches A. `pscore` is deliberately unused -- it is 0 even for NFKB1-NFKBIA, which carries `escore = 0.999` (2026-10-25, Sun)
+ The session cache directory is no longer forced to `tempdir()` on attach, which overwrote a directory set in `.Rprofile` and made downloaded collections miss their cache (2026-10-25, Sun)
+ **Switching ontology no longer evicts the GO annotation cache.** The cache was keyed on the ontology and stored only the last call's annotation, so BP then MF then BP rebuilt it every time: five alternating `enrichGO()` calls took 85.5 s of which the enrichment itself is 0.4 s, and now take 11.6 s, with identical results (2026-10-27, Tue)
+ **A `gseaResult` is now tested with GSEA, not re-tested with ORA.** Taking the union of `core_enrichment` and running `enricher()` on it discarded the ranking and the sign, so the "GSEA path" was the ORA path with a different input and the two agreed against ground truth because they were the same test (2026-10-27, Tue)
+ The two front ends can legitimately disagree: on the airway list `GR_Q6` has ORA padj 0.0002 and GSEA padj 0.166, both correct, because one asks about over-representation among significant genes and the other about moving monotonically with the phenotype. Use ORA for the first question, GSEA for the second (2026-10-27, Tue)
+ An interrupted Enrichr download no longer poisons the cache: the union libraries count as cached only when a `.meta` marker sits beside the GMT, as the MSigDB files already did. A truncated file used to be reused, shrinking the collection and grading every regulator outside it `absent` instead of `failed` (2026-10-28, Wed)
+ `tft` combined with `upstream = "tfunion"` is an error rather than a silent no-op. `tft = "gtrd"` is how a caller escapes the family and consensus labels no perturbation dataset can answer for, and dropping it turned a 504-set search into a 2869-set one without a word (2026-10-28, Wed)
+ `enrichDAVID()` applies `qvalueCutoff` row by row: an `NA` q-value means Storey estimation failed for that row, so it bypasses this cutoff instead of switching the cutoff off for every term (clusterProfiler#710, #584, #614, #332, 2026-10-27, Tue)
+ `gseKEGG()` reports the same `category` / `subcategory` columns as `enrichKEGG()`; a reader that filtered on `category` silently behaved differently depending on which function produced the result (#787, 2026-10-26, Mon)

## interpret()

+ `$pathways` carries `pathway_padj`, the pathway's own enrichment p-value, beside `padj_within_pathway`, which is the attribution's. The two answer different questions and without the first the narrative context could only order by attribution strength: on the hypoxia contrast that put "carbohydrate metabolic process" above "response to hypoxia" (2026-10-28, Wed)
+ An adjudicated run now adds one fixed sentence to its system prompt: only `supported` findings may be stated as mechanisms, `single_source` is a hypothesis, `absent` is not "does not exist", and `untested` is not an absence. It is appended rather than substituted, so the four legacy input types keep the prompt they had byte for byte -- a test compares the strings rather than trusting it (2026-10-28, Wed)
+ The pathway layer now reaches the narrative. `interpret(adjudication)` lists which regulator explains each enriched pathway, in its own prompt section marked as context, and only for entities the adjudication already graded -- so an attribution can never introduce a regulator or make a claim sound stronger. The "context, not evidence" instruction travels inside that section rather than in the system prompt, which would have changed the bytes of the four legacy input types (2026-10-28, Wed)
+ `network = TRUE` is documented as reference material: `$ppi` carries the graded STRING edges but nothing turns an edge into a corroboration, so it changes no verdict. A test pins that, so the assumption cannot creep back (2026-10-28, Wed)
+ fix `interpret()` reverting a caller's adjudication. It re-ran `adjudicate()` with the default `min_independent` in three places -- the prompt builder and the two post-hoc evidence tables -- so a caller who adjudicated at a stricter setting had it silently overridden, which is the one thing making adjudication a separate step exists to prevent. The carried verdict table now travels through the whole call (2026-10-28, Wed)
+ `interpret()` gains `evidence_status()`, a deterministic grounding check: every gene symbol the interpretation names is classified as `observed` or `out_of_evidence`. No model is called, so unlike a second model the check is exact (2026-09-26, Sat)
+ `interpret()`, `interpret_agent()` and `interpret_hierarchical()` stamp their results with an execution record, read back with `provenance()`: the model actually addressed, the decoding parameters and whether the response confirmed them, digests of the prompts and response, token usage and `finish_reason`. A remote model is not a versionable dependency, and without the record a transport failure arrives as an ordinary result (2026-09-03, Thu)
+ `interpret()` reads per-gene statistics through `gene_stats`, a differential-expression table whose gene, log fold change and adjusted p-value columns are recognised by name or row names, and lists the genes of the enriched terms with those values (2026-10-18, Sun)
+ `interpret()` no longer truncates the genes it sends: every input gene is listed, every fold change of the sent terms' genes is given, and an agent cluster without enriched terms lists every input gene. `n_pathways` stays `NULL`, because on five public lists a lower cap named fewer themes; lower `reduce_threshold` instead (2026-10-11, Sun)
+ `interpret()` now interprets an enrichment result and nothing else, so its `task` and `prior` arguments are removed; combining enrichment evidence with other tools moves to `interpret_agent(prior = )`, whose Synthesizer returns `conclusion` and `prior_verdict` (2026-09-30, Wed)
+ `reduce_redundancy()` and `reduce_go_redundancy()` are exported and `interpret()` takes `reduce` / `reduce_threshold`, so every significant theme reaches the model once: GO by Wang semantic similarity, other collections by Jaccard overlap (2026-09-20, Sun)
+ `interpret()` and `interpret_agent()` take `reasoning`, passed through to the provider and recorded in the execution stamp; left `NULL` the request is byte-identical to one made before the argument existed (2026-09-20, Sun)
+ `interpret_agent()` gains an optional Reviewer agent (`review = TRUE`) that audits the draft, reporting grounding separately from it (2026-09-03, Thu)
+ `interpret()` and `interpret_agent()` warn when `add_ppi = TRUE` and the STRING query fails, instead of silently interpreting without PPI evidence (2026-09-30, Wed)
+ fix `interpret()` to read gene names from the object's symbols rather than from a mismatched identifier column; an unparsed model response is reported instead of returning `NULL` (2026-09-26, Sat)
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
