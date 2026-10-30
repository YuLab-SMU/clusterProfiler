# clusterProfiler 4.21.4.038

## GO annotation cache

+ **The built GO gene-set object is cached per (organism, ont, keytype) instead of being rebuilt on every call.** After the annotation cache landed, what a warm call still paid was reshaping that annotation into the object the engine consumes: 1.4 s of a 2.8 s warm BP call on human ENTREZID, and `compareCluster()` pays it once per cluster, so 13 clusters took 39.5 s against 19.7 s with the object cached, results identical. Up to three objects are kept, least recently used first out, so BP/MF/CC of one organism cost about 35 MB and a session walking through organisms evicts rather than accumulates; an `evidence`-filtered annotation is never cached, the same rule the annotation cache already applied. This is the second level of the two-level cache proposed in clusterProfiler#71, the first level being the annotation cache that no longer evicts on an ontology switch (2026-10-29, Thu)

## pathway-level mechanism evidence

+ `pathway_evidence()` detects a numeric member space (ReactomePA returns Entrez) and translates to symbols itself, so a Reactome enrichment result can be passed straight in. Verified: "Metabolism of carbohydrates and carbohydrate derivatives" and `Glycolysis` reach level A on HIF1A from `enrichPathway` output with no manual conversion (2026-10-30, Fri)
+ `pathway_evidence()` refuses to run when no pathway member matches the `signed` gene list, instead of returning every pathway as untestable. Passing a ReactomePA result alongside symbol fold changes hit exactly this: Entrez members matched nothing, no network call was made, and the output looked like a normal empty finding. The message names the likely cause. A partial match is not an error -- a few members missing from the differential table is ordinary (2026-10-30, Fri)
+ Reactome enrichment output can now be used directly once its `geneID` column is translated to symbols, which removes the GO-to-Reactome mapping step and its 8% coverage. Verified on the hypoxia pilot: "Metabolism of carbohydrates and carbohydrate derivatives" reaches level A on HIF1A at padj 1.5e-08 (2026-10-30, Fri)
+ `pathway_evidence(reactome = )` takes a second annotation source instead of the PPI comparison. `reactome.db` also carries an official GO cross-reference, and for `GO:0071456` that maps to "Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha" -- an exact mechanistic match where PPI density had nothing to say. Overlap scoring is reported in `reactome_overlap` with a three-way label; note that the member-overlap fallback picks IL-4/IL-13 signalling for every hypoxia pathway, because that pathway shares HIF1A and glycolytic genes, so the overlap alone is not corroboration and the threshold reflects that (2026-10-30, Fri)
+ `pathway_evidence(ppi = TRUE)` compares a pathway's PPI density against size-matched random sets as a second evidence type, reported as `agrees` / `no-additional-signal` / `disagrees` and never used to change a grade. The three-way label matters: a tie is not a contradiction, and calling one "disagrees" would read as the sources opposing each other when the second source simply has nothing to say. On the hypoxia pilot all six pathways come back `no-additional-signal` -- GO terms are not protein complexes, so their members are no more connected than random (2026-10-30, Fri)
+ `pathway_evidence()` no longer re-tests datasets for pathways that share the same member set. The hypoxia result has three pathways over one identical set of 33 genes and three more over another, so the same query was being sent repeatedly; results are now memoised on (dataset, member signature) and the call count is exposed as the `network_calls` attribute (2026-10-30, Fri)
+ `pathway_evidence()` now requires **replication across two cell lines** before granting level A. It tests every dataset available for each candidate driver rather than keeping the single most significant one, and reports `n_datasets_tested`, `n_cell_lines_supported` and a `single-dataset` verdict. Evidence from one cell line is not weak evidence about one cell line -- it is unevaluated evidence about the others, and on the hypoxia pilot EPAS1 support existed in a single line while HIF1A support held in five (2026-10-30, Fri)
+ fix `pathway_evidence()` granting level B on a tautology: when the pathway's direction is derived from its own members, the fraction moving together is at least 50% by construction, so `min_support = 0.5` was always met and every pathway with a direction earned a B. Member support is an independent signal only when the direction comes from elsewhere (an explicit `pathway_direction` or a GSEA `NES`); otherwise it is reported but cannot promote (2026-10-30, Fri)
+ fix `pathway_evidence()` grading level A on significance alone. A perturbation that moves the pathway the same way it is already running reads as repression, not drive; it is now reported as `direction-conflict` and held at B (2026-10-30, Fri)
+ `pathway_evidence()` tests enriched pathways against measured perturbation data instead of hunting for transcription factors. Give it candidate driver genes and it resolves each to its knockdown datasets, asks whether that perturbation moves the pathway's members, and reports A/B/C with a `statement` column holding the claim a perturbation experiment could refute. Candidates are supplied by the caller on purpose: stage zero measured what happens when the tool tries to pick a driver itself and the ordering carried no usable signal (2026-10-30, Fri)

## mechanism analysis: collect, adjudicate, write

+ `print()` now says how much of the input the measured evidence covers. A perturbation dataset's `padj_within_source` is the dataset's own best p-value, not evidence that it explains the caller's list, so a level-A `supported` can rest on a few percent of the genes -- the airway example graded GR at A on 23 of 896, against a median of 17% across the datasets queried. A `coverage:` line on the bundle and a `thin:` line on the adjudication make that visible without changing any verdict (2026-10-29, Thu)
+ `options(clusterProfiler_coverage_note = )` sets the threshold, default 0.05; 0 switches the report off (2026-10-29, Thu)
+ fix `unsupported_claims()` treating every single-letter-plus-digits token as a panel label. `P53`, `P63`, `P68`, `P300` and `F10` are real symbols and were dropped with `A4` and `Q4`, so a narrative naming p53 as the mechanism went unseen. Only a letter and one digit is suppressed now (2026-10-29, Thu)
+ fix `unsupported_claims()` letting a rejected regulator through whenever its symbol was also an input gene. Being an input gene exempted a symbol, so "ESR1 drives this program" passed untouched if ESR1 was differentially expressed -- which regulators routinely are, ESR1 at lfc -1.94 in the MCF-7 example. A symbol the adjudication graded and did not admit now keeps no exemption (2026-10-29, Thu)
+ fix the same check's family exemption using the whole vocabulary instead of the admitted labels, so a rejected regulator was excused whenever its symbol prefixed an unrelated gene: MYB was exempted because MYBL1 happened to be an input gene (2026-10-29, Thu)
+ fix `print()` reporting the perturbation sweep as "31 of 30 candidate regulator(s) queried". The numerator counted the names the lookup sent -- a regulator is queried under every alias it is known by -- and the denominator counted candidate regulators, so the two halves of the fraction measured different things. Both are now counted in regulators (2026-10-29, Thu)
+ fix `perturbation_effect` reporting an activator as a repressor when the caller's own contrast is a knockdown of that same factor. The expectation was flipped for a programme that runs down, but the caller has already said which genes to look at, so the question is only whether removing the factor moves those genes down -- true for an activator either way. Measured on MCF-7 cells with ESR1 (an activator) knocked out against GPSAdb's own ESR1 knockdown, which moved the caller's down-genes down 39 of 39: a reproduced experiment, called a repressor (2026-10-28, Wed)
+ fix `measured_direction` reporting the up half while `perturbation_effect` read the half the programme actually runs in, so on a down programme the two columns described different sets of genes. Both now come from one helper (2026-10-28, Wed)
+ `consensus_findings()` and `mechanism_significant()` are exported. Both answer questions a caller asks immediately after `adjudicate()` -- which findings may be discussed, and which rows reached significance -- and both are easy to get subtly wrong by hand: filtering on a p-value ignores a source's own `significant` flag, and treating `absent` as fit to discuss contradicts the verdict (2026-10-28, Wed)
+ `mechanism_alias_regulator()` is exported. The collections say `HIF1` where the database says `HIF1A`, and `GR` where it says `NR3C1`; a caller comparing either side against an outside database hits the same mismatch, and a hand-written alias table drifts from this one (2026-10-28, Wed)
+ A `<enrichResult>` handed to `pathway =` now gets a `pathway_direction` inferred from the caller's own up and down halves, instead of reaching the narrative as "direction not stated". A pathway split evenly between the halves stays `NA` rather than being given the majority side (2026-10-28, Wed)
+ `evidence()` grades a gene list against upstream collections -- transcription-factor targets, microRNA targets, cell-state and tissue signatures -- and returns an `evidenceBundle` whose `$regulators` table carries a per-regulator level instead of a bare p-value. The column is `padj_within_source`: C3:TFT and C3:MIR hold different numbers of sets, so pooling would let a large collection dilute a small one (2026-10-25, Sun)
+ `$diagnostics` keeps separate counts whose denominators differ: `n_input` is what was handed in, `n_tested` is how many could enter the test at all, `n_sig_covered` is how many a significant regulator covers. `$diagnostics$per_gene` labels each gene `supported`, `absent` or `untested`, because the absence of evidence is not the evidence of absence (2026-10-25, Sun)
+ `$regulators` carries `evidence_kind` (`measured`, `predicted`, `observed` or `none`), and `supported` requires at least one measured source corroborated by a second. Counting sources is not counting evidence: two predictions agreeing is one observation confirmed by the same reasoning twice. On the airway list `tf+mirna` gave 177 findings and no `supported`, `tf+perturbagen` 22 and 22 (2026-10-25, Sun)
+ A failed query is no longer reported as an empty answer. A source that fails is recorded in `$diagnostics$failures` and shown by `print()`, and a run in which every consulted source fails aborts rather than reporting `absent` for everything -- "no dataset moves these genes" and "the database rejected the request" lead to opposite conclusions (2026-10-25, Sun)
+ `upstream = "perturbagen"` asks GPSAdb, over its MCP endpoint, which measured perturbation moves the input genes -- the only source that reaches level A, because a knockdown is an experiment rather than an annotation. Candidates are the regulators whose target sets were significant, capped at 30 (`options(clusterProfiler_perturb_candidates =)`); querying every named regulator of a 900-gene list never finished (2026-10-25, Sun)
+ `upstream = "tfunion"` adds ENCODE ChIP-seq, ChEA 2022, TF_Perturbations and TRRUST to C3:TFT: 2869 sets naming 315 factors C3:TFT has no set for, including `EPAS1`, `NFKB1`, `TP53`, `CTCF`, `ESR1`, `JUN` and `FOS`. Pair it with a measured source -- the union is one source however many sets it holds, so on its own it reaches only `single_source` (2026-10-27, Tue)
+ `evidence()` takes a differential-expression table as its canonical input and reads the direction from `design$logfc_col`; a named numeric vector does the same. Each half is tested separately and both are corrected together, so a regulator's row says whether its programme rose or fell. A table without that column warns rather than running undirected without saying so (2026-10-28, Wed)
+ An `enrichResult` is no longer accepted as `genes`: reading genes back out of one drops every gene that appears only in a pathway below its p-value cutoff (897 of them came back as 483 on the airway list). Pass the DE table as `genes` and the result as `pathway =` (2026-10-28, Wed)
+ New `pathway` layer, `pathway = TRUE` or a supplied result: the conventional enrichment is run on each half, and each pathway is then broken down into the regulators whose target sets explain it. Its rows live in `$pathways` and never in `$regulators`, because they derive from the same collections and counting them as a source would let one collection corroborate itself (2026-10-28, Wed)
+ The perturbation layer decides direction with a sign test restricted to one half of the caller's own contrast, instead of requiring every moved gene to agree -- a requirement a knockdown can never meet, which left `measured_direction` at `mixed` for every dataset (2026-10-28, Wed)
+ `adjudicate()` reports `perturbation_effect` (`activator`, `repressor` or `NA`) from the collection's direction together with the perturbation's. It is a report, not a verdict condition: "a knockdown pushes the up-genes further up" means the factor represses them, which is as much a mechanism as activating them (2026-10-28, Wed)
+ `consensus()` is now `adjudicate()`. It returns an `adjudication` carrying `verdicts`, the evidence it was derived from, `diagnostics`, `parameters` and a `provenance` fingerprint. It decides which claims survive a change of evidence type, so "consensus" -- which means agreement between clusterings -- named the wrong operation (2026-10-28, Wed)
+ `interpret()` now requires an `adjudication` and refuses an `evidenceBundle`, so the decision about which claims survive cannot be skipped on the way to prose (2026-10-28, Wed)
+ `$verdicts` separates `supported` (at least `min_independent` sources agree on the same entity and one is graded B or above) from `single_source`, `absent` and `untested`. Agreement counts only for an identical regulator label, and a label that names different things in different sources can never be `supported` (2026-10-25, Sun)
+ `mechanism_writable()` turns a verdict and a level into the sentence the evidence permits, so a level-C finding is offered wording that describes a correlation and never wording that asserts a mechanism (2026-10-25, Sun)
+ `mechanism_sets()` returns an MSigDB collection as a `GSON` for `enricher(gson =)` or `GSEA(gson =)`: `"tf"`, `"tfunion"`, `"mirna"`, `"cellstate"`, `"tissue"` or `"perturbagen"`. The GMT is fetched once and reused only when its recorded release matches, so changing the release cannot silently reuse an old collection (2026-10-25, Sun)
+ `mechanism_xref()` publishes the level rules: A for a measured perturbation, B for a recorded mechanism such as a ChIP-supported target set, C for an association. An unnamed motif set is reported as C whatever rule matched, because a consensus sequence must not be quotable as a mechanism (2026-10-25, Sun)
+ `collapse_findings()` groups findings whose gene sets largely overlap and keeps the strongest of each group. It is a separate step because the need is not uniform: on the airway list transcription-factor findings barely overlap and only tissue signatures are genuinely redundant (2026-10-25, Sun)
+ `mechanism_graph()` turns an `evidenceBundle` into a `mechgraph` object and `plot_mechanism()` draws it, colouring regulators by the grade they earned and naming the collection behind each edge. A link means "this collection lists this gene", not "this factor regulates it", and the caption says so (2026-10-26, Mon)
+ `options(clusterProfiler.progress =)` reports what a slow run is doing: `"major"` (the default) prints one line per source, `"all"` also times each step, `"none"` is silent (2026-10-26, Mon)
+ `simplify()` removes redundancy from any collection, not only GO: `measure = "JC"` compares term gene sets by Jaccard index, `measure = "overlap"` by shared genes over the smaller set. Both read `@geneSets`, so KEGG, Reactome and MSigDB results work without `GOSemSimDATA` (#743, 2026-10-26, Mon)
+ Collections are parsed once per session and the entry is keyed on each source file's path, mtime and size, so a re-downloaded library invalidates it instead of serving a stale one (2026-10-28, Wed)

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
