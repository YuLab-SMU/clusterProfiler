# Mechanism evidence: three-stage completion record

Development version **4.21.5.003**, formal baseline **4.21.5**. Completed 2026-11-11.
This is a correctness, reproducibility and disclosure change, not a new causal
identification algorithm or a validated accuracy improvement.

## Completed changes

1. **Correctness:** replication belongs to one candidate; member support gates
   significance-based promotion; mapping covers every pathway; finite-member and
   known-context counts cannot be inflated by missing or repeated identities.
   Coverage lookup failures remain unknown. Explicit adjudication coverage floors
   are effective, actual evidence changes its fingerprint, shared dataset witnesses
   are deduplicated, and measured-only inputs reach interpretation.
2. **Fixed-denominator audit:** all17 profiles are retained, full rankings distinguish
   unscorable/missing labels from low ranks, and activation/KO plus study/driver
   strata are reported. Frozen inputs, checksums and an independently extracted
   standalone reproduction archive are provided. Earlier denominator/PPV claims
   are explicitly superseded, not silently re-labelled as current results.
3. **Bounded disclosure:** every retained candidate/dataset test is available with
   support, conflict, nonsignificance, failure or non-computability status. Reports
   expose per-candidate summaries, requested-but-unresolved candidates and coverage
   gaps. Nonsignificance is not exclusion; annotation overlap and broad DE are
   descriptive flags. Adjudicated prompts serialize exact identity, directions
   and writable strength outside pathway truncation.

Independent final review found three further edge cases, all fixed and covered:
requested candidates disappeared when no datasets resolved; row subsets imported
unrelated pathway ledgers; blank cell identities falsely established replication.
Conflicting dataset/cell identities are also rejected with a warning. A bounded
live replicate then exposed one more false negative: dataset resolution required
the literal method `"knockdown"`, so a factor recorded only as CRISPRi or knockout
(`STAT1`) read as "no dataset" while the report inventory listed its data. The
resolver now falls back to the other recorded spellings when the literal
`"knockdown"` query is empty (it stops at the first spelling that yields datasets,
so methods are not united), overexpression stays excluded, and the same profile
moved from zero calls to five real non-significant tests.

A second, independent review of both commits found three more defects, all fixed
here: the lookup error above; a row with direction supplied but no `signed` list
claiming its members were missing from that list (reason `no-members`) when the
missing input was the fold changes (now `no-fold-changes`, again with the ledger
retained); and a pure column reorder (`pe[, rev(names(pe))]`) silently stripping
the retained ledger while keeping the class, so `mechanism_report()` saw no
candidate evidence. The review also found a pre-existing order dependence in
`adjudicate()`: with two measured rows for one regulator, the first row alone
decided `perturbation_effect`, so reversing the same bundle flipped activator to
repressor while the fingerprint stayed identical. Every measured row is now read
and disagreement is reported as `mixed` with no single effect. The same review
confirmed that a failed lookup produced
byte-identical output to a successful empty one at the previous commit, which is
the P1 fixed by `dataset_lookup`/`lookup-failed`.

Reviewing that change exposed the mirror-image error it could amplify: the resolver
swallowed every lookup error, so a network outage or a machine without the suggested
`gpsadb` package produced "No usable perturbation data found in GPSAdb" and reason
`no-dataset` -- a false statement about the database. Resolution now reports
`available`, `empty` or `failed` per candidate (retained as the `dataset_lookup`
attribute); failures produce reason `lookup-failed` and an explicit warning, and a
partial failure is disclosed even when another candidate resolves.

## Actual validation

| Gate | Observed result |
|---|---|
| Project documentation | `make rd` completed; Rd generated from R/roxygen, not hand-edited |
| Eight focused test files | **559 passed, 0 failures, 0 errors**; 4 pre-existing network-gated skips (recorded warnings are the expected empty-input warnings) |
| Full project check | `make check`: **0 errors, 0 warnings, 0 notes**, 15m30.4s; tests, examples and vignettes passed |
| Earlier check | A prior tree (4.21.5.001) reported one NOTE, `Unable to verify current time`; the final tree did not |
| Frozen benchmark | `--self-test` passed; complete expected-label ranks and denominators asserted |
| Standalone reproduction | Archive extracted separately; base-R rerun produced byte-identical output manifest `ba060060f4776ff9ba5626a2f5f031aa` |
| Real-data/live-service smoke | Hypoxia profile, HIF1A/EPAS1, two pathways, maximum two datasets per candidate; **8 retained outcomes, 4 test calls**, report/ledger agreement asserted |
| Whitespace | `git diff --check` passed |

The online smoke found consistent HIF1A evidence in two known cell lines for both
related pathways. The report also retained EPAS1 outcomes instead of hiding them
behind the selected HIF1A row. The two pathways share their member set, so these
are not two independent biological confirmations. No live language-model call was
used as a compliance or causal-validation test.

## Fixed exploratory corpus results

| Endpoint | Result |
|---|---:|
| Expected labels scorable | 17/17 |
| Default signed hit@1 / hit@5 / hit@10 | 6/17 / 10/17 / 13/17 |
| Magnitude hit@1 / hit@5 / hit@10 | 4/17 / 8/17 / 10/17 |
| Historical top1 target-overlap set coverage | 7/17, mean size 1.882 |
| Top1-anchored, equal-size random expectation | 6.5/17 |
| Historical candidate verification | 7 expected-label matches among 37 historical passes; not current grades or causal PPV |

The fixed corpus has 16 activation profiles and one KO, 13 expected driver labels,
13 known study groups plus one unknown-study profile. Repeated contrasts and
condition-informed labels do not establish independent or exclusive causal truth.
The module endpoint provides no basis for claiming a substantial validated gain.
No target-count weighting parameter was introduced.

## Bounded current-grading replicate (not a replacement for the frozen table)

Because the frozen verification table was produced by pre-fix code, the expected
driver of every frozen profile was re-graded with the fixed logic at
deliberately capped depth: at most **5 pathways** and **2 datasets per driver**,
`alpha = 0.05`. This run is a pipeline and reason-code check, not a full re-scan.

| Observed | Count |
|---|---:|
| Profiles graded | 15/17 |
| Profiles with no enrichment result | 2/17 (`GSE311462` oestradiol, `GSE276150` cigarette extract) |
| Perturbation test calls / failed | 102 / 0 |
| Expected driver with any supporting test | 8/17 |
| Expected driver with at least one A pathway | 4/17 |
| Datasets where the historical table had an expected-label verdict | 10 |
| Historical A/B, still any support | 4 |
| Historical A/B, no support at this cap | 3 |
| Historical C, any support | 0 |

The losses are explained by the cap as much as by the fix: `best` is chosen among
the datasets actually queried, so asking two instead of six can change both the
`es` and the significant/conflicting split. For example both NR3C1 profiles now
record three or four significant *conflicts* rather than support, and SMAD3 records
five. These rows are not evidence that the old grades were right; they are evidence
that a two-dataset cap cannot adjudicate them. Only a full re-scan under the fixed
logic could, and that has not been run.

Live data availability is not guaranteed. The stored clean run has no failed
calls, but an earlier attempt of this replicate returned no usable statistics for
the hypoxia profile in any of its ten calls; that attempt was not retained, and
the stored artifact is the clean rerun. Back-to-back repeats of the same profile in
one session now reproduce exactly (8 supporting outcomes, A3/B2, minimum padj
7.30e-08; `live_repeat.tsv`), whereas the stored run graded that profile A2/B3
with minimum padj 1.27e-08. The two observations differ in code version as well as
time, so the difference cannot be attributed to the service alone. What the
artifacts establish is narrower: a stored live grade is one sample, not a
reproducibility guarantee, and the historical table must not be replaced by it.

Artifacts: [summary](<../../adjudicate-interpret/results/mechanism-4215-regrade/summary.tsv>),
[retained tests](<../../adjudicate-interpret/results/mechanism-4215-regrade/perturbation-tests.tsv>),
[back-to-back repeat](<../../adjudicate-interpret/results/mechanism-4215-regrade/live_repeat.tsv>),
[script](<../../adjudicate-interpret/scripts/64_regrade_4215.R>).

## Reproduce and inspect

- [Workflow and interpretation boundaries](<mechanism-evidence-workflow.md>)
- [Benchmark methodology/results](<mechanism-benchmark-4215.md>)
- [Benchmark runner](<mechanism-benchmark-4215.R>)
- [Frozen input manifest](<../../adjudicate-interpret/results/audit-benchmark-4215/input-manifest.tsv>)
- [Live per-dataset smoke results](<../../adjudicate-interpret/results/mechanism-4215-live-smoke/perturbation-tests.tsv>)
- [Standalone benchmark archive](<../../adjudicate-interpret/results/mechanism-audit-4215.tar.gz>)

Historical verification was not fully rerun against the fixed grading logic; the
bounded replicate above covers the expected driver at reduced depth. Original DE
was not recomputed. Structured prompt constraints are auditable input restrictions,
not a formal guarantee of model compliance. External independent validation and
intervention-based causal testing remain outside this completion.
