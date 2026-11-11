# Mechanism evidence: three-stage completion record

Development version **4.21.5.004**, formal baseline **4.21.5**. Updated 2026-11-11.
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
repressor while the fingerprint stayed identical. Every measured row is now read;
raw response directions and method-normalized effects are summarized separately.
Opposing normalized effects yield no single effect. The .004 follow-up retains
explicit `mixed` measurements, does not choose a programme half when collection
directions disagree, and returns unknown for incomplete half-statistics instead of
crashing. Opposite methods can have opposite raw responses yet imply the same
known effect; the interpretation instructions now preserve that distinction.
The same review
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

## Current .004 validation gate

| Gate | Observed result |
|---|---|
| Documentation | `make rd` completed from R/roxygen |
| Eight focused test files | **663 passed, 0 failures, 0 errors**; 2 expected empty-input warning records and 4 existing network-gated skips |
| Full project check | `make check`: **0 errors, 0 warnings, 0 notes**, 16m2.9s |
| Independent scoped review | Mixed measurements, mixed collection directions, missing statistics, row order and connected interpretation paths verified; no remaining defect found in the reviewed changes |
| Frozen benchmark | `--self-test` passed; signed hit@1/5/10 remains **6/17, 10/17, 13/17** |
| Whitespace | `git diff --check` passed |

This gate checks implementation and reporting, not improved discovery accuracy.
No new live grading or model-compliance evaluation was run for .004; the stored
.002 regrade below remains unchanged and explicitly historical.

## Historical .003 validation gate

The following gate records development **4.21.5.003**. It is retained as history
and must not be substituted for the separate .004 gate above.

| Gate | Observed result |
|---|---|
| Project documentation | `make rd` completed; Rd generated from R/roxygen, not hand-edited |
| Eight focused test files | **559 passed, 0 failures, 0 errors**; 4 pre-existing network-gated skips (recorded warnings are the expected empty-input warnings) |
| Full project check | `make check`: **0 errors, 0 warnings, 0 notes**, 15m30.4s; tests, examples and vignettes passed |
| Earlier check | A prior tree (4.21.5.001) reported one NOTE, `Unable to verify current time`; the .003 gate did not |
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
| Historical candidate verification | 7 expected-label matches among 37 historical A/B classifications; not current grades, validated passes or causal PPV |

The fixed corpus has 16 activation profiles and one KO, 13 expected driver labels,
13 known study groups plus one unknown-study profile. Repeated contrasts and
condition-informed labels do not establish independent or exclusive causal truth.
The module endpoint provides no basis for claiming a substantial validated gain.
No target-count weighting parameter was introduced.

## Bounded historical .002 regrade (not a replacement for the frozen table)

The stored [regrade object](<../../adjudicate-interpret/results/mechanism-4215-regrade/regrade.rds>)
records version **4.21.5.2** and timestamp **2026-11-11 10:02:39.08213**:
this is a historical **.002** sample, not a rerun of HEAD .003 or .004.
It re-graded the expected driver of every frozen profile at at most **5 pathways**
and **2 datasets per driver**, `alpha = 0.05`. The historical
[script 60](<../../adjudicate-interpret/scripts/60_discover_then_verify.R>)
used **8 pathways and 4 datasets per driver**, not six datasets, and tested only
the discovered top-five candidates. This bounded regrade is a pipeline and
reason-code check, not a full re-scan or a controlled comparison of code fixes.

| Observed in the stored .002 sample | Count |
|---|---:|
| Profiles graded | 15/17 |
| Profiles with at least one finite p-adjusted value | 12/17 |
| Graded profiles with only non-computable outcomes | 3/17 (ESR1 KO, TP53, VDR) |
| Profiles with no enrichment result | 2/17 (`GSE311462` oestradiol, `GSE276150` cigarette extract) |
| Retained ledger outcomes | 135 |
| Supporting / conflicting / non-significant / non-computable outcomes | 40 / 12 / 42 / 41 |
| Computable outcomes | 94 |
| Network perturbation test calls / failed outcomes | 102 / 0 |
| Expected driver with any supporting test | 8/17 |
| Expected driver with at least one A pathway | 4/17 |
| Datasets where the historical table had an expected-label verdict | 10 |
| Historical expected-driver A/B, retaining some A/B | 7/7 |
| Historical expected-driver A/B, any support | 4/7 |
| Historical expected-driver A/B, zero support | 3/7 |
| Historical C, any support | 0 |

Of ten historical expected-label cases, nine are graded in this stored sample;
the remaining AHR case has no enrichment. The historical C group therefore has
two graded cases and one ungraded case, not three measured negatives.

Ledger outcomes are pathway/dataset records, not network calls: shared queries
can serve multiple outcomes. The summary's `n_tests` does not count all retained
outcomes; inspect the ledger, including its 41 non-computable records.

All seven historical expected-driver A/B cases retain some A/B classification,
but only four have any supporting outcome. The other three have zero support:
NR3C1 in airway has **3 conflicts**, NR3C1 in GSE286890 has **4**, and SMAD3 has
**5**. B includes coherent-only and conflicting evidence; it is not a pass or
validation. Inspect the verdict, reason and retained ledger rather than treating
A/B as a support endpoint. These conflicts are evidence in the queried scope,
not proof that the candidate has no biological role or that the old grade was right.

Four additional supported expected drivers were outside the historical top five
(NFE2L2, STAT3, HEK293T XBP1 and RELA). They represent **added candidate coverage**,
not rescue of historical C grades. The different candidate scope, pathway/dataset
caps, code and live observations preclude causal attribution of changes to the cap
versus fixes. A controlled comparison has not been run; neither the cap nor the
fixes have been established as the explanation. A full fixed-logic re-scan also
has not been run and would not by itself isolate those effects.

Live data availability is not guaranteed. The stored clean run has zero failed
outcomes, but an earlier, unretained hypoxia attempt reported **10 failed ledger
outcomes from 4 network test calls**, not ten calls. The stored artifact is the
clean rerun. The [back-to-back repeat table](<../../adjudicate-interpret/results/mechanism-4215-regrade/live_repeat.tsv>)
shows **aggregate agreement** in two runs (4 calls, 8 supporting outcomes, A3/B2,
minimum padj 7.30e-08). It contains no full per-test results and therefore does not
establish identical ledgers, statistics or full reproducibility. The stored .002
sample instead has A2/B3 and minimum padj 1.27e-08 for that profile. Code version
and time differ as well, so the difference cannot be attributed to the service
alone. A stored live grade is one sample, not a reproducibility guarantee, and
must not replace the historical table.

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
