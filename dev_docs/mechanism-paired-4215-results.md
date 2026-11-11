# Evidence-first output and frozen paired comparison (.005)

## What changed in the interface

Pathway printing and mechanism reports now put **support / conflict / mixed /
undetermined / noncomputable** before the legacy grade. These states summarize
all retained candidate/dataset tests for each pathway; `mixed` means both support
and conflict are present. Candidate-specific assessments remain separately
available. Availability, incomplete coverage and attempted/computable counts are
not collapsed into that state. Missing ledgers cannot establish support, and
noncomputable results are not nonsignificant results.

A/B/C, `class`, `verdict` and the grading/ranking algorithms remain unchanged.
Ledger support can coexist with a failed member-support floor; it is not a grading
pass or a causal conclusion. Serialized .004 reports still print their saved rows
in an explicit legacy section without inventing unavailable evidence states.

## Design and frozen observations

The [design](<mechanism-paired-4215-plan.md>) fixes all **17 condition-declared
expected drivers**, preserving two unenriched inputs in every denominator. This
is not a top-five precision or causal accuracy experiment. Original DE and
enrichment are reused. The four conditions cross:

- historical source `978eabaee502350e12e5b6a3988c54091ac79fe3`;
- fixed source `fe7c060cc26a94178bb512433a019b092c00b00b`;
- caps **8 pathways / 4 datasets** and **5 pathways / 2 datasets**.

Both implementations receive the same first retained response for each exact
request within one snapshot. Source changes include dataset lookup/selection
policy, not only final classification. Identical request keys retain original
argument order and multiplicity. There is no network access in replay; missing,
corrupt or malformed cache entries fail outside either grader's error handler.

### Collection quality and operational amendment

The preserved pilot `mechanism-paired-4215` issued 235 unique requests and retained
85 HTTP 429 failures (83 test requests, 2 lookups). It is not used as biological
negative evidence. Its original manifest did not bind the complete design plan,
and its captured runner predates the strengthened malformed-cache guard. No pilot
response was overwritten or selectively substituted into the final snapshot.

The separate final
[paced snapshot](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/inputs/design.rds>)
recollected **all requests afresh**, with at least two seconds between HTTP attempt
starts, including initialization. The operational amendment is explicit in its
checksummed plan; this is not external timestamped preregistration.

| Unique requests in the final snapshot | Count |
|---|---:|
| Dataset lookups | 13 (12 nonempty, 1 successful empty) |
| Perturbation tests | 222 |
| Finite ES/padj and nonzero ES | 165 |
| Noncomputable test responses | 57 |
| Failed requests | **0** |
| Total unique requests | **235** |

Zero request failures does not mean complete biological evidence. Noncomputable
responses are preserved. Shared queries can generate multiple pathway ledger
records; none of these counts are independent biological replication counts.

## Paired implementation comparison: hold the cap fixed

| Cap | Common pathways | Grade changes | Verdict changes | Reason changes |
|---|---:|---:|---:|---:|
| 8 / 4 | 120 | **0** | **0** | 8 |
| 5 / 2 | 75 | **0** | **0** | 5 |

Selected driver/dataset/ES/padj are also unchanged at each fixed cap. The reason
changes concern STAT1: historical `no-dataset` becomes `member-derived`, with
8 (wide) or 5 (narrow) real nonsignificant test outcomes recorded by the fixed
implementation. The repair fixes a false coverage diagnosis; it does not turn
STAT1 into a supported candidate.

Historical native full ledgers and lookup attributes do not exist, so their counts
are **NA, not zero**. The shared endpoint below uses only each returned selected
row's ES, adjusted p < .05 and known pathway direction. It does not reconstruct a
historical full ledger. All finite/failed raw-query coverage is reported separately.

| Implementation / cap | A / B / C pathways | Profiles with selected-row support | Profiles with selected-row conflict |
|---|---:|---:|---:|
| Historical 8 / 4 | 19 / 39 / 62 | 8 / 17 | 3 / 17 |
| Fixed 8 / 4 | 19 / 39 / 62 | 8 / 17 | 3 / 17 |
| Historical 5 / 2 | 12 / 29 / 34 | 8 / 17 | 3 / 17 |
| Fixed 5 / 2 | 12 / 29 / 34 | 8 / 17 | 3 / 17 |

Each condition retains 17 profiles, 15 graded and 2 without enrichment, and has no
implementation error. The unchanged grades on this reused corpus do not negate
regression-tested edge-case repairs; they also do not establish improved accuracy.

## Paired cap comparison: hold the implementation fixed

For either implementation, reducing 8/4 to 5/2 removes **45 pathway rows**. Among
75 common pathway rows, **3 change A -> B** (3 verdict changes, 3 reason changes):

- Hypoxia: `cellular response to hypoxia`;
- Hypoxia: `cellular response to decreased oxygen levels`;
- DHT/AR: `alcohol metabolic process`.

Their selected best rows stay the same, but the narrow selection no longer
includes enough additional supporting contexts to establish two-cell-line
replication. Separately, three other hypoxia pathways (`response to hypoxia`,
`response to decreased oxygen levels`, `response to oxygen levels`) switch their
selected dataset from D23437 to D22946 and change ES/padj while retaining grade A.
Selected-row switches and grade changes are therefore distinct endpoints.
This contrast jointly changes pathway and dataset caps;
it does not estimate each cap's separate contribution. Profile-level any A/B and
selected-row support stay at 11/17 and 8/17, respectively, for both caps.

The fixed implementation's **full native ledger** reveals information the best
row alone misses:

| Fixed-code endpoint | 8 / 4 | 5 / 2 |
|---|---:|---:|
| Retained pathway/dataset outcomes | 328 | 135 |
| Support | 75 | 41 |
| Conflict | 15 | 12 |
| Nonsignificant | 122 | 41 |
| Noncomputable | 116 | 41 |
| Failed | 0 | 0 |
| Profiles with any support | **9 / 17** | **8 / 17** |
| Profiles with any conflict | 3 / 17 | 3 / 17 |
| Profiles with any finite ES/padj | 12 / 17 | 12 / 17 |

SMAD3 accounts for the extra wide-cap supported profile: **2 supports and 7
conflicts** at 8/4, versus **0 supports and 5 conflicts** at 5/2. It is mixed
profile-level evidence, not a newly confirmed causal driver. No profile had both
outcomes inferred solely from a grade. The wide-cap ledger contains more overlapping
pathways and datasets; higher outcome totals are not a measured gain in accuracy.

## Reproduction and integrity

From the package root, before running R: `export R_LIBS_USER=~/Rlib`.

```sh
# No network: requires the frozen snapshot and its recorded R/package environment.
Rscript dev_docs/mechanism-paired-4215-self-test.R
Rscript dev_docs/mechanism-paired-4215-summary-self-test.R \
  ../adjudicate-interpret/results/mechanism-paired-4215-paced
# Write another analysis without overwriting the stored one:
Rscript dev_docs/mechanism-paired-4215-summary.R \
  ../adjudicate-interpret/results/mechanism-paired-4215-paced /tmp/new-paired-analysis
```

The final archive includes its own `replay-runner.R`, complete inputs, source
snapshots, raw responses, manifests and saved replay. To independently regenerate
replay, copy the archive's inputs/sources/raw/manifests into a fresh directory
without its existing `replay` or `analysis` directories, then run:

```sh
Rscript /path/to/archive/replay-runner.R replay /absolute/path/to/fresh-snapshot
```

The runner verifies exact R and package versions, so those dependencies must be
available; this is an offline data/source archive, not a container image. The
summary runner is included separately in the distribution. Its manifest is
relative to the analysis output directory. One harmless frozen-runner message
counts 236 finalized files as “requests” because it includes the events table;
the true count is 235 request records plus one events table. The checked-in
runner corrects that message; the captured runner remains byte-for-byte frozen.
Failure, noncomputability and missing
native historical evidence remain visible in every replay.

Sources: [condition summary](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/analysis/condition-summary.tsv>),
[profile endpoints](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/analysis/profile-endpoints.tsv>),
[paired pathway counts](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/analysis/paired-pathway-summary.tsv>),
[selected dataset/statistic transitions](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/analysis/paired-selected-transitions.tsv>),
[full selected rows](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/replay/rows.tsv>),
[fixed native ledger](<../../adjudicate-interpret/results/mechanism-paired-4215-paced/replay/ledger.tsv>).

## Limits and next decision

No score formula, target-count weighting, grade rule or threshold was tuned. No
original DE recomputation, external independent validation, full historical
candidate-set rescan, or live model-compliance assessment was performed. The 17
labels and profile dependence limit generalization. This snapshot supports
implementation diagnostics and cap sensitivity—not unique causal attribution,
precision estimates, or a claim that fixes increased discovery accuracy.

The immediate practical recommendation is to inspect state, availability and the
full ledger before grades. If changing the grading method is next, first decide
how explicit support-plus-conflict should be communicated; the frozen responses
can then test that proposed rule without fresh service variation or selecting a
rule for favorable counts.
