# Frozen mechanism-discovery audit (4.21.5 work)

## Scope and reproducibility

This is an **exploratory, reused development corpus, not independent validation**.
It audits ranking denominators, repeated studies, module-set size and the meaning
of historical verification. It neither changes the discovery score nor runs live
GPSAdb/network requests. Baseline repository commit:
`978eabaee502350e12e5b6a3988c54091ac79fe3`.

Run from the `clusterProfiler` package directory:

```sh
export R_LIBS_USER=~/Rlib
Rscript dev_docs/mechanism-benchmark-4215.R --self-test
# Optional positional root: a complete copy of the frozen audit directory.
Rscript dev_docs/mechanism-benchmark-4215.R /path/to/audit-benchmark-4215 --self-test
```

Default artifacts are in the sibling
`../adjudicate-interpret/results/audit-benchmark-4215/`:

- `inputs/profiles.rds`: only the 17 signed vectors, original labels and declared
  drivers, without enrichment tables or redundant gene universes (2,889,700 bytes).
- `inputs/metadata.tsv`: expected labels, declared drivers, activation/KO design,
  manually evidenced study groups and generation-script line references.
- `inputs/collectri_net.tsv`, `inputs/p6_discover_verify.tsv` and
  `inputs/p7_module_all17.tsv`: byte-identical frozen source tables.
- `inputs/original-manifest.tsv`: MD5/size of the five original source inputs;
  `input-manifest.tsv`: MD5/size of every frozen input. The latter is checked on
  every run. Inputs are never rewritten during normal execution.
- `outputs/ranks.tsv`: full-universe expected-label ranks and explicit status;
  `full_rankings.tsv`: every scored candidate; `metrics.tsv`: counts/rates for
  all17, design, study and driver; `grouped_macro.tsv`: equal-group averages.
- `outputs/modules.tsv`, `module_summary.tsv`: candidate sets, module size,
  label coverage and exact size-matched random expectation.
- `outputs/historical_verification.tsv`, `historical_summary.tsv`: frozen old
  candidate calls, not rerun current grading.
- `outputs/implementation.tsv`: actual sourced implementation MD5 and baseline
  commit; `output-manifest.tsv`: checksums for generated tables.

There are no `/tmp` dependencies. Base R alone is sufficient: the script sources
`R/mechanism_discover.R` directly, avoiding the package's load-time dependencies.
Makefile was read; it has no benchmark target. No package installation/build,
`make rd`, or `make check` is required for this source-only audit. Do not use this
script as evidence that the package-wide test/check workflow passed.

The cache was initially created with `--freeze --self-test`; `--freeze` reads
`all11_de.rds`, `p3/p3_new_de.rds` and the three TSVs from the parent results
directory and **refuses to overwrite an existing inputs directory**. For routine
reproduction, use the already frozen inputs and omit `--freeze`. The manifest
makes accidental changes detectable; MD5 is an integrity check, not an authenticity
signature. Source-generation scripts are provenance references, not execution
requirements for rerunning the cache.

## Labels, study dependence and denominator

The expected label is the first declared driver, with `GR` manually normalized
to `NR3C1`. Hypoxia originally declares `HIF1A,EPAS1,ARNT,VHL`; only `HIF1A` is
the prespecified exact-match label here. These condition-informed labels are
**not exhaustive causal gold standards**: a nonmatching TF can be a real parallel
or downstream response factor. Thus “expected-label match” is used rather than
claiming causal precision/PPV.

Study groups come from the actual local generation scripts, not a guess that each
contrast is independent. `36_all11_de.R:110–119` identifies the ESR1 KO as
GSE246363. `50_p3_new_datasets.R:74–91` places both Oci-AML3 timepoints and HEK293T
in GSE245918; lines 93–107 place IFN-gamma and LPS in GSE348821 (shared controls);
lines 109–122 identify A549 TGF-beta1 as GSE351288. Other explicit GSE accessions
are retained from their generation blocks. Airway has no accession in the saved
metadata/generation block and is **unknown**, not assigned a guessed unique study.
The metadata table records each mapping's evidence.

There are 17 profiles, 13 known study groups plus one unknown-study profile, and
13 expected-label groups. GSE245918 contributes three profiles; GSE348821 two.
Even different GSE IDs do not prove biological independence or independence from
CollecTRI's literature. Known-study macro rates exclude the unknown profile and
disclose that exclusion; all17 micro rates always include it. Grouped rates are
descriptive, not held-out/group-cross-validated performance.

Every call requests at least the number of network TFs, then asserts that returned
rows equal the scorer's actual `n_factors_tested`. The label status is one of
`scored`, `missing_from_network`, or `unscorable`; `truncated` is explicitly false.
Missing/unscorable labels count as misses, never disappear from the denominator.
**All 17 labels in these frozen inputs are actually scorable.** Earlier top20
absence must not be recast as network absence or an unavailable truth.

## Fixed ranking design and results

Settings are unchanged: all measured signed fold changes, `min_targets=5`, removal
of self-targets, no size weighting, original signed Pearson/t scoring. Two
prespecified rankings are shown: descending signed score and descending magnitude
(`abs`). A design-aware rule uses signed score for the 16 activation/condition
profiles and magnitude for the one ESR1 KO. No rank-dependent choice or threshold
tuning is made in this audit.

| Rule / stratum | Profiles | hit@1 | hit@5 | hit@10 |
|---|---:|---:|---:|---:|
| Signed, full corpus | 17 | 6/17 (35.3%) | 10/17 (58.8%) | 13/17 (76.5%) |
| Signed, activation | 16 | 6/16 | 10/16 | 13/16 |
| Signed, KO | 1 | 0/1 | 0/1 | 0/1 |
| Magnitude, full corpus | 17 | 4/17 | 8/17 | 10/17 |
| Magnitude, activation | 16 | 4/16 | 8/16 | 10/16 |
| Magnitude, KO | 1 | 0/1 | 0/1 | 0/1 |
| Prespecified design-aware, full corpus | 17 | 6/17 | 10/17 | 13/17 |

One KO profile cannot validate an activation-versus-KO recommendation. Magnitude
ranking moves ESR1 from signed rank 639 to magnitude rank 41 (653 scored TFs),
but does not recover its expected label in the top10 here. The signed
known-study macro rates are **34.6%, 55.1%, 69.2%** at 1/5/10; driver-macro rates
are **38.5%, 59.0%, 73.1%**. These averages prevent repeated XBP1 contrasts from
being treated as independent replications but do not remove corpus selection bias.

## Modules: report set coverage together with set size

The historical module definition is preserved, not optimized: the signed top1
plus its direct Jaccard-neighbours at **0.10**, restricted to that dataset's
historical **top-five candidate set**. This is not transitive clustering and not
a search over the whole network. The script reconstructs each set using the frozen
network, verifies equality to the frozen modules, and verifies that current
signed top5 equals the historical candidate order. The threshold was explored
on this same corpus historically; fixing it now does not make it independent.

Across all17, singleton top1 matches **6/17**; module label coverage is **7/17
(41.2%)**, with mean set size **1.882**. Labels occur within the restricted top5
in 10/17 profiles. Restricting analysis to those ten would make the same coverage
look like 7/10; it is not the headline denominator.

The size-matched random baseline **always retains the same top1**, then adds a
uniform `(module_size - 1)`-subset from the **same remaining four candidates**.
Its exact expected label coverage is 1 when top1 matches, otherwise
`I(label in remaining four) * (module_size - 1) / 4`. This is label-aware
**evaluation of an expectation**, not a random method allowed to select the
reference label. It is not a Monte Carlo realization. Self-tests exhaustively
enumerate all top1-anchored size-matched subsets to check the formula. This
baseline conditions on discovery's shortlist and singleton anchor, not random TF
selection from the network. Its expected coverage is **6.5/17 (38.2%)**, compared
with observed module coverage **7/17 (41.2%)**: only 0.5 expected matches more.
Coverage rises partly because sets are larger. This is **set coverage, not PPV
gain, causal validation, calibrated confidence, or an independent module test**.
The old “14/17 correct or honestly refused” composite is not a benchmark endpoint.

## Historical verification is separate from current fixed grading

The frozen table contains 85 candidate rows (17 × 5). Its old A/B calls pass
37 rows; **7/37** are exact expected-label matches. This is a historical
label-agreement description, not causal PPV, and nonmatches are not established
false mechanisms. It is **not** a rerun of today's corrected direction,
provenance, study-group, reproduction or preflight logic. The benchmark performs
no new verification and cannot establish how many of those historical calls would
pass current grading. Frozen DE/global-effect fields are retained for provenance
but are not used for ranking, exclusion, refusal credit or endpoint tuning.

## Checks and limits

`--self-test` asserts full17 denominators, all labels scorable, complete untruncated
rankings, signed hit counts 6/10/13, module sets and 7 covered labels, known shared
study groups, the single KO, historical 7/37 agreement, and exhaustive random-set
expectations. Regenerating outputs twice must leave `output-manifest.tsv`
byte-identical (including the sourced implementation MD5).

Not checked here: original DE recomputation, live source inventory, external
validation, causal adjudication, or package-wide `make check`. The corpus and
historically selected module threshold remain exploratory; no new performance or
independence claim is justified by this audit.
