# How accurate is the upstream half?

An assessment of whether `evidence()` finds the regulator that actually caused a
gene list, measured against known answers rather than asserted.

## What was measured

GPSAdb records which gene each dataset perturbed. That is a ground truth: if a
dataset knocked down `X` and produced a set of differentially expressed genes,
then `X` is the regulator of that list **by construction**.

So the test is:

1. Take a GPSAdb dataset whose perturbed gene is a real gene in C3:TFT.
2. Feed its significant DE genes to `evidence(upstream = "tf")`.
3. Ask whether the perturbed gene comes back.

Thirty datasets, drawn at random from the 54 that qualify. The ground truth sits
in GPSAdb; the recovery runs against C3:TFT, which is built from GTRD and
TRANSFAC — a different database, so the test is not circular.

## Result

| Metric | Value |
| --- | --- |
| Test cases | 30 |
| True regulator among the significant entities | **37%** (11/30) |
| recall@1 | **3%** |
| recall@5 | 3% |
| recall@10 | 3% |
| recall@20 | 7% |
| MRR | **0.039** |
| Median rank of the true regulator | **157** |

Being right about one in three, and usually far down the list, is not a
ranking you would act on.

## What the controls say

**Specificity is not the problem.** Random gene sets of the same size, through
the same code:

| | Real perturbation signatures | Random gene sets |
| --- | --- | --- |
| Significant entities, median | **275** | **0** |
| True regulator recovered | 37% | **0%** |

Random lists produce nothing at all. The channel detects real structure; it is
not a generator of plausible-looking noise. On the airway dexamethasone list the
same contrast holds: 225 significant target sets against 0 for each of three
random lists.

**Missing data is not the problem either.** The true factor's own target set
does overlap the input — median 13 genes, and only 2 of 30 cases had zero
overlap.

## The actual cause: statistical power

The true factor's target set *is* enriched in its own knockdown signature:

- median observed/expected overlap ratio **11.1×**
- **93%** of cases above 1×

What separates the cases that came back from the ones that did not is the
absolute overlap, not the fold enrichment:

| | Overlap | Fold | p | padj | Significant |
| --- | --- | --- | --- | --- | --- |
| Recovered (n=5 with a testable set) | 28 genes | 14.9× | 0.0033 | 0.021 | **80%** |
| Not recovered (n=9) | **6 genes** | 9.1× | **0.187** | 0.318 | **0%** |

A target set holding six of the input genes is nine-fold enriched and still
fails the hypergeometric test, because the collections list hundreds of genes
each and the test has no power at that overlap. For 16 of the 30 cases the true
entity produced no testable result at all.

This is a property of over-representation analysis on large collections, not of
the grading on top of it. The signal is real and invisible at the same time.

## Bugs the assessment found

Measuring accuracy is what surfaced these; they had been invisible.

1. **Concurrent queries returned partial data silently.** Twenty parallel
   lookups of one dataset gave the complete answer 16 times and a truncated one
   4 times — same request, same status, no warning. The supported set flipped
   between runs as a result. The sweep is now serial by default.
2. **`GR` was classified as a family label.** TRANSFAC names a nuclear receptor
   by its protein; none of `GR`, `PXR`, `FXR`, `HIF2A` are HGNC symbols, so each
   was called a family — and families are never queried against a perturbation
   database. On the airway dexamethasone list this left the glucocorticoid
   receptor at `single_source` level B, while GPSAdb held two `NR3C1` datasets
   moving the list at padj 1e-22.
3. **The lookup used the collection's spelling, not the database's.** C3:TFT
   says `GR`; GPSAdb says `NR3C1`; querying `GR` returns nothing, which is
   indistinguishable from "never perturbed".

After the fixes the airway list carries 24 supported regulators instead of 22,
`GR` among them at level A — the answer a glucocorticoid experiment is
guaranteed to involve, previously missed for two unrelated name reasons.

## What this means in practice

**Sound:**

- Specificity. The channel does not fire on noise, and `absent` / `untested` /
  `failed` are kept apart so a negative is not a silent one.
- The `supported` tier means what it says: two independent sources agreeing, one
  of them measured.
- Family labels are separated from genes rather than reported as candidates.

**Not sound:**

- The output is **not a ranking to act on for "the" regulator**. The true factor
  is usually not in the head of the list.
- Recall inherits the power limit above. Raising the candidate cap raises the
  supported count (22 → 39 at cap 60 on the airway list) but there is no
  evidence that the extra names are truer.

**Suggestion:** report the fold enrichment beside the p-value. In this benchmark
a nine-fold enrichment failed significance for want of overlap, and a reader
shown only `padj = 0.19` would discard a real signal.

**A caveat on the test itself.** "Recover the knocked-down factor from its own
DE list" is a hard task, and harder than the intended use. A knockdown's
signature is dominated by downstream consequences, so the direct regulator need
not dominate it. The benchmark bounds what the method can be trusted for — it
does not show the mechanism reports on well-powered, on-target comparisons are
wrong. On the airway list the expected answer was recovered.

## GSEA instead of ORA

`evidence()` accepts a `gseaResult` and tests each collection **with GSEA**, not
with ORA. The first implementation took `core_enrichment` and ran `enricher()`
on those genes, which silently discarded the ranking and the sign — the "GSEA
path" was the ORA path with a different input, and measuring the two against
each other compared nothing. That is fixed: the ranked statistic travels from
`@geneList` into the collection stage, and each set is tested with `GSEA()`.

What the correct version adds is **direction**, which ORA cannot supply: a
binary gene list has no sign, so on the ORA path every collection row carries
`direction = NA`. On the airway list GSEA labels all 1095 collection rows, 770
up and 325 down. A regulator that a histone deacetylase inhibitor *induces* and
one it *represses* were previously the same row.

### The two methods can disagree, and both are right

Same list, same collection, opposite conclusions:

| Set | ORA padj | GSEA padj | GSEA NES |
| --- | --- | --- | --- |
| GR_Q6 | **0.0002** | 0.166 | 1.27 |
| GR_01 | 0.0031 | 0.285 | 1.27 |

The reason is measurable rather than mysterious. GR's set holds 236 genes, of
which **27 are among the 897 significant genes** — a strong over-representation,
which is exactly what ORA tests. Across the full 22899-gene ranking those same
236 genes sit at a **median rank of 10571**, barely ahead of the 11450 expected
by chance: the set is over-represented at the top *and* spread thinly through
the middle, and GSEA scores the whole distribution.

Neither is wrong. ORA asks "are this set's members over-represented among what
changed", GSEA asks "does this set move together". For a glucocorticoid
experiment the first is the question you want, which is why the ORA path finds
the receptor and the GSEA path does not.

### What this means for choosing a front end

- Use ORA when the question is "which regulator explains my significant genes".
  That is the deviation-from-threshold question, and it is the one that
  recovers GR here.
- Use GSEA when the question is "which regulator's programme moves
  monotonically with my phenotype", and when direction matters more than a
  head-of-list hit.
- Do not read a GSEA `absent` as evidence against a regulator that ORA
  supports. The two disagree on real signals, as `GR` shows.

## What the two front ends actually deliver on the airway list

Same data, same collection, both run through the full pipeline:

| | ORA | GSEA |
| --- | --- | --- |
| Wall time | 283 s | 562 s |
| Entities with a verdict | 168 | 64 |
| `supported` | **24** | **12** |
| Rows carrying a direction | 25 | 64 |
| Overlap | 4 of 32 | |

`GR` is `supported` under ORA and **absent** under GSEA.

### The direction GSEA adds is lost where it matters

The case for GSEA was direction: ORA's binary list has no sign, GSEA's NES
does. At the collection level that holds — 61 sets up, 19 down.

At the entity level it arrives intact (49 up, 17 down) **until a regulator is
corroborated by a perturbation**, at which point it is destroyed:

- the collection rows say `up`
- the perturbation rows say `mixed` — a knockdown moves genes both ways
- `consensus()` takes the union of the two and reports `mixed`

So every one of the 8 regulators GSEA found and ORA did not came back
`mixed`. The direction is real, and it survives exactly on the regulators
with the least evidence.

The cause is that two different questions share one column. A collection's
direction is *the regulator's inferred activity*; a perturbation's direction is
*which way the responding genes moved*. Merging them makes the second overwrite
the first. They should be two fields.

## A second example is needed, and the obvious candidates fail

**C3:TFT covers about 36% of human transcription factors** (591 gene entities
against ~1639). The gaps are not marginal — `NFKB1`, `RELA`, `JUN`, `FOS`,
`E2F1`, `TP53`, `CTCF`, `NOTCH1`, `ESR1` and `HIF1A` have **no target set at
all**. Three independently chosen examples (`NFKB1`, `JUN`, `E2F1` knockdowns)
returned zero findings for this reason alone.

**GPSAdb datasets cannot substitute for a full transcriptome.** `gpsadb_de()`
caps at 500 rows and reports no tested-gene total, so:

- as ORA input the background is unrecoverable. Substituting the whole
  annotation made the 500 genes 0.26% of the universe, against 3.92% for
  airway, and no collection reached significance — 0 `supported` across five
  knockdown datasets.
- as GSEA input the ranking is a truncated top-500, so 3 of 5 datasets produced
  no significant set at all (0, 13, 9, 0, 0).

A second example therefore needs **a second full expression dataset with a
known expected answer**. `airway` is currently the only one on this machine;
`pasilla` carries metadata but not the counts, and the `airway` package's `gse`
object is the same eight dexamethasone samples.
