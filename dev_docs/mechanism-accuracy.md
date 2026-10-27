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

A fair question: GSEA uses the whole ranked list, so it should carry more than
a threshold-filtered gene list. Measured, it does not.

`evidence()` now accepts a `gseaResult` and takes:

- **genes** from `core_enrichment` — the genes that produced each enriched
  set's ranking, union across the sets the caller kept. Not `@geneList`, which
  is every gene including the ten thousand that barely moved: asking which
  regulator explains all 22899 returned **0 supported** regulators on the
  airway list, because no factor stands out against the whole transcriptome.
  On that list `core_enrichment` covers 4594 genes, against the 897 the ORA
  path uses.
- **universe** from the ranked list — the background GSEA itself tested
  against. (`gseaResult` has no `@universe` slot, unlike `enrichResult`.)
- **keytype** from `@keytype`, mapped as for `enrichResult`.

Same 25 ground-truth datasets, both paths:

| Metric | ORA | GSEA |
| --- | --- | --- |
| True regulator appears | 44% | **44%** |
| recall@10 | 4% | **4%** |
| recall@20 | 8% | **8%** |
| MRR | 0.047 | **0.046** |
| Median rank | 157 | **155** |

Per case they are near-identical (FOXO4 162 vs 161, ADNP 168 vs 163, ATF6 157
vs 155). Both land at the same rank for the cases that fail.

That is the useful finding: **the limitation is not the threshold.** Feeding
GSEA's continuous ranking or ORA's binary selection into the same
over-representation test gives the same answer, because the test is what has no
power at an overlap of six genes. Changing the front end cannot fix the back
end.

What differs between the two inputs on the airway list is the candidate pool:
`core_enrichment` makes 384 regulators significant against 55 for the 897-gene
ORA list, so a fixed candidate cap of 30 admits a different, smaller fraction of
the truthful ones. `GR` ranks 52nd under GSEA and 57th under ORA — outside the
cap either way, which is why it needed the cap raised rather than a better front
end.
