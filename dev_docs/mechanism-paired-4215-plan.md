# Frozen paired perturbation comparison: prespecified design

## Operational amendment after the first collection

The first immutable pilot snapshot (`mechanism-paired-4215`) retained 235 unique
requests, including 85 HTTP 429 responses. It is preserved, not repaired or used
as a biological negative reference. A separate `mechanism-paired-4215-paced`
snapshot will repeat the same four conditions and all 17 inputs, with at least
two seconds between HTTP requests (including MCP session initialization). Every
request in that second snapshot is collected afresh; no successful pilot result
is selectively substituted for a failed request. Failures in each snapshot stay
frozen. Do not compare the snapshots to estimate a code effect: only the four
conditions within each identical frozen cache support that contrast.

The revised runner additionally treats unreadable/malformed/mismatched cached
records as fatal infrastructure errors, not service failures that either grader
may swallow. The entire plan is included in the revised frozen manifest. The
pilot's original frozen manifest did not bind this full plan, so the archive is
not evidence of external or cryptographically timestamped preregistration.
Endpoint definitions, profiles, pinned code and caps are unchanged by this
operational amendment.

Design recorded before response collection on 2026-11-11. This is a deterministic
implementation and coverage experiment on a reused exploratory corpus, not an
independent accuracy or causal-identification benchmark.

## Inputs and scope

- Retain all 17 original profiles and their prespecified expected labels from the
  frozen benchmark metadata. The expected label is condition-informed and is not
  an exhaustive causal gold standard.
- Grade the expected driver of each profile. Do not choose candidates using the
  observed perturbation results. This does not measure shortlist precision.
- Retain profiles with no enrichment or unusable measurements in the denominator
  and explicitly record why no result can be interpreted.
- Reuse the original signed vectors and enrichment objects; do not recompute DE
  or enrichment during this comparison.

## Four paired conditions

| Code | Pathway cap | Datasets per driver |
|---|---:|---:|
| Historical baseline `978eabaee502350e12e5b6a3988c54091ac79fe3` | 8 | 4 |
| Historical baseline | 5 | 2 |
| Fixed implementation `fe7c060` | 8 | 4 |
| Fixed implementation | 5 | 2 |

The 8/4 condition matches the historical verification script, not the package's
current defaults. Both implementations receive the same frozen raw lookup and
test responses for each identical request. Code differences may change which
requests are selected; that is part of the code-policy contrast, not a claim that
the selected datasets are identical. The snapshot records the union of requests.

Compare code at a fixed cap, and caps within a fixed implementation. Do not
attribute the historical-versus-bounded difference to one factor without showing
these separate contrasts. An interaction between code and caps is possible.

## Freeze and replay contract

Store original source snapshots, complete input objects, exact query keys and
responses (including successful empty responses, failures and non-computable
statistics), collection timestamps, provenance and file checksums. Identical
requests reuse the first retained response; never replace an observed response
merely because another response looks better. Collection may resume unfinished
requests, but finalized snapshots are immutable.

Replay is offline. A missing request is an incomplete snapshot and must fail
rather than query the service or fabricate an empty response. Save all four
condition results, query traces and paired transition tables. A fresh replay must
produce identical deterministic result files; exclude collection timestamps from
those files. Distinguish successful collection/replay from biologically usable
evidence.

## Endpoints and reporting

1. **Coverage:** all 17 profiles; presence of enrichment; candidate lookup status;
   finite test statistics; failed and non-computable outcomes.
2. **Directional evidence:** support, conflict, both, non-significance and unknown,
   using the exact frozen query response and pathway direction. Old code without
   a ledger must not be credited with a complete native evidence disclosure.
3. **Implementation outputs:** pathway identity, selected driver/dataset,
   A/B/C grade, original verdict/reason and paired transitions. These are not
   accuracy or passing-rate endpoints.
4. **Sensitivity:** profile-level any support/conflict and evidence completeness,
   separately from any A/B. Overlapping pathways and shared datasets are not
   independent biological confirmations.

Do not alter score formulas, A/B/C rules, thresholds or dataset-selection policies
to improve these results. Evidence-first presentation is a separate compatible
change; it should expose the distinction rather than regrade the observations.
