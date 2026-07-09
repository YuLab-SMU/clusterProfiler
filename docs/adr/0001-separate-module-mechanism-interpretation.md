# 0001 - Separate Module Summaries, Mechanism Inference, and LLM Interpretation

## Status

Accepted

## Context

clusterProfiler currently provides enrichment analysis and an `interpret()` workflow that sends top enriched terms, optional context, PPI, and fold-change information to an LLM for structured biological interpretation.

The long-term product direction is to make clusterProfiler a mechanism interpretation tool rather than only an enrichment analysis tool. That requires more than LLM narration: users need traceable summaries, multiple evidence sources, auditable confidence, and eventually testable mechanism hypotheses.

## Decision

We will expose two structured interpretation layers before narrative reporting:

1. `Module Summary`: groups enriched terms into coherent biological themes.
2. `Mechanism Inference`: turns module summaries and optional evidence adapters into testable mechanism hypotheses.

`interpret()` remains the user-facing narrative and report interface. It may accept enrichment results directly for compatibility, but the documented main workflow should prefer explicit structured inputs:

```r
modules <- summarize_terms(x)
mechanisms <- infer_mechanisms(modules, evidence = ...)
interpret(mechanisms)
```

LLMs may propose mechanism hypotheses, but those hypotheses must be represented as structured mechanism objects with evidence, confidence components, and rationale. They must not exist only inside free-text narrative output.

## Consequences

- Existing `interpret(enrichResult)` users can remain supported through an automatic compatibility path.
- New engineering work should prioritize structured objects and contracts over prompt-only improvements.
- Confidence scoring can include AI judgment, but it must expose component scores and rationale.
- Evidence adapters can extend interpretation quality without becoming mandatory for the offline default path.
- Tests must cover both structural contracts and curated golden cases for biological usefulness.
