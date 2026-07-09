# Domain Context

## Glossary

### Mechanism Interpretation Tool
A workflow that turns enrichment analysis outputs into progressively stronger biological explanations, starting with traceable mechanism summaries and evolving toward testable mechanism hypotheses.

### Mechanism Summary
A structured, evidence-linked synthesis of enriched terms into higher-level biological themes. It explains what the enrichment results collectively suggest, without claiming causal proof.

### Mechanism Hypothesis
A testable biological claim that proposes a directional relationship among regulators, processes, genes, phenotypes, or experimental context. It must expose its supporting evidence and uncertainty.

### Evidence Adapter
A pluggable source of supporting biological evidence, such as pathway hierarchy, gene annotation, PPI, literature, perturbation databases, or model-derived signals. Evidence adapters may enrich the interpretation workflow but are not required by the offline default path.

### Statistical Evidence Adapter
The default evidence adapter that derives support from enrichment result fields, such as adjusted p-values, NES, gene ratios, term overlap, shared genes, and direction consistency.

### Knowledge Hierarchy Adapter
An optional evidence adapter that uses ontology or pathway structure, such as GO hierarchy, KEGG pathway categories, or WikiPathways relationships, to support module grouping and biological interpretation.

### PPI Network Adapter
An optional evidence adapter that uses protein-protein or functional interaction networks to support hubs, modules, or regulator-process relationships.

### PubMed Literature Adapter
An explicit opt-in evidence adapter that retrieves literature support from PubMed and records cached evidence, including query, retrieval date, PMID, title, and supporting abstract snippets.

### Offline Default Path
The stable, reproducible mechanism interpretation path that uses only the enrichment result and user-provided analysis context. It must remain usable without network access or external services beyond the optional LLM interpretation step.

### Module Summary
A public intermediate object that groups enriched terms into coherent biological themes. It is produced by term summarization and is consumed by interpretation, reporting, and mechanism inference.

### Mechanism Inference
A public workflow that turns module summaries and optional adapter evidence into mechanism hypotheses. It is separate from term summarization because it may involve biological interpretation, directionality, uncertainty, and testability.

### Comprehensive Confidence Score
A mechanism-level confidence assessment that may combine statistical support, adapter evidence, literature or knowledge signals, and AI judgment. It must remain auditable by exposing the evidence and rationale that contributed to the score.

### Golden Case
A curated real or realistic enrichment-analysis example with an expected biological interpretation. Golden cases are used to evaluate whether mechanism summaries and mechanism hypotheses are biologically useful, not just structurally valid.
