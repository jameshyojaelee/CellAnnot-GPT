# Why GPT Cell Annotator

Single-cell RNA-seq studies routinely generate dozens of clusters that must be mapped to known cell types. The manual process—cross-referencing marker tables, literature, and ontologies—can take days and often produces inconsistent labels. GPT Cell Annotator combines curated knowledge bases, retrieval-augmented prompting, and validation guardrails to deliver faster, auditable annotations.

## State of the Field (2024)

| Approach | Strengths | Gaps | Representative Tools |
| --- | --- | --- | --- |
| **Reference-based matching** | Fast, reproducible, deterministic; works well when query data matches reference atlases. | Struggles with rare or novel cell states; heavily dependent on curated references; limited explanations. | Seurat label transfer, SingleR, CellTypist |
| **Latent-transfer / deep embedding** | Captures complex expression patterns; can learn cross-dataset invariances. | Requires expression matrices and GPU resources; latent spaces can be opaque; confidence estimates often ad hoc. | scANVI, scArches, scVI-tools |
| **LLM + validation (GPT Cell Annotator)** | Flexible prompts digest marker summaries, ortholog mapping, and literature snippets; explicit rationales; validation downgrades hallucinated calls. | Dependent on high-quality marker rankings; live LLM access introduces cost/latency without offline fallback; requires guardrail tuning. | GPT Cell Annotator, emerging academic prototypes |

Traditional pipelines either rely on close matches to existing references or require training sophisticated latent models. GPT Cell Annotator complements these approaches by:

- Translating marker lists (including orthologs and synonyms) into ontological suggestions.
- Surfacing supporting evidence and validation warnings for each cluster.
- Falling back to heuristic annotations when API access is unavailable (mock mode).

See [`docs/benchmarks.md`](docs/benchmarks.md) for empirical comparisons against marker-overlap baselines and [`docs/operations.md`](docs/operations.md#validation-guardrails) for details on the guardrail design.

## Differentiators

- **Knowledge-aware prompts**: Retrieval pulls evidence from the marker database before each LLM call, keeping responses anchored in curated markers.
- **Structured validation**: Results are checked against configurable overlap thresholds, ontology expectations, and species mappings, ensuring “Unknown or Novel” is surfaced when evidence is weak.
- **Offline-first packaging**: Bundled assets and CLI commands (`gca annotate`, `gca build-db`) let users run demos or smoke tests without external dependencies.
- **Explainable outputs**: Every annotation includes rationale sentences, supporting markers, ontology IDs, and mapping notes, which can be exported alongside Scanpy pipelines.

## Roadmap Snapshot

Highlights from the near-term roadmap (see [`docs/roadmap.md`](docs/roadmap.md) for the full list):

- Expand the marker knowledge base with tissue-specific sources and literature scraping.
- Introduce active-learning loops that propose validation experiments or additional markers.
- Support lightweight local LLMs for organisations that cannot send data to external APIs.

With these pillars, GPT Cell Annotator delivers a balanced workflow: fast automated suggestions, transparent evidence, and configurable safeguards that keep biologists in control.
