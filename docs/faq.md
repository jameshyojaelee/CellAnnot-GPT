# Frequently Asked Questions

## How does GPT Cell Annotator handle gene synonyms and species differences?

The annotator normalises marker lists before retrieval and prompting. Synonyms and orthologs are
defined in `config/gene_synonyms.json`. By default, aliases (e.g., `CD20` → `MS4A1`) and select
mouse→human ortholog mappings are resolved, so prompts and validation operate on canonical IDs.
Set `SYNONYM_CONFIG_PATH` to point to a custom JSON file or disable ortholog expansion with
`SYNONYM_ENABLE_ORTHOLOGS=false`.

## Can I disable retrieval-augmented prompting?

Yes. Set `RAG_ENABLED=false` to fall back to the legacy prompts. You can also adjust
`RAG_TOP_K` and `RAG_MIN_OVERLAP` to limit the number of candidates surfaced from the marker
database.
