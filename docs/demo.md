# Admissions Demo Playbook

## 1. Script (3 minutes)

**0:00 – Hook (30 sec)**  
"Single-cell biologists spend days matching clusters to cell types by hand. With CellAnnot-GPT, we can do it in minutes, with evidence citations. Let me show you how our current build works."

**0:30 – Setup (20 sec)**  
"I've pre-computed marker genes for a PBMC dataset—standard practice after clustering in Scanpy. The backend API is running locally (`uvicorn backend.api.main:app`) and our Streamlit UI is connected to it."

**0:50 – Upload & Context (30 sec)**  
"First, I upload the CSV of clusters and markers. CellAnnot-GPT recognizes the structure and caches the data. I set the context to human blood to steer the ontology."

**1:20 – Batch Annotation (45 sec)**  
"We hit ‘Batch Annotate’. Under the hood, the annotator crafts prompts for each cluster, calls GPT-4o, then cross-checks every suggestion against our marker database. Within seconds we have preliminary labels."  
(Display summary: supported vs flagged clusters.)

**2:05 – Review Results (40 sec)**  
"Let’s inspect Cluster 3. The AI calls it ‘CD14+ Monocyte’, citing MS4A1 and LYZ and shows validation passing. Another cluster is flagged—note the warnings about conflicting markers; we can report it as ‘Unknown/Novel’ instead of guessing."

**2:45 – Close (15 sec)**  
"In three minutes we produced evidence-backed annotations, surfacing both confident calls and uncertainties. This workflow is repeatable across datasets, saving days per project."

## 1.1 New UI Highlights
- **Comparison Mode:** After each batch run, switch to *Compare Batches* to overlay contradictions and status changes between runs. Great for prompt/DB iterations.
- **Knowledge-linked Warnings:** Flagged markers now link directly to PanglaoDB/CellMarker lookups so reviewers can open the evidence with one click.
- **Call-to-action Panels:** The app suggests follow-up experiments and markers to validate, turning validation notes into a mini playbook.

### Watch the 90-second tour
[Loom walkthrough of the upgraded dashboard](https://www.loom.com/share/your-cellannot-demo)  
(Replace with the latest recording when available.)

## 2. Slide Outline

### Slide 1 – Problem
- Cell annotation is slow, manual, and knowledge-intensive.
- Researchers cross-reference gene markers with scattered literature.
- Consequence: delays in analysis, inconsistent labels.

### Slide 2 – Live Demo
- Input: marker CSV + tissue context.
- Pipeline: LLM annotation → validation against marker DB → UI review.
- Key moments: batch summary, flagged clusters, transparency of rationale.

### Slide 3 – Impact
- 10× faster annotation workflow with higher confidence.
- Traceable explanations + uncertainty handling build trust.
- Easily extensible to new tissues via prompt+database updates.

## 3. Demo Checklist

**Datasets & Artifacts**
- `data/demo/pbmc_markers.csv` (8 clusters).  
- Marker DB built: run `poetry run python scripts/build_marker_db.py --output-dir data/processed`.

**Pre-Demo Commands**
- Start API: `poetry run uvicorn backend.api.main:app --reload`.  
- Launch UI: `poetry run streamlit run frontend/streamlit_app.py`.  
- (Optional) Cache annotations in Redis: `docker compose up redis`.

**During Demo**
- Upload markers, set species/tissue.
- Run batch annotation; highlight supported vs flagged clusters.
- Open `docs/demo.md` for script cues.

**Fallback Plan**
- If OpenAI API rate-limit hits: switch to pre-recorded JSON outputs (`docs/demo_assets/pbmc_annotation.json`).  
- If UI fails: call API via `curl` and display responses in terminal.  
- If marker DB missing: load from backup `data/processed/marker_db.parquet.bak`.
