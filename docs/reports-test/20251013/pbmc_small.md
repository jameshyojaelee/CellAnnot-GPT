# CellAnnot-GPT Benchmark Report

- Accuracy: 50.00%
- Macro F1: 33.33%

## Per-class Metrics

- **B cell**: Precision 100.00%, Recall 100.00%, F1 100.00%
- **T cell**: Precision 0.00%, Recall 0.00%, F1 0.00%

## Prediction Details

- Cluster 0: predicted B cell (truth B cell)
- Cluster 1: predicted CD8 T cell (truth T cell)

Confusion Matrix (rows=truth, cols=prediction)
 	B cell	CD8 T cell	T cell
B cell	1	0	0
CD8 T cell	0	0	0
T cell	0	1	0