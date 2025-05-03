# CNA Inferencer for scRNA-seq Data

This Python package provides a method to infer **Copy Number Alterations (CNAs)** from single-cell RNA-sequencing (scRNA-seq) data. CNAs are identified as genomic gains or losses and assigned to individual cells or cell groups using a transparent sliding-window approach. The tool is designed to work with `AnnData` objects and outputs structured CNA annotations for visualization and downstream analysis.

---

## Project Tasks

This tool addresses the following project objectives:

### Task 1: CNA Inference

* **Input**: an `AnnData` object containing:

  * Raw expression counts in `adata.layers["counts"]`.
  * Gene genomic coordinates in `adata.var["chromosome"]`, `adata.var["start"]`, `adata.var["end"]`.
  * Genotype labels in `adata.obs["simulated_cnvs"]`.
* **Output**: CNA calls as a pandas DataFrame (and stored in `adata.uns["cna_calls"]`) with columns:

  * `chrom`: chromosome of the event
  * `start`, `end`: genomic coordinates (basepairs)
  * `type`: "gain" or "loss"
  * `cells`: list of cell barcodes harboring the event
  * `obs_label`: genotype group for which the call was made
* **API**:

  ```python
  from mycnatool.infer import infer_all_cnas
  import scanpy as sc

  adata = sc.read_h5ad("PBMC_simulated_cnas_041025.h5ad")
  cna_calls = infer_all_cnas(
      adata,
      window_size=5,
      threshold=0.2
  )
  adata.uns['cna_calls'] = cna_calls.to_dict('records')
  ```
* **Workflow Overview**:

  1. **Reference detection**: identifies the most frequent genotype as the diploid baseline
  2. **Log₂‐fold-change**: computes per‐gene fold‐changes vs. that baseline (with +1 pseudocount)
  3. **Genome ordering**: builds and sorts a table of gene positions + mean log₂FC
  4. **Sliding-window smoothing**: applies a moving average across adjacent genes to reduce noise
  5. **Threshold & segmentation**: finds contiguous windows above/below threshold, merges them into CNA regions
  6. **Cell assignment**: records which cells cross threshold within each region
* **Visualization**: extract the inferred regions from `adata.uns['cna_calls']`, subset the AnnData to overlapping genes, and use the built-in heatmap utilities to plot log₂FC grouped by genotype.
* **Novelty**: a pure-Python, transparent method requiring only `scanpy`, `numpy`, and `pandas`; no external HMM or PCA dependencies; and automatic reference selection.

### Task 2: Performance Assessment

* Evaluate method accuracy on test datasets using precision, recall, and AUC metrics.
* Explore the effect of read depth via downsampling.

### Task 3: Application to PSC scRNA-seq

* Apply CNA inference to selected pluripotent stem cell (PSC) datasets to detect known and novel CNAs.

### Optional (Extra Credit)

* Simulate gold-standard data for benchmarking (Task 2B).
* Predict downstream functional impacts of CNAs on PSCs.

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/RuihanYuan/scRNA_CNA_detection.git
cd scRNA_CNA_detection
pip install -e .
```

Dependencies:

* Python 3.8+
* `scanpy`
* `numpy`
* `pandas`

---

## Usage Example

```python
import scanpy as sc
from mycnatool.infer import infer_all_cnas

adata = sc.read_h5ad("PBMC_simulated_cnas_041025.h5ad")
cna_calls = infer_all_cnas(adata, window_size=5, threshold=0.2)
adata.uns['cna_calls'] = cna_calls.to_dict('records')
print(cna_calls)
```

---

## License

MIT License
