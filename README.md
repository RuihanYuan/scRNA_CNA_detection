# cnatool

A lightweight Python package for inferring copy-number alterations (CNAs) from single-cell RNA-seq (scRNA-seq) data via a clustering-driven workflow.

---

## Table of Contents

1. [Overview](#overview)  
2. [Installation](#installation)  
3. [Data Preparation](#data-preparation)  
4. [Quick Start / Usage](#quick-start--usage)  
5. [API Reference](#api-reference)  
6. [Simulation / Gold-standard generation](#simulation--gold-standard-generation)  
7. [Benchmarking](#benchmarking)  
8. [Project Structure](#project-structure)  

---

## Overview

`cnatool` implements a three-step pipeline to detect CNAs from scRNA-seq:

1. **Clustering**: cluster reference (“diploid”) and query cells separately.  
2. **Cluster comparison**: flag clusters whose transcriptional dispersion or cell-overlap deviates beyond thresholds.  
3. **Per-cluster CNA calling**: for each flagged cluster, compute gene-level log₂-fold-change vs reference, smooth along the genome by sliding window, then call gains/losses.

By focusing only on “abnormal” clusters, the method reduces noise and computational cost while retaining sensitivity to subclonal events.

---

## Installation

```bash
# Clone the repo
git clone https://github.com/yourusername/cnatool.git
cd cnatool

# Create & activate environment
conda create -n cnatool-env python=3.10 -y
conda activate cnatool-env

# Install dependencies
pip install -r requirements.txt

# Install package in editable mode
pip install -e .
```

---

## Dependencies

The following Python packages are required:
	•	scanpy
	•	numpy
	•	pandas
	•	scipy
	•	scikit-learn
	•	leidenalg

---

## Data Preparation

Before running CNA inference, your AnnData must include genomic coordinates in .var:

```bash
import scanpy as sc
from cnatool.io import annotate_with_coords

adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()

# Annotate with Ensembl Biomart
adata = annotate_with_coords(
    adata,
    gene_id_col="gene_symbol",
    host="http://www.ensembl.org",
    dataset="hsapiens_gene_ensembl"
)

adata.write_h5ad("GSE243112_with_coords.h5ad")
```

That yields .var columns:

```bash
ensembl_gene_id | chromosome | start | end | strand
```

## Quick start:

```bash
import scanpy as sc
from cnatool.infer import infer_cna

# 1) Load your query data
adata_query = sc.read_h5ad("GSE243112_with_coords.h5ad")

# 2) (Optional) Provide your own diploid reference
# adata_ref = sc.read_h5ad("healthy_reference_with_coords.h5ad")

# 3) Run CNA inference
dclusters, summary_df = infer_cna(
    adata=adata_query,
    adata_ref=None,              # or adata_ref
    overlap_threshold=0.5,
    min_cells=10,
    dispersion_ratio_threshold=2.0,
    min_genes=5,
    z_threshold=1.0,
    verbose=True
)

# 4) Inspect results
print("Flagged clusters:", dclusters)
print(summary_df.head())
```
Results are saved in adata_query.uns['cna_summary'], a DataFrame with:

```
cluster | region | chromosome | start | end | strand | gain | loss | neutral
```
---
## API Reference

- **`infer_cna(...)`**  
  Detect abnormal clusters and call CNAs. Returns `(dclusters, summary_df)`.

- **`call_cna_by_region(...)`**  
  Binned, z-score CNA caller. Builds region bins if missing.

- **`add_genome_bins(var_df, bin_size)`**  
  Assigns genes to fixed-size genomic bins (`region`, `bin_start`, `bin_end`).

- **`easy_cluster(adata, ...)`**  
  Filter → normalize → log1p → PCA → KNN → Leiden clustering.
---
## Simulation / Gold-standard generation

Generate synthetic CNAs of arbitrary size/frequency:

```
from cnatool.simulation import simulate_gold_standard_cnas

regions      = ["1:1e7-2e7", "6:2.5e7-3.5e7", "X:0-1.55e8"]
frequencies  = [0.10, 0.15, 0.05]
copy_numbers = [0, 1, 4]

adata_gold = simulate_gold_standard_cnas(
    adata_query,
    regions,
    frequencies,
    copy_numbers,
    layer="counts"
)

# 'adata_gold.obs["gold_cna"]' now holds ground-truth labels per cell
```

## Benchmarking

See task2_gold_standard.ipynb for a full walk-through:
	1.	Flatten predictions (adata_query.uns['cna_summary']) and truth (adata_gold.obs['gold_cna']).
	2.	Compute precision, recall, F1 per region.
	3.	Sweep parameters (e.g. overlap_threshold, z_threshold) and display an F1 heatmap.

Launch it with:
```
jupyter notebook task2_gold_standard.ipynb
```

## Project Structure
```
Project Package/
├── cnatool/
│   ├── __init__.py
│   ├── io.py
│   ├── infer.py
│   ├── cna_call.py
│   ├── clustering.py
│   ├── comparison.py
│   ├── simulation.py
│   └── reference/         # default reference .h5ad files
├── GSE243112_with_coords.h5ad
├── PBMC_simulated_cnas_041025.h5ad
├── PBMC_gold_standard.h5ad
├── task2_gold_standard.ipynb
├── Project Run.ipynb
├── test.ipynb
├── requirements.txt
└── README.md
```
