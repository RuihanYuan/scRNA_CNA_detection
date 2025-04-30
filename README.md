# CNA Inferencer for scRNA-seq Data

This Python package provides a method to infer **Copy Number Alterations (CNAs)** from single-cell RNA-sequencing (scRNA-seq) data. CNAs are identified as genomic gains or losses and assigned to individual cells or cell groups using a combination of gene expression changes and optional allele imbalance. The tool is designed to work with `AnnData` objects and outputs structured CNA annotations for visualization and downstream analysis.

---

## Project Tasks

This tool addresses the following project objectives:

### Task 1: CNA Inference
- Input: AnnData object with raw gene expression counts
- Output: List of CNAs defined by genomic region and type ('gain' or 'loss'), and per-cell CNA assignments
- Compatible with visualization and downstream phylogenetic analysis

### Task 2: Performance Assessment
- Evaluate method accuracy on test datasets using metrics such as precision, recall, and AUC
- Explore effects of read depth using downsampling

### Task 3: Application to PSC scRNA-seq
- Apply this method to selected pluripotent stem cell (PSC) datasets to detect previously known and potentially novel CNAs

### Optional (Extra Credit)
- Simulate gold-standard data for benchmarking (Task 2B)
- Predict the downstream functional impact of CNAs on PSCs (Task 4)

---

## Installation

Clone the repository:
```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME# scRNA_CNA_detection
