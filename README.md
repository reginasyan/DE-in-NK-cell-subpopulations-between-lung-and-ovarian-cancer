# Project Overview

The aim of this project is to compare NK cell subsets in hot and cold tumors to potentially identify differences in their gene expression signatures.

## Datasets

Two publicly available single-cell RNA sequencing datasets were used:
- **Ovarian Tumor**: [GEO Accession GSE192898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192898)
- **Lung Cancer (Cold and Hot Tumors)**: [Lung Cancer Data](https://lungcancer.chenlulab.com/#/download)

## Methodology

1. **Data Download and Preprocessing**:
   - The datasets were downloaded, cleaned, and preprocessed.
   
2. **Dimensionality Reduction and Clustering**:
   - Performed using the `Seurat` package.
   
3. **Cell Annotation**:
   - Conducted using the `SingleR` package.
   
4. **Differential Expression and Functional Enrichment Analyses**:
   - Conducted to identify and interpret differences in gene expression signatures.

## Key Packages Used

- `data.table`
- `tidyverse`
- `Seurat`
- `SingleR`
- `Matrix`
- `ComplexHeatmap`
- `ggplot2`
- `clusterProfiler`

