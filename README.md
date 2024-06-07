# Project Overview

The aim of this project is to compare NK cell subsets in hot and cold tumors to potentially identify differences in their gene expression signatures.

## Datasets

Two publicly available single-cell RNA sequencing datasets were used:
- **Ovarian Tumor**: [GEO Accession GSE192898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192898)
- **Lung Cancer (Cold and Hot Tumors)**: [Lung Cancer Data](https://lungcancer.chenlulab.com/#/download)

## Methodology

The datasets were downloaded, cleaned, and preprocessed. Dimensionality reduction and clustering were performed using the `Seurat` package. Cell annotation was conducted using the `SingleR` package. Finally, differential expression and functional enrichment analyses were done to identify and interpret differences in gene expression signatures.

## Key Packages Used

- `data.table`
- `tidyverse`
- `Seurat`
- `SingleR`
- `Matrix`
- `ComplexHeatmap`
- `ggplot2`
- `clusterProfiler`

