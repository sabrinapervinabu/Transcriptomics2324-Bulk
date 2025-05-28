# Bulk RNA-Seq Differential Expression Analysis (GTEx Data)

This repository contains an R script for differential gene expression analysis using bulk RNA-seq data from the GTEx project. The analysis compares gene expression in **brain**, **heart**, and **kidney** tissues, leveraging `recount3` and the `edgeR` pipeline.

---

## Contents

- `bulk trascrittomica.R`: R script for:
  - Loading and processing GTEx datasets
  - Sample quality control and selection
  - Data normalization (TMM)
  - Differential expression analysis with edgeR
  - Identification of tissue-specific genes
  - Visualization and export of DE genes

---

## Requirements

Make sure to install the required R and Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "recount3", "recount", "GenomicRanges", "limma", "edgeR", "DESeq2",
  "regionReport", "clusterProfiler", "org.Hs.eg.db", "gplots",
  "derfinder", "GenomicState", "bumphunter", "derfinderPlot", "sessioninfo"
))
```
## How to Run

Clone or download this repository.
Place the following GTEx-derived files in your working directory:
- rse_brain.RDS
- rse_heart.RDS
- rse_kidney.RDS
Run the script in R or RStudio:
```r
source("bulk trascrittomica.R")
```
## Key Steps
Sample Quality Control
Filters samples based on:
- RNA Integrity Number (RIN ≥ 6)
- % rRNA content (≤ 10%)
- % uniquely mapped reads (≥ 85%)
### Normalization
- Uses TMM (Trimmed Mean of M-values) normalization with edgeR.
- Pre/post normalization comparison using cpm() and boxplot().
### Differential Expression
Linear model design using model.matrix(~0 + group).
Pairwise comparisons:
- Brain vs Heart
- Brain vs Kidney
- Heart vs Kidney
### Visualization
- Multi-Dimensional Scaling (MDS) plots
- Biological Coefficient of Variation (BCV) plot
- TPM expression boxplots for key DEGs
### DEG Identification
Genes are filtered by:
- FDR < 0.01
- |logFC| ≥ 1
Gene lists are extracted for:
- Up-regulated in brain vs (heart + kidney)
- Up-regulated in heart vs (brain + kidney)
- Up-regulated in kidney vs (brain + heart)

## Output
- resultsHB.txt, resultsKB.txt, resultsKH.txt: full DE tables
- CSVs of the top 500 upregulated genes per tissue:
    - brain_upregulated_first_500.csv
    - heart_upregulated_first_500.csv
    - kidney_upregulated_first_500.csv

## References
- GTEx Project
- recount3 documentation
- edgeR User Guide
- Professor's resources:
    - http://159.149.160.56/GT_2022_GTEX/recount_DEG.html
    - http://159.149.160.56/GT_2022_GTEX/recount3.html
