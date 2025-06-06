# miRNA RPM Differential Expression Analysis Pipeline

This repository contains a comprehensive R-based pipeline for the analysis of miRNA-seq data, focusing on normalization, filtering, descriptive statistics, differential expression, and normalization assessment. The workflow is modularized into separate scripts for ease of use and reproducibility.

---

## Repository Structure

├── data/  

│     ├── counts/

│     ├── metadata/

│     └── Metrics/

└── R/

│     ├── 00_packages_setup.R

│     ├── 01_data_loading_and_filtering.R

│     ├── 02_descriptive_analysis.R

│     ├── 03_DE_analysis.R
  
│     └── 04_DANA_analysis.R



## Script Descriptions

### [00_packages_setup.R](R/00_packages_setup.R)
**Purpose:** Automates the installation and loading of all required R packages from CRAN, Bioconductor, and GitHub to prepare the analysis environment.  
**Details:** Run this script first to ensure all dependencies are installed and loaded before executing any other scripts.  
**Output:** Displays session info and silently loads packages.

---

### [01_data_loading_and_filtering.R](R/01_data_loading_and_filtering.R) 
**Purpose:** Loads raw miRNA-seq data, normalizes counts to RPM, filters miRNAs based on expression/read thresholds, and saves processed data for downstream analyses.  
**Inputs:** Raw count matrix (`counts.tsv`), sample metadata (`metadata.xlsx`), sequencing metrics (`Metrics.xlsx`).  
**Outputs:** Filtered RPM matrix, scaled data, and serialized S4 objects saved under `results/`.

---

### [02_descriptive_analysis.R](R/02_descriptive_analysis.R)  
**Purpose:** Performs exploratory data analysis including PCA, distribution plots, normality and variance tests, and heatmaps of miRNA expression.  
**Inputs:** Filtered and scaled expression data, sample metadata.  
**Outputs:** PCA plots, histograms, heatmaps, and statistical summaries saved under `figures/` and `results/`.

---

### [03_DE_analysis.R](R/03_DE_analysis.R)  
**Purpose:** Executes non-parametric differential expression tests (Kruskal-Wallis, Dunn’s post-hoc, Mann-Whitney U) with FDR correction, generates plots and result tables.  
**Inputs:** Filtered miRNA expression, metadata, scaled data.  
**Outputs:** DE results tables and visualization plots saved in `results/` and `figures/`.

---

### [04_DANA_analysis.R](R/04_DANA_analysis.R)  
**Purpose:** Assesses and compares normalization methods using correlation metrics and visualizations (PCA, distributions, control behavior).  
**Inputs:** Raw counts, RPM-normalized counts, metadata, miRNA clusters, control definitions.  
**Outputs:** Normalized data, robust correlation coefficients, and comparative plots saved under `results/DANA_results/` and `figures/DANA_figures/`.

---
### [DE_pipeline.R](R/DE_pipeline.R)  
This script integrates all the major steps of the miRNA differential expression analysis into a single, streamlined pipeline. It automates the workflow from raw count data loading, normalization, filtering, through descriptive statistics, differential expression testing, and result visualization. By running this script, users can perform the entire analysis process with minimal manual intervention, ensuring reproducibility and consistency across datasets.

---

## Usage Instructions

1. Clone the repository to your local machine.  
2. Run `00_packages_setup.R` to prepare the R environment.  
3. Sequentially run the scripts from `01_data_loading_and_filtering.R` through `04_DANA_analysis.R`.  
4. Inspect the outputs in the `results/` and `figures/` directories.  

---

## Dependencies

- CRAN packages: tidyverse, ggplot2, dplyr, readxl, tibble, stats, ggrepel, car, DescTools, dunn.test, FSA, plotly, ComplexHeatmap, circlize  
- Bioconductor packages: DESeq2, RUVSeq  
- GitHub packages: DANA, PoissonSeq, miRPM (custom)  

---

## Contact

For questions or contributions, please contact Sergio Pérez Oliveira(sergio30po@gmail.com).


*This pipeline is developed as part of a doctoral thesis project analyzing miRNA expression and differential expression in neurological disorders.*

