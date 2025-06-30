# Connecting HTT intermediate alleles and microRNA dysregulation to enhanced tauopathy in Late-Onset Alzheimer's Disease

Juan Castilla-Silgado ^1^,^2^,^6^ *, Sergio Perez-Oliveira ^2^,^4^ *, Paola Pinto-Hernandez ^1^,^2^, Manuel Fernandez-Sanjurjo ^1^,^2^, Maria Daniela Corte-Torres ^2^,^7^, Eduardo Iglesias-Gutierrez ^1^,^2^, Manuel Menendez-Gonzalez ^2^,^3^,^5^, Cristina Tomas-Zapico ^1^,^2^,#,*, Victoria Alvarez ^2^,^4^,#

1 Department of Functional Biology (Physiology), University of Oviedo, 33006, Oviedo, Spain.
2 Instituto de Investigación Sanitaria del Principado de Asturias (ISPA), 33011, Oviedo, Spain.
3 Department of Neurology, Hospital Universitario Central de Asturias (HUCA), 33011, Oviedo, Spain. 
4 Laboratory of Genetics Hospital Universitario central de Asturias, 33011, Oviedo, Spain.
5 Department of Medicine, University of Oviedo, 33006, Oviedo, Spain.
6 Asociación Parkinson Asturias 33011, Oviedo, Spain.
7 Biobank of Principado de Asturias, Hospital Universitario Central de Asturias (HUCA), 33011, Oviedo, Spain.

**Abstract**
Late-onset Alzheimer´s disease (LOAD) is a heterogenous disorder influenced by genetic factors, where we have previously described intermediate alleles (IAs) in the huntingtin (HTT) gene as potential modifiers. We investigated the impact of HTT IAs on LOAD progression and neuropathology using a comprehensive approach, genotyping HTT CAG repeats in 323 LOAD patients and 335 healthy controls, and performing multi-omic analyses on caudate nucleus samples. Clinically, HTT IAs carriers exhibited decreased survival after disease onset, suggesting accelerated progression. Histopathologically, while LOAD patients showed increased soluble HTT levels and altered tau pathology compared to controls, these changes were consistently and markedly exacerbated in HTT IA carriers, characterized by heightened diffuse HTT immunoreactivity, a pronounced tau 3R isoform imbalance, and increased 3R tau-enriched ghost tangles. Mechanistically, this pathological exacerbation was supported by alterations in key splicing factors, including decreased SRSF6 and increased FUS-SFPQ complex formation. Our pioneering microRNA (miRNA) profiling in the caudate nucleus revealed a LOAD-associated miRNA dysregulation that was significantly amplified in HTT IA carriers. In silico analyses, supported by network modeling and direct target validation, demonstrated that these altered miRNAs specifically target components of the nuclear spliceosome machinery, including MAPT and HTT genes, suggesting a direct link to the observed tau 3R/4R imbalance. These findings underscore that HTT IAs act as critical accelerators of LOAD progression through a complex interplay involving miRNA-mediated dysregulation of splicing. Identifying HTT IAs through routine genetic screening offers a practical, non-invasive biomarker for patient stratification, paving the way for personalized therapeutic strategies in LOAD.

---

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

│     ├── 04_DANA_analysis.R
  
│     └── DE_pipeline.R



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
2. Run [00_packages_setup.R](R/00_packages_setup.R) to prepare the R environment.  
3. Sequentially run the scripts from  [01_data_loading_and_filtering.R](R/01_data_loading_and_filtering.R) through [04_DANA_analysis.R](R/04_DANA_analysis.R).  
4. Inspect the outputs in the `results/` and `figures/` directories.  

---

## Dependencies

  **CRAN packages**: tidyverse, ggplot2, dplyr, readxl, tibble, stats, ggrepel, car, DescTools, dunn.test, FSA, plotly, ComplexHeatmap, circlize

  **Bioconductor packages:**
  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)  
  - [RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
  
  **GitHub packages:**
  - [DANA](https://github.com/LXQin/DANA)  
  - [PoissonSeq](https://github.com/lsy1056/PoissonSeq)  
  - [miRPM (custom)](https://github.com/sergio30po/miRPM)  


---

## Contact

For questions or contributions, please contact Sergio Pérez Oliveira(sergio30po@gmail.com).


*This pipeline is developed as part of a doctoral thesis project analyzing miRNA expression and differential expression in neurological disorders.*

