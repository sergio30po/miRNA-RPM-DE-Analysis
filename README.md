# Connecting HTT intermediate alleles and microRNA dysregulation to enhanced tauopathy in Late-Onset Alzheimer's Disease

Juan Castilla-Silgado<sup>1,2,6</sup> *, Sergio Perez-Oliveira<sup>2,4</sup> *, Paola Pinto-Hernandez<sup>1,2</sup>, Manuel Fernandez-Sanjurjo<sup>1,2</sup>, Maria Daniela Corte-Torres<sup>2,7</sup>,Olaya Fernández-Álvarez<sup>2,4,6</sup>, Eduardo Iglesias-Gutierrez<sup>1,2</sup>, Manuel Menendez-Gonzalez<sup>2,3,5</sup>, Cristina Tomas-Zapico<sup>1,2</sup><sup>#</sup>, Victoria Alvarez<sup>2,4</sup><sup>#</sup>

<sup>1</sup> Department of Functional Biology (Physiology), University of Oviedo, 33006, Oviedo, Spain  
<sup>2</sup> Instituto de Investigación Sanitaria del Principado de Asturias (ISPA), 33011, Oviedo, Spain  
<sup>3</sup> Department of Neurology, Hospital Universitario Central de Asturias (HUCA), 33011, Oviedo, Spain  
<sup>4</sup> Laboratory of Genetics, Hospital Universitario Central de Asturias (HUCA), 33011, Oviedo, Spain  
<sup>5</sup> Department of Medicine, University of Oviedo, 33006, Oviedo, Spain  
<sup>6</sup> Asociación Parkinson Asturias, 33011, Oviedo, Spain  
<sup>7</sup> Biobank of Principado de Asturias, Hospital Universitario Central de Asturias (HUCA), 33011, Oviedo, Spain  

**\*** These authors contributed equally  
**#** Joint senior authors

**Abstract**

**Background:** Late-onset Alzheimer´s disease (LOAD) is a heterogeneous disorder influenced by complex genetic factors. We previously described intermediate alleles (IAs; 27-35 CAG repeats) in the huntingtin (HTT) gene as potential modifiers in around 6% of AD population. The caudate nucleus, the most affected region in Huntington's disease, is highly sensitive to these HTT CAG repeats. We hypothesized that HTT IAs induce gene expression deregulation, including altered microRNA (miRNA) profiles, leading to altered disease progression.

**Methods:** We investigated the impact of HTT IAs on LOAD progression by genotyping HTT CAG repeats in a cohort of 323 LOAD patients and 335 healthy controls. Comprehensive histopathological and molecular analyses were performed on caudate nucleus samples from a matched subcohort (6 healthy controls, 14 LOAD non-HTT IA carriers, and 13 LOAD HTT IA carriers). 

**Results:** HTT IAs carriers patients exhibited decreased survival after disease onset compared to non-carriers. Histopathologically, while LOAD patients showed increased soluble HTT levels and altered tau pathology compared to controls, these changes were consistently and markedly exacerbated in HTT IA carriers. This phenotype was characterized by heightened diffuse HTT immunoreactivity and an advanced maturation of tau pathology, specifically a pronounced increase in neuronal tau 3R burden and 3R tau-enriched ghost tangles. Interestingly, this pathological state was associated with alterations in key splicing factors, including decreased SRSF6 levels and increased nuclear FUS-SFPQ complex assembly. Analysis of microRNA (miRNA) profiling in the caudate nucleus revealed that the LOAD-associated miRNA dysregulation was significantly amplified in HTT IA carriers, identifying a signature of five miRNAs (miR-100-5p, miR-218-5p, miR-27b-3p, miR-487-3p, and miR-9-3p). In silico modeling and target validation demonstrated that these miRNAs target components of the nuclear spliceosome machinery, such as SRSF family, along with MAPT and HTT genes, suggesting a direct link to the observed tauopathy. 

**Conclusions:** Our findings underscore that HTT IAs as critical modifiers in LOAD progression through an miRNA-mediated dysregulation of splicing and proteostasis. Thus, identifying HTT IAs through routine blood genetic screening offers a practical, non-invasive biomarker for patient stratification, taking a step forward to personalized therapeutic strategies in LOAD."

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

