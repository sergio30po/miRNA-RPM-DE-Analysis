# 00_packages_setup.R
# ==============================================================================
# miRNA-seq Analysis Pipeline — Package Installation and Loading
# ==============================================================================

#' @title Package Setup for miRNA-seq Analysis
#'
#' @description
#' Automates the installation (if necessary) and loading of all required R packages
#' from CRAN, Bioconductor, and GitHub for the miRNA-seq analysis pipeline.
#'
#' @details
#' This script should be run before executing any other scripts in the analysis.
#' It ensures that all required dependencies are available in the current R session.
#'
#' @sections
#' Input:
#'   - None (manual execution)
#'
#' Output:
#'   - Installs and loads required packages silently
#'   - Displays `sessionInfo()` at the end
#'
#'
#' @dependencies
#' - CRAN packages (e.g., tidyverse, ggplot2, etc.)
#' - Bioconductor packages (e.g., DESeq2, RUVSeq)
#' - GitHub packages (e.g., DANA, PoissonSeq, miRPM)
#'
#' @usage
#' Run this script once before executing the analysis pipeline to prepare the R environment.
#' 
# ==============================================================================

# 1. CRAN Packages ------------------------------------------------------------
# Standard R packages available on CRAN repository
cran_packages <- c(
  # Core Data Handling
  "tidyverse",       # Essential suite (dplyr, ggplot2, tidyr, etc.)
  "readxl",          # Excel file import/export
  "here",            # Robust file path management
  "rio",             # Unified interface for data I/O in multiple formats
  "dplyr",
  
  # Visualization
  "ggplot2",         # Advanced plotting system
  "ggrepel",         # Prevent label overlap in plots
  "plotly",          # Interactive visualizations
  "pheatmap",        # Publication-quality heatmaps
  "corrplot",        # Correlation matrix visualization
  "circlize",        # Circular layout visualizations
  "RColorBrewer",    # Color palettes for visualizations
  
  # Statistical Analysis
  "rstatix",         # Pipe-friendly statistical tests
  "car",             # Companion to Applied Regression
  "performance",     # Model performance assessment
  "factoextra",      # PCA and multivariate analysis
  "FactoMineR",      # Exploratory multivariate analysis
  "xlsx",            # Read/write Excel files
  
  # Specialized Tools
  "rgl",             # 3D visualization
  "knitr",           # Dynamic report generation
  "rmarkdown"        # R Markdown document creation
)

# 2. GitHub Packages ----------------------------------------------------------
# Packages requiring installation from GitHub repositories
github_packages <- list(
  # DANA: Data-driven Normalization Assessment
  list(repo = "LXQin/DANA", package_name = "DANA"),
  
  # PoissonSeq: Poisson-based differential expression
  list(repo = "cran/PoissonSeq", package_name = "PoissonSeq"),
  
  # miRPM: miRNA Reads Per Million normalization
  list(repo = "sergio30po/miRPM", package_name = "miRPM")
)
library(miRPM)
library(DANA)
library(PoissonSeq)

# 3. Bioconductor Packages ----------------------------------------------------
# Bioinformatics packages from Bioconductor
bioc_packages <- c(
  # RUVSeq: Remove Unwanted Variation from RNA-seq data
  "RUVSeq",
  
  # clusterProfiler: Functional enrichment analysis
  "clusterProfiler",
  
  # org.Hs.eg.db: Human genome annotation database
  "org.Hs.eg.db",
  
  # enrichplot: Visualization for functional enrichment
  "enrichplot",
  
  # Advanced heatmap customization
  "ComplexHeatmap",
  
  # miRNA DE
  "DESeq2"  
)

# 4. Setup Function -----------------------------------------------------------
#' Install and Load Required Packages
#' 
#' @description
#' Checks for installed packages, installs missing ones from appropriate sources
#' (CRAN, Bioconductor, GitHub), and loads all packages silently.
#'
#' The function:
#' 1. Identifies missing CRAN packages and installs them
#' 2. Installs specified GitHub packages using remotes
#' 3. Handles Bioconductor package installation via BiocManager
#' 4. Loads all packages with suppressed startup messages
#' 5. Provides detailed error reporting if issues occur

setup_packages <- function() {
  # Install missing CRAN packages
  missing_cran <- cran_packages[!cran_packages %in% installed.packages()]
  if(length(missing_cran) > 0) {
    message("\n[1/3] Installing ", length(missing_cran), " CRAN packages:")
    message(paste("  -", missing_cran, collapse = "\n"))
    install.packages(missing_cran)
  }
  
  # Install GitHub packages
  if(length(github_packages) > 0) {
    if(!requireNamespace("remotes", quietly = TRUE)) {
      message("\nInstalling remotes package for GitHub installations")
      install.packages("remotes")
    }
    
    missing_github <- sapply(github_packages, function(pkg) {
      !requireNamespace(pkg$package_name, quietly = TRUE)
    })
    
    if(any(missing_github)) {
      message("\n[2/3] Installing ", sum(missing_github), " GitHub packages:")
      for(pkg in github_packages[missing_github]) {
        message("  - ", pkg$package_name, " (", pkg$repo, ")")
        remotes::install_github(pkg$repo)
      }
    }
  }
  
  # Install Bioconductor packages
  if(length(bioc_packages) > 0) {
    missing_bioc <- bioc_packages[!bioc_packages %in% installed.packages()]
    if(length(missing_bioc) > 0) {
      if(!requireNamespace("BiocManager", quietly = TRUE)) {
        message("\nInstalling BiocManager for Bioconductor packages")
        install.packages("BiocManager")
      }
      message("\n[3/3] Installing ", length(missing_bioc), " Bioconductor packages:")
      message(paste("  -", missing_bioc, collapse = "\n"))
      BiocManager::install(missing_bioc)
    }
  }
  
  # Load all packages silently
  suppressPackageStartupMessages({
    # Load CRAN and Bioconductor packages
    lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)
    
    # Load GitHub packages
    lapply(sapply(github_packages, `[[`, "package_name"), library, character.only = TRUE)
  })
  
  # Success message
  message("\n=== Package setup completed successfully ===")
  message(length(c(cran_packages, bioc_packages, 
                   sapply(github_packages, `[[`, "package_name"))), 
          " packages available for use")
}

# 5. Execute Setup with Error Handling ----------------------------------------
message("Starting package setup for miRNA-seq analysis...")
start_time <- Sys.time()

tryCatch(
  {
    setup_packages()
    
    # Display session information
    message("\nSession information:")
    print(sessionInfo())
    
    # Display completion time
    end_time <- Sys.time()
    message("\nSetup completed in ", 
            round(difftime(end_time, start_time, units = "secs"), 2), 
            " seconds")
  },
  error = function(e) {
    message("\nERROR during package setup:")
    message(e$message)
    message("\nTroubleshooting steps:")
    message("1. Check internet connection")
    message("2. Verify package names are correct")
    message("3. For Bioconductor issues, try: BiocManager::install()")
    message("4. For GitHub packages, ensure you have proper access")
    message("5. Check R version compatibility (>= 4.0.0 recommended)")
    
    # Suggest specific solutions for common errors
    if(grepl("Bioconductor", e$message)) {
      message("\nBioconductor-specific solution:")
      message("Try: BiocManager::install(version = 'devel')")
    }
    
    if(grepl("GitHub", e$message)) {
      message("\nGitHub-specific solution:")
      message("1. Ensure you have 'remotes' installed")
      message("2. Check repository exists: https://github.com/sergio30po/miRPM")
    }
  }
)

# 01_data_loading_and_filtering.R
# ==============================================================================
# miRNA-seq Data Loading and Filtering Script
# ==============================================================================

#' @title Data Loading and Filtering for miRNA-seq
#'
#' @description 
#' Automates the initial data preprocessing for miRNA-seq analysis. This includes loading input files,
#' normalizing counts to RPM, filtering miRNAs by expression and read metrics, and generating structured 
#' data objects for downstream analyses.
#'
#' @details 
#' The script expects a standardized directory structure and uses input files located in the `data/` folder.
#' Normalization is performed using the `miRPM` package, and filtering is based on user-defined thresholds.
#' It concludes by constructing two S4 class objects: one with raw RPM values and another with scaled data.
#'
#' @sections
#' Input:
#'   - `data/counts.tsv`: Raw miRNA read counts (miRNAs × samples)
#'   - `data/metadata.xlsx`: Sample metadata with clinical/experimental info
#'   - `data/Metrics.xlsx`: Sequencing statistics per sample (e.g., total reads)
#'
#' Output:
#'   - `results/miRNA_ftd.csv`: Filtered RPM-normalized matrix
#'   - `results/scaled_data.csv`: Z-score scaled matrix
#'   - `results/RPM_object_scaled.rds`: S4 object with scaled expression
#'   - `results/RPM_object_raw.rds`: S4 object with raw RPM data
#'
#'
#' @dependencies 
#' - readxl: for reading `.xlsx` metadata and QC files
#' - dplyr: for manipulating metadata tables
#' - tibble: for tidy data compatibility (implicitly via dplyr)
#' - stats: for `scale()` function
#' - miRPM: custom package with functions:
#'     - `normalize_rpm()`: RPM normalization
#'     - `filter_mirnas()`: expression/read count-based filtering
#'     
#' @usage 
#' 1. Execute after running `00_packages_setup.R`.
#' 2. Ensure data files are located in `data/` folder at project root.
#' 3. Outputs are saved in `results/` for use in downstream steps.
#' 
# ==============================================================================

# 1. Utility Function ------------------------------------------------------------

create_if_missing <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory at: ", path)
  }
}

# 2. Set Paths -------------------------------------------------------------------

project_root <- here::here()
data_dir     <- file.path(project_root, "data")
results_dir  <- file.path(project_root, "results")
figures_dir  <- file.path(project_root, "figures")

create_if_missing(results_dir)
create_if_missing(figures_dir)

# Validate directory structure
if (!dir.exists(data_dir)) {
  stop("Data directory not found at: ", data_dir, "\n",
       "Expected project structure:\n",
       "project/\n",
       "├── data/\n",
       "│   ├── counts.tsv\n",
       "│   ├── metadata.xlsx\n",
       "│   └── Metrics.xlsx\n",
       "└── R/\n",
       "    └── [this script]\n",
       "└── figures/\n",
       "└── results/\n")
}

# 3. Load Data -------------------------------------------------------------------

required_files <- c(
  counts   = file.path(data_dir, "counts.tsv"),
  metadata = file.path(data_dir, "metadata.xlsx"),
  metrics  = file.path(data_dir, "Metrics.xlsx")
)

missing_files <- !file.exists(required_files)
if (any(missing_files)) {
  stop("Missing required files:\n",
       paste("-", required_files[missing_files], collapse = "\n"))
}

message("Loading data from: ", data_dir)
Metrics  <- readxl::read_excel(required_files["metrics"])
counts   <- read.delim(required_files["counts"], row.names = 1)
metadata <- readxl::read_excel(required_files["metadata"])

# 4. Process and Align Sample Names ---------------------------------------------

rownames(metadata) <- metadata$Sample

# Attempt to remove common suffix pattern (e.g., "_SM385B") if present
colnames(counts) <- gsub("_SM385B", "", colnames(counts))

# Keep only matched samples
counts <- counts[, rownames(metadata), drop = FALSE]

if (!identical(rownames(metadata), colnames(counts))) {
  mismatches <- which(rownames(metadata) != colnames(counts))
  stop("Sample mismatch between metadata and counts.\n",
       "First mismatches:\n",
       paste("  ", rownames(metadata)[mismatches][1:5], "!=", 
             colnames(counts)[mismatches][1:5], collapse = "\n"))
}

message("Successfully loaded:\n",
        "- Counts: ", nrow(counts), " miRNAs × ", ncol(counts), " samples\n",
        "- Metadata: ", nrow(metadata), " samples with ", ncol(metadata), " variables\n",
        "- Metrics: ", nrow(Metrics), " QC measurements")

# 5. RPM Normalization & Filtering -----------------------------------------------

RPM_counts <- normalize_rpm(counts, Metrics, 'Reads')

miRNA_ftd <- filter_mirnas(
  RPM_counts, 
  metadata, 
  threshold = 1, 
  min_reads = 1000, 
  threshold_comparison = ">=", 
  read_comparison = ">="
)

write.csv(miRNA_ftd, file = file.path(results_dir, "miRNA_ftd.csv"), row.names = TRUE)

raw_ftd_counts <- counts[rownames(counts) %in% rownames(miRNA_ftd), ]

# 6. Scaled Data & S4 Class Object Construction ----------------------------------

setClass(
  "RPM_object",
  slots = list(
    assay = "matrix",
    colData = "data.frame",
    rowData = "data.frame"
  )
)

metadata_df <- metadata %>%
  as.data.frame() %>%
  dplyr::select(Sample, Condition, Pathology)
rownames(metadata_df) <- metadata_df$Sample
metadata_df <- metadata_df[order(metadata_df$Condition), ]

scaled_data <- t(scale(t(miRNA_ftd)))
scaled_data <- as.matrix(scaled_data)
storage.mode(scaled_data) <- "numeric"
scaled_data <- scaled_data[, rownames(metadata_df)]

write.csv(scaled_data, file = file.path(results_dir, "scaled_data.csv"), row.names = TRUE)

RPM_object_scaled <- new(
  "RPM_object",
  assay = scaled_data,
  colData = metadata_df,
  rowData = data.frame(miRNA = rownames(scaled_data))
)

RPM_object_raw <- new(
  "RPM_object",
  assay = as.matrix(miRNA_ftd),
  colData = metadata_df,
  rowData = data.frame(miRNA = rownames(miRNA_ftd))
)

saveRDS(RPM_object_scaled, file = file.path(results_dir, "RPM_object_scaled.rds"))
saveRDS(RPM_object_raw,   file = file.path(results_dir, "RPM_object_raw.rds"))

# 7. Completion Message ----------------------------------------------------------

message("\nData loading and filtering complete.")
message("Output files saved to: ", results_dir)
message(paste("- miRNA_ftd.csv",
              "- scaled_data.csv",
              "- RPM_object_scaled.rds", 
              "- RPM_object_raw.rds", sep = "\n"))

# 02_descriptive_analysis.R
# ==============================================================================
# miRNA-seq Descriptive Analysis Script
# ==============================================================================

#' @title Descriptive analyisis before DE study
#'
#' @description
#' Performs a comprehensive descriptive analysis of miRNA expression data,
#' including PCA visualization, distribution analysis, formal statistical
#' assessment of normality and homogeneity of variances, and generation of
#' expression heatmaps.
#' 
#' @details
#' This script produces:
#' 1. PCA plots colored by experimental conditions.
#' 2. Histograms with density curves for each miRNA.
#' 3. Formal evaluation of normality (Shapiro-Wilk test) and homoscedasticity (Levene's test).
#' 4. Heatmaps of normalized expression data (Z-score).
#' 5. Statistical summary output in CSV format.
#' 
#' @sections 
#' Input Requirements:
#' - `miRNA_ftd`: filtered expression matrix (rows: miRNAs, columns: samples).
#' - `metadata`: dataframe with sample information, including columns 'Sample', 'Condition', and 'Pathology'.
#' - `scaled_data`: Z-score normalized expression matrix (columns ordered as in metadata).
#' - `figures_dir`, `results_dir`: paths to save figures and results.
#' 
#' Outputs:
#' - `figures/PCA/`: PCA plots by condition and pathology.
#' - `figures/distributions/`: histograms of expression for each miRNA.
#' - `figures/`: heatmap of expression and summary plots of normality and variance.
#' - `results/descriptive_stats.csv`: table with statistical tests and dispersion metrics.
#' 
#' @dependencies
#' Required R packages:
#' - `ggplot2`, `ggrepel` (visualization)
#' - `ComplexHeatmap`, `circlize` (heatmaps)
#' - `car` (Levene's test)
#' 
#' @usage
#' 1. Run after loading and filtering data in 'data_loading_and_filtering.R'.
#' 2. Ensure `miRNA_ftd`, `metadata`, and `scaled_data` objects are loaded in the environment.
#' 3. Outputs are automatically saved in the configured directories.
#' 
# ==============================================================================

# 0. Setup ------------------------------------------------------------------------
# Create output directories
pca_dir <- file.path(figures_dir, "PCA")
distributions_dir <- file.path(figures_dir, "distributions")
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(distributions_dir, showWarnings = FALSE, recursive = TRUE)

# 1. PCA Analysis ------------------------------------------------------------

#' @section PCA Analysis:
#' Performs dimensionality reduction and visualizes sample clustering
#' based on miRNA expression patterns. Key parameters:
#' - Data scaling: TRUE (unit variance)
#' - Components plotted: PC1 vs PC2
#' - Visualization includes:
#'   - Sample labels with repel
#'   - Group ellipses (95% CI)
#'   - Condition-specific coloring  

# Transpose and scale data
miRNA_expr_t <- t(miRNA_ftd)
pca <- prcomp(miRNA_expr_t, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by = "Sample")

# Function to create PCA plots
create_pca_plot <- function(color_var, fill_var, shape_var, 
                            color_values, fill_values, shape_values,
                            color_labels, fill_labels, shape_labels,
                            filename) {
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                          color = {{color_var}}, 
                          fill = {{fill_var}},
                          shape = {{shape_var}})) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.15, color = NA) +
    stat_ellipse(geom = "path", linewidth = 0.6) +
    geom_point(size = 3) +
    geom_text_repel(
      aes(label = Sample),  
      size = 3,
      max.overlaps = Inf,    
      show.legend = FALSE,
      box.padding = 0.2,     
      force = 1              
    ) +
    scale_color_manual(values = color_values, labels = color_labels) +
    scale_fill_manual(values = fill_values, labels = fill_labels) +
    scale_shape_manual(values = shape_values, labels = shape_labels) +
    labs(x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.line = element_line(color = "gray", linewidth = 0.2),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "grey90")
    )
  
  ggsave(file.path(pca_dir, filename), 
         plot = p, bg = "white", dpi = 300, height = 8, width = 8)
  return(p)
}

# Create PCA plots
pca_condition <- create_pca_plot(
  color_var = Condition, fill_var = Condition, shape_var = Condition,
  color_values = c("C" = "#757575", "A" = "#1874CD", "B" = "#CD2626"),
  fill_values = c("C" = "#AAAAAA", "A" = "#C6DBEF", "B" = "#FEE0D2"),
  shape_values = c("C" = 16, "A" = 17, "B" = 15),
  color_labels = c("A" = expression("LOAD non-"*italic("HTT")*" IAs"),
                   "B" = expression("LOAD "*italic("HTT")*" IAs"),
                   "C" = "Control"),
  fill_labels = c("A" = expression("LOAD non-"*italic("HTT")*" IAs"),
                  "B" = expression("LOAD "*italic("HTT")*" IAs"),
                  "C" = "Control"),
  shape_labels = c("A" = expression("LOAD non-"*italic("HTT")*" IAs"),
                   "B" = expression("LOAD "*italic("HTT")*" IAs"),
                   "C" = "Control"),
  filename = "PCA_Condition.png"
)

pca_pathology <- create_pca_plot(
  color_var = Pathology, fill_var = Pathology, shape_var = Pathology,
  color_values = c("Control" = "#757575", "AD" = "#CD2626"),
  fill_values = c("Control" = "#AAAAAA", "AD" = "#FEE0D2"),
  shape_values = c("Control" = 16, "AD" = 17),
  color_labels = c("AD" = "LOAD", "Control" = "Control"),
  fill_labels = c("AD" = "LOAD", "Control" = "Control"),
  shape_labels = c("AD" = "LOAD", "Control" = "Control"),
  filename = "PCA_Pathology.png"
)

# 2. Heatmap Visualization ----------------------------------------------------
#' @section Heatmap Analysis:
#' Visualizes z-score normalized expression patterns across samples.
#' Key features:
#' - Rows: miRNAs
#' - Columns: Samples grouped by Condition
#' - Color scale: Navy (low) to Red (high) expression
#' - Includes sample grouping annotations

message("\nGenerating expression heatmap...")

# Prepare metadata for heatmap
metadata_df <- metadata_df %>%
  mutate(
    Condition = factor(Condition, levels = c("C", "A", "B"))  #Define the order
  ) %>%
  arrange(Condition) 
scaled_data <- scaled_data[, rownames(metadata_df)]

library(ComplexHeatmap)
png(file.path(figures_dir, "expression_heatmap.png"), width = 1200, height = 1000, res = 150)
Heatmap(
  mat = scaled_data,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  column_split = metadata_df$Condition,
  show_heatmap_legend = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "gray", lwd = 0.2, fill = NA))
  }
)
dev.off()

# 3. Distribution Analysis ----------------------------------------------------
#' @section Distribution Analysis:
#' Evaluates expression distribution for each miRNA:
#' - Generates histogram with density overlay for each miRNA
#' - Saves individual plots to distributions/ subfolder
#' - Provides visual assessment of normality

# Create histograms with density plots
for (i in 1:nrow(miRNA_ftd)) {
  RPM <- as.numeric(miRNA_ftd[i, ])
  file_name <- file.path(distributions_dir, paste0(rownames(miRNA_ftd)[i], "_distribution.png"))
  
  png(filename = file_name, width = 800, height = 600)
  hist(RPM, breaks = 20, 
       main = paste("Distribution of", rownames(miRNA_ftd)[i]), 
       xlab = "RPM", probability = TRUE, col = "lightblue")
  lines(density(RPM), col = "red", lwd = 2)
  dev.off()
}

# 4. Statistical Assessment ---------------------------------------------------
#' @section Statistical Assessment:
#' Performs formal tests of data characteristics:
#' - Normality: Shapiro-Wilk test
#' - Homoscedasticity: Levene's test
#' - Dispersion: Coefficient of Variation
#' - Saves comprehensive results to CSV

# Normality tests
normality_results <- data.frame(
  miRNA = rownames(miRNA_ftd),
  p_value = sapply(1:nrow(miRNA_ftd), function(i) {
    shapiro.test(as.numeric(miRNA_ftd[i, ]))$p.value
  }),
  normality = NA_character_
)

normality_results$normality <- ifelse(normality_results$p_value > 0.05, "Normal", "Not-normal")

# Homoscedasticity tests
if (!all(metadata$Sample %in% colnames(miRNA_ftd))) {
  stop("Sample names in metadata and counts do not match!")
}

group <- metadata$Condition[match(colnames(miRNA_ftd), metadata$Sample)]
p_values_homoscedasticity <- sapply(1:nrow(miRNA_ftd), function(i) {
  leveneTest(as.numeric(miRNA_ftd[i, ]) ~ group)$`Pr(>F)`[1]
})

# Dispersion metrics
variance_values <- apply(miRNA_ftd, 1, var)
cv_values <- apply(miRNA_ftd, 1, function(x) sd(x)/mean(x))

# Combine all results
combined_results <- data.frame(
  miRNA = rownames(miRNA_ftd),
  shapiro_p = normality_results$p_value,
  normality = normality_results$normality,
  levene_p = p_values_homoscedasticity,
  homoscedasticity = ifelse(p_values_homoscedasticity > 0.05, "Yes", "No"),
  variance = variance_values,
  coefficient_of_variation = cv_values,
  dispersion_level = cut(cv_values, breaks = c(-Inf, 1, Inf), labels = c("low", "high"))
)

# Save results
write.csv(combined_results, file.path(results_dir, "descriptive_stats.csv"), row.names = FALSE)

# Create summary plot of distribution characteristics
png(file.path(figures_dir, "normality_variance_summary.png"), width = 800, height = 600)
par(mfrow = c(1,2))
hist(normality_results$p_value, main = "Normality test p-values", 
     xlab = "Shapiro-Wilk p-value", col = "lightblue")
abline(v = 0.05, col = "red", lty = 2)
hist(cv_values, main = "Coefficient of Variation", 
     xlab = "CV", col = "lightgreen")
abline(v = 1, col = "red", lty = 2)
dev.off()

# 5. Completion Message ----------------------------------------------------------
message("Descriptive analysis complete. Results saved to:")
message("- figures/PCA/: PCA plots")
message("- figures/distributions/: miRNA distribution plots")
message("- figures/: Summary heatmap and plots")
message("- results/descriptive_stats.csv: Statistical test results")

# 03_DE_analysis.R
# ==============================================================================
# miRNA Differential Expression Analysis Script
# ==============================================================================

#' @title Differential expression analysis
#'
#' @description 
#' This script performs non-parametric differential expression analysis using:
#' - Kruskal-Wallis test with FDR correction followed by Dunn’s post-hoc test (for >2 groups)
#' - Mann-Whitney U test with FDR correction (for 2-group comparisons)
#' - False Discovery Rate (FDR) correction for multiple testing
#' - DE analysis is performed using `miRPM` package functions
#'
#' @details
#' This script generates:
#' 1. Differential expression result tables
#' 2. Individual miRNA expression plots
#' 3. Interactive and static dot plots
#' 4. Heatmaps summarizing group-specific expression
#'
#' @sections 
#' Input Requirements:
#' - `miRNA_ftd`: Filtered miRNA expression matrix
#' - `metadata`: Sample metadata including Condition and Pathology columns
#' - `scaled_data`: Z-score normalized expression matrix
#' 
#' Outputs:
#' - `results/`: Differential expression results (CSV/XLSX)
#' - `figures/DE_expression_plots/`: Individual plots and heatmaps
#'
#' @dependencies
#' - `dunn.test`, `FSA`: For non-parametric tests
#' - `ggplot2`, `plotly`: For visualization
#' - `ComplexHeatmap`: For heatmap creation
#' - `miRPM`: DE analysis
#'
#' @usage
#' 1. Run after descriptive_analysis.R
#' 2. Make sure required objects are loaded
#' 3. Outputs are automatically saved to organized folders
#' 
# ==============================================================================

# 0. Setup ------------------------------------------------------------------------

# Create output directories
move_folder_recursively <- function(source_dir, destination_dir) {
  
  # Check if source exists
  if (!dir.exists(source_dir)) {
    stop("❌ Source folder does not exist: ", source_dir)
  }
  
  # Create destination if it doesn't exist
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get all items (files + folders)
  items <- list.files(source_dir, full.names = TRUE)
  
  # Copy each item (recursively if folder)
  for (item in items) {
    if (file.info(item)$isdir) {
      # If item is a directory, recursively copy
      destination_subdir <- file.path(destination_dir, basename(item))
      dir.create(destination_subdir, showWarnings = FALSE, recursive = TRUE)
      file.copy(from = list.files(item, full.names = TRUE),
                to = destination_subdir,
                recursive = TRUE,
                overwrite = TRUE)
    } else {
      # If item is a file, copy normally
      file.copy(from = item, to = destination_dir, overwrite = TRUE)
    }
  }
  
  # Force delete the original source folder
  unlink(source_dir, recursive = TRUE, force = TRUE)
  
  message("✅ All files and folders moved from '", source_dir, "' to '", destination_dir, "' and source deleted successfully!")
}

heatmap_dir <- file.path(figures_dir, "DE_Heatmaps")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)


# 1. Statistical Analysis -----------------------------------------------------

#' @section Kruskal-Wallis + FDR + Dunn Test Analysis:
#' Performs non-parametric comparison across multiple groups
#' and identifies differentially expressed miRNAs.

# Perform Kruskal-Wallis + Dunn’s Test for multiple groups
kruskal_wallis_results <- perform_statistical_tests(
  miRNA_ftd,
  metadata,
  condition = "Condition",
  output_file = file.path("Dunn_test_results.xlsx")
)

# Perform Mann-Whitney U Test for two-group comparison with FDR correction
perform_statistical_tests(
  miRNA_ftd,
  metadata,
  condition = "Pathology",
  output_file = file.path("M-W_test_results.xlsx")
)

# Organize the test results in the 'results' folder
move_folder_recursively(
  source_dir = "Tests_results",
  destination_dir = file.path(results_dir, "DE_results")
)


# 2. Individual miRNA Expression Plots ----------------------------------------

#' @section Individual Expression Plots:
#' Generate expression plots for each significant miRNA per comparison.

message("\nGenerating individual expression plots...")

plot_list <- list()

# A vs B
plot_list$A_vs_B <- individual_miRNA_plot(
  filtered_results = dunn_test_A_vs_B,
  rpm_matrix = miRNA_ftd,
  metadata = metadata,
  condition_column = "Condition",
  sample_column = "Sample",
  groups_to_include = c("A", "B"),
  condition_colors = c("A" = "blue", "B" = "red"),
  adjusted_pvalue_column = "p_value"
)

# A vs C
plot_list$A_vs_C <- individual_miRNA_plot(
  filtered_results = dunn_test_A_vs_C,
  rpm_matrix = miRNA_ftd,
  metadata = metadata,
  condition_column = "Condition",
  sample_column = "Sample",
  groups_to_include = c("A", "C"),
  condition_colors = c("A" = "blue", "C" = "green"),
  adjusted_pvalue_column = "p_value"
)

# B vs C
plot_list$B_vs_C <- individual_miRNA_plot(
  filtered_results = dunn_test_B_vs_C,
  rpm_matrix = miRNA_ftd,
  metadata = metadata,
  condition_column = "Condition",
  sample_column = "Sample",
  groups_to_include = c("B", "C"),
  condition_colors = c("B" = "red", "C" = "green"),
  adjusted_pvalue_column = "p_value"
)

# AD vs Control
plot_list$AD_vs_Control <- individual_miRNA_plot(
  filtered_results = mann_whitney_results,
  rpm_matrix = miRNA_ftd,
  metadata = metadata,
  condition_column = "Pathology",
  sample_column = "Sample",
  groups_to_include = c("AD", "Control"),
  condition_colors = c("Control" = "green", "AD" = "purple"),
  adjusted_pvalue_column = "FDR"
)

# Organize the plots in the 'figures' folder

move_folder_recursively(
  source_dir = "Individual_plots",
  destination_dir = file.path(figures_dir, "DE_individual_plots")
)
# Save individual plots
# (optional: could be saved automatically if your function does that)

# 3. Interactive Dot Plots ----------------------------------------------------

#' @section Interactive Dot Plots:
#' Visualize differentially expressed miRNAs interactively.

message("\nGenerating interactive plots...")

# Example
output <- miRNA_expression_plot(
  miRNAs_DE = dunn_test_A_vs_B,
  miRNA_ftd = miRNA_ftd,
  metadata = metadata,
  condition_column = "Condition",
  groups = c("A", "B"),
  colors = c("A" = "blue", "B" = "red"),
  output_name = "miRNA_expression_A_B",
  plot_title = "Differentially Expressed miRNAs: A vs B"
)

output$ggplot
output$interactive

# (Repeat for other comparisons A vs C, B vs C, AD vs Control)

# Organize the plots in the 'figures' folder

move_folder_recursively(
  source_dir = "Interactive_plots",
  destination_dir = file.path(figures_dir, "DE_interactive_plots")
)

# 4. Heatmaps -----------------------------------------------------------------

#' @section Heatmap Analysis:
#' Summarizes expression patterns of differentially expressed miRNAs.

generate_heatmap <- function(comparison_results, metadata, scaled_data,
                             condition_column, condition_order,
                             output_dir, output_file,
                             sorting_column) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Process miRNA results - sort by specified column
  miRNA_pval_df <- comparison_results %>%
    arrange(desc(.data[[sorting_column]])) %>%
    dplyr::select(miRNA, all_of(sorting_column))
  
  # 2. Filter and order samples
  metadata_df <- as.data.frame(metadata)
  rownames(metadata_df) <- metadata_df$Sample
  metadata_df$Sample <- NULL
  
  # Keep only samples with conditions we want
  samples_to_keep <- metadata_df[[condition_column]] %in% condition_order
  filtered_metadata <- metadata_df[samples_to_keep, ]
  
  # Order by specified condition order
  filtered_metadata <- filtered_metadata %>%
    mutate(
      !!condition_column := factor(.data[[condition_column]], levels = condition_order)
    ) %>%
    arrange(.data[[condition_column]])
  
  # Filter and order the expression data
  filtered_data <- scaled_data[miRNA_pval_df$miRNA, rownames(filtered_metadata)]
  
  # Create full output path
  full_path <- file.path(output_dir, output_file)
  
  # Create heatmap with fixed dimensions
  png(filename = full_path, width = 1200, height = 1000, res = 150)
  ht <- Heatmap(
    mat = filtered_data,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    column_split = filtered_metadata[[condition_column]],
    
    # Add cell borders
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, 
                gp = gpar(col = "gray", lwd = 0.2, fill = NA))
    }
  )
  draw(ht)
  dev.off()
  
  message(paste("Heatmap saved to:", full_path))
  return(ht)
}
# Heatmap: A vs B

generate_heatmap(
  comparison_results = dunn_test_A_vs_B,
  metadata = metadata,
  scaled_data = scaled_data,
  condition_column = "Condition",
  condition_order = c("A", "B"),
  output_file = "heatmap_A_B.png",
  output_dir = heatmap_dir,
  sorting_column = "p_value"
)

# Heatmap: A vs C
generate_heatmap(
  comparison_results = dunn_test_A_vs_C,
  metadata = metadata,
  scaled_data = scaled_data,
  condition_column = "Condition",
  condition_order = c("C", "A"),
  output_file = "heatmap_A_C.png",
  output_dir = heatmap_dir,
  sorting_column = "p_value"
)

# Heatmap: B vs C
generate_heatmap(
  comparison_results = dunn_test_B_vs_C,
  metadata = metadata,
  scaled_data = scaled_data,
  condition_column = "Condition",
  condition_order = c("C", "B"),
  output_file =  "heatmap_B_C.png",
  output_dir = heatmap_dir,
  sorting_column = "p_value"
)

# Heatmap: AD vs Control
generate_heatmap(
  comparison_results = mann_whitney_results,
  metadata = metadata,
  scaled_data = scaled_data,
  condition_column = "Pathology",
  condition_order = c("Control", "AD"),
  output_file = "heatmap_AD_Control.png",
  output_dir = heatmap_dir,
  sorting_column = "FDR"
)

# 5. Completion Message ----------------------------------------------------------
message("Differential expression analysis complete. Results saved to:")
message("- results/DE_results: Differential expression results")
message("- figures/DE_individual_plots/: Expression plots of each miRNA and comparison")
message("- figures/DE_interactive_plots/: Interactive expression plots")
message("- figures/DE_heatmaps/: Heatmaps of all the comparisons")


# 04_DANA_analysis.R
# ==============================================================================
# miRNA-seq Normalization Assessment Script
# ==============================================================================
#'
#' @description 
#' Performs a comprehensive assessment of normalization methods applied to miRNA expression data. 
#' This script includes:
#' - Application of various normalization methods (e.g., TC, UQ, TMM, DESeq2, etc.)
#' - Evaluation of normalization quality using correlation coefficients and robustness measures
#' - Visualizations comparing normalization efficiency and control behavior
#'
#' @details
#' This script generates:
#' 1. Normalized data using multiple methods (e.g., Total Count (TC), Upper Quartile (UQ), DESeq2, etc.)
#' 2. Robust correlation coefficient (CC) calculations for each normalization method
#' 3. Visualizations for distribution analysis, negative and positive control correlations, and PCA
#' 4. Summary of the normalization effectiveness across methods
#' 
#' All outputs are saved to organized subdirectories within the project folder.
#'
#' @section 
#' 
#' Input Requirements:
#' - `counts`: miRNA count matrix (raw)
#' - `RPM_counts`: miRNA count matrix normalized by Reads Per Million (RPM)
#' - `metadata`: Sample metadata with Condition column
#' - `clusters`: List of miRNA clusters (from defineClusters function)
#' - `control definitions`: Negative and positive controls (from defineControls function)
#' 
#' Outputs:
#' - `figures/DANA_figures/`: Plots including count distributions, control correlations, and PCA
#' - `results/DANA_results/`: CSV results of normalization assessment including robust correlation coefficients (cc)
#'
#' @dependencies
#' Requires the following packages:
#' - `ggplot2`, `ggrepel`: For visualization (e.g., PCA, distribution plots)
#' - `DescTools`: For robust correlation coefficient calculations
#' - `ComplexHeatmap`: For heatmap generation
#' - `car`: For statistical tests on distribution and variance
#'
#' @usage
#' 1. Run after data loading and filtering (e.g., in data_loading_and_filtering.R)
#' 2. Ensure required objects (counts, RPM_counts, metadata, clusters, controls) are loaded
#' 3. All outputs (normalized data, plots, and results) will be saved automatically to the specified directories
# ==============================================================================

# 0. Setup --------------------------------------------------------------------
# Create necessary output directories to store results and figures
dir.create(file.path(results_dir, "DANA_results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "DANA_figures"), showWarnings = FALSE, recursive = TRUE)
library(DANA)

# 1. Data Preparation ---------------------------------------------------------
# Prepare groups and counts matrices
groups <- metadata$Condition
counts <- as.matrix(counts)
RPM_counts <- as.matrix(RPM_counts)

# 2. Control Definition -------------------------------------------------------
# Define clusters and identify miRNAs excluded from clustering
clusters <- defineClusters(rownames(counts))
message(paste(length(setdiff(rownames(counts), clusters)), 
              " miRNAs excluded from clustering"))

# Define positive and negative controls with optimal thresholds
controls <- defineControls(
  counts, 
  tZero = 2,     # Lower bound for negative controls
  tPoor = 5,     # Upper bound for negative controls
  tWell = 10,    # Lower bound for positive controls
  clusters
)

# Save the control information to a CSV file
write.csv(
  data.frame(
    miRNA = c(controls$negControls, controls$posControls),
    Control_Type = c(rep("Negative", length(controls$negControls)), 
                     rep("Positive", length(controls$posControls)))
  ),
  file.path(results_dir, "DANA_results/control_miRNAs.csv"),
  row.names = FALSE
)

# 3. Diagnostic Analysis ------------------------------------------------------
# Plot distributions of miRNA counts for diagnostic purposes

# Count distribution plot
png(file.path(figures_dir, "DANA_figures/count_distribution.png"), 
    width = 800, height = 600)
plotCountHist(counts, binwidth = 0.1, 2, 5, 10, 
              title = "miRNA Count Distribution")
dev.off()

# Mean vs SD plot
png(file.path(figures_dir, "DANA_figures/mean_vs_sd.png"), 
    width = 800, height = 600)
plotMeanSD(counts, 2, 5, 10, 
           title = "Mean vs Standard Deviation")
dev.off()

# Negative control correlations plot
neg_cor <- cor(t(counts[controls$negControls, ]))
png(file.path(figures_dir, "DANA_figures/neg_control_correlations.png"), 
    width = 800, height = 600)
hist(neg_cor, breaks = 50, 
     main = "Negative Control Correlations",
     xlab = "Pearson Correlation", 
     col = "lightblue")
dev.off()

# Positive control clusters plot
posClusters <- table(clusters[controls$posControls])
png(file.path(figures_dir, "DANA_figures/pos_control_clusters.png"), 
    width = 800, height = 600)
barplot(posClusters, 
        main = "Positive Control Clusters",
        xlab = "Cluster", 
        ylab = "Number of miRNAs",
        col = "lightgreen", 
        border = "black")
dev.off()

# 4. Normalization Methods ----------------------------------------------------
# Apply multiple normalization methods to the raw count data
normalized <- applyNormalization(
  counts,
  groups,
  method = c("TC", "UQ", "median", "TMM", "DESeq", "PoissonSeq", "QN", "RUV")
)
normalized$RPM <- RPM_counts  # Add RPM normalization to the list

# 5. Robust Assessment Functions ----------------------------------------------
# Custom function to calculate robust correlation coefficient (CC)
compute.robust.cc <- function(raw, norm, posControls, clusters) {
  # Attempt 1: Cluster-based calculation
  cluster_cors <- lapply(clusters, function(clust) {
    clust_pos <- intersect(clust, posControls)
    if(length(clust_pos) >= 2) {
      raw_cor <- suppressWarnings(cor(t(raw[clust_pos, ]), use = "pairwise.complete.obs"))
      norm_cor <- suppressWarnings(cor(t(norm[clust_pos, ]), use = "pairwise.complete.obs"))
      if(!any(is.na(raw_cor))) {
        return(list(
          raw = raw_cor[upper.tri(raw_cor)],
          norm = norm_cor[upper.tri(norm_cor)]
        ))
      }
    }
    NULL
  })
  
  # Combine valid cluster correlations and calculate the final robust CC
  valid_cors <- cluster_cors[!sapply(cluster_cors, is.null)]
  if(length(valid_cors) > 0) {
    combined_raw <- unlist(lapply(valid_cors, function(x) x$raw))
    combined_norm <- unlist(lapply(valid_cors, function(x) x$norm))
    if(length(combined_raw) >= 3) {
      ccc <- tryCatch(
        DescTools::CCC(combined_raw, combined_norm)$rho.c$est,
        error = function(e) NA
      )
      if(!is.na(ccc)) return(ccc)
    }
  }
  
  # Fallback: Global correlation of positive controls
  common_pos <- intersect(posControls, rownames(norm))
  if(length(common_pos) >= 2) {
    global_cc <- suppressWarnings(
      cor(
        as.vector(raw[common_pos, ]),
        as.vector(norm[common_pos, ]),
        use = "pairwise.complete.obs"
      )
    )
    if(!is.na(global_cc)) return(global_cc)
  }
  
  # Final fallback: Global correlation between raw and normalized counts
  suppressWarnings(
    cor(
      as.vector(raw),
      as.vector(norm),
      use = "pairwise.complete.obs"
    )
  )
}

# 6. Run the Assessment -------------------------------------------------------
# Run the DANA normalization assessment with robust correlation coefficient
res <- assessNormalization(
  raw = counts,
  normalized = normalized,
  negControls = controls$negControls,
  posControls = controls$posControls,
  clusters = clusters
)

# Replace CC values with robust calculations for each method
for(method in names(normalized)) {
  res[method, "cc"] <- compute.robust.cc(
    counts,
    normalized[[method]],
    controls$posControls,
    clusters
  )
}

# 7. Save Results -------------------------------------------------------------
# Save the assessment results to CSV
write.csv(
  res,
  file.path(results_dir, "DANA_results/normalization_assessment.csv"),
  row.names = TRUE
)

# 8. Visualization -----------------------------------------------------------
# Plot DANA assessment results
png(file.path(figures_dir, "DANA_figures/DANA_assessment.png"), 
    width = 900, height = 600)
plotDANA(res)
dev.off()

# 9. Final Summary -----------------------------------------------------------
# Print final messages and paths to results and figures
message("\nDANA Assessment Complete")
message(paste("Results saved to:", file.path(results_dir, "DANA_results")))
message(paste("Figures saved to:", file.path(figures_dir, "DANA_figures")))
message("\nAssessment Results:")
print(res)
