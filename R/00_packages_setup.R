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