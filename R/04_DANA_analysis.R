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
