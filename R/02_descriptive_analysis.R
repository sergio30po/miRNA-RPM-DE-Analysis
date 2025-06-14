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