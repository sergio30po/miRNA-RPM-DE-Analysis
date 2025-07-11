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
