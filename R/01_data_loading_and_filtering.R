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