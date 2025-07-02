# ==============================================================================
# SINGLE CELL RNA-SEQ ANALYSIS: COVID-19 vs CONTROL LUNG SAMPLES
# Script 2: Metadata Alignment and Quality Control
# ==============================================================================
# Purpose: 
# 1. Align published metadata (cell types) to our merged Seurat object
# 2. Perform quality control filtering
# 3. Save filtered object for downstream analysis
# ==============================================================================

# Load required packages ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)  # For plot_annotation function

print("Starting metadata alignment and QC...")

# SECTION 1: LOAD DATA ----
print("Loading merged Seurat object...")
seurat_obj <- readRDS("results/merged_C51_L01.rds")

print(paste("Loaded object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes"))

# SECTION 2: LOAD AND PROCESS METADATA ----
print("Loading published metadata...")

# Load metadata file with cell type annotations
meta_data <- read.delim("data/GSE171524_lung_metaData.txt.gz", 
                        sep = "\t", header = TRUE)

# Remove artifact header row if present (common in GEO downloads)
meta_data <- meta_data[-1, ]

# Rename first column for clarity
colnames(meta_data)[1] <- "barcode"

print(paste("Loaded metadata for", nrow(meta_data), "cells"))

# SECTION 3: BARCODE ALIGNMENT (CRITICAL STEP) ----
print("Aligning barcodes between Seurat object and metadata...")

# Filter metadata for our specific samples (C51 control, L01 COVID)
meta_cleaned <- meta_data %>%
  filter(biosample_id %in% c("C51ctr", "L01cov")) %>%
  mutate(
    # Step 1: Remove sample-specific suffixes from original barcodes
    barcode_stripped = case_when(
      biosample_id == "C51ctr" ~ str_remove(barcode, "-1_1$"),
      biosample_id == "L01cov" ~ str_remove(barcode, "-1_8$"),
      TRUE ~ barcode
    ),
    # Step 2: Create new barcodes matching Seurat object format
    clean_barcode = case_when(
      biosample_id == "C51ctr" ~ paste0("C51_", barcode_stripped, ".1_1"),
      biosample_id == "L01cov" ~ paste0("L01_", barcode_stripped, ".1_8"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(clean_barcode))

# Standardize condition names (COVID-19 -> COVID for consistency)
meta_cleaned$group <- ifelse(meta_cleaned$group == "COVID-19", "COVID", meta_cleaned$group)

print(paste("Cleaned metadata for", nrow(meta_cleaned), "cells"))

# SECTION 4: VALIDATE ALIGNMENT ----
# Check how many barcodes match between datasets
matched_barcodes <- colnames(seurat_obj)[colnames(seurat_obj) %in% meta_cleaned$clean_barcode]
alignment_rate <- length(matched_barcodes) / ncol(seurat_obj) * 100

cat("✅ Barcode alignment results:\n")
cat("   Matched:", length(matched_barcodes), "out of", ncol(seurat_obj), "cells\n")
cat("   Alignment rate:", round(alignment_rate, 1), "%\n")

if(alignment_rate < 80) {
  warning("Low alignment rate - check barcode formatting!")
}

# SECTION 5: ADD METADATA TO SEURAT OBJECT ----
print("Adding metadata columns to Seurat object...")

# Align metadata rows to match Seurat object column order
meta_aligned <- meta_cleaned[match(colnames(seurat_obj), meta_cleaned$clean_barcode), ]

# Safely add metadata columns (check if they exist first)
metadata_columns <- c("cell_type_main", "cell_type_intermediate", "cell_type_fine", 
                      "group", "donor_id")

for(col in metadata_columns) {
  if(col %in% colnames(meta_aligned)) {
    if(col == "cell_type_main") seurat_obj$celltype_main <- meta_aligned[[col]]
    if(col == "cell_type_intermediate") seurat_obj$celltype_intermediate <- meta_aligned[[col]]
    if(col == "cell_type_fine") seurat_obj$celltype_fine <- meta_aligned[[col]]
    if(col == "group") seurat_obj$published_condition <- meta_aligned[[col]]
    if(col == "donor_id") seurat_obj$donor_id <- meta_aligned[[col]]
    
    cat("   Added:", col, "\n")
  }
}

# SECTION 6: VALIDATE METADATA INTEGRATION ----
print("Validating metadata integration...")

cat("\n📊 Published condition distribution:\n")
print(table(seurat_obj$published_condition, useNA = "ifany"))

cat("\n📊 Original condition distribution:\n")
print(table(seurat_obj$condition, useNA = "ifany"))

if("celltype_intermediate" %in% colnames(seurat_obj@meta.data)) {
  cat("\n📊 Cell types by condition:\n")
  celltype_table <- table(seurat_obj$celltype_intermediate, seurat_obj$condition, useNA = "ifany")
  print(celltype_table)
}

# SECTION 7: QUALITY CONTROL CALCULATIONS ----
print("Calculating QC metrics...")

# Calculate mitochondrial gene percentage (indicator of cell stress/death)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Calculate ribosomal gene percentage (optional - indicator of metabolic activity)
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

cat("\n📊 QC Summary Statistics:\n")
cat("Features per cell (genes):\n")
print(summary(seurat_obj$nFeature_RNA))
cat("\nCounts per cell (UMIs):\n")
print(summary(seurat_obj$nCount_RNA))
cat("\nMitochondrial percentage:\n")
print(summary(seurat_obj$percent.mt))

# SECTION 8: QC VISUALIZATIONS ----
print("Creating QC plots...")

# Violin plots for key QC metrics
qc_violin <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "condition",
  pt.size = 0.1,
  ncol = 3
)

# Add title using patchwork
qc_violin <- qc_violin + plot_annotation(title = "Quality Control Metrics by Condition")

# Scatter plots to identify outliers
qc_scatter1 <- FeatureScatter(seurat_obj, 
                              feature1 = "nCount_RNA", 
                              feature2 = "percent.mt",
                              group.by = "condition") +
  ggtitle("Total UMIs vs Mitochondrial %")

qc_scatter2 <- FeatureScatter(seurat_obj, 
                              feature1 = "nCount_RNA", 
                              feature2 = "nFeature_RNA",
                              group.by = "condition") +
  ggtitle("Total UMIs vs Number of Genes")

# Save plots
ggsave("figures/qc_violin_plots.png", qc_violin, width = 12, height = 4, dpi = 300)
ggsave("figures/qc_scatter_mito.png", qc_scatter1, width = 8, height = 6, dpi = 300)
ggsave("figures/qc_scatter_features.png", qc_scatter2, width = 8, height = 6, dpi = 300)

# SECTION 9: QUALITY FILTERING ----
print("Applying quality filters...")

# Define filtering criteria (adjust based on your data)
min_features <- 200    # Minimum genes per cell
max_features <- 6000   # Maximum genes per cell (remove potential doublets)
max_mito <- 5          # Maximum mitochondrial % (remove dying cells)

cat("\nFiltering criteria:\n")
cat("  Minimum features per cell:", min_features, "\n")
cat("  Maximum features per cell:", max_features, "\n")
cat("  Maximum mitochondrial %:", max_mito, "\n")

# Apply filters
seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA > min_features & 
    nFeature_RNA < max_features & 
    percent.mt < max_mito
)

# Report filtering results
cells_removed <- ncol(seurat_obj) - ncol(seurat_obj_filtered)
removal_rate <- cells_removed / ncol(seurat_obj) * 100

cat("\n📊 Filtering Results:\n")
cat("  Before filtering:", ncol(seurat_obj), "cells\n")
cat("  After filtering:", ncol(seurat_obj_filtered), "cells\n")
cat("  Cells removed:", cells_removed, "(", round(removal_rate, 1), "%)\n")

# Check condition balance after filtering
cat("\n📊 Condition balance after filtering:\n")
print(table(seurat_obj_filtered$condition))

# SECTION 10: SAVE RESULTS ----
print("Saving filtered object...")

# Save filtered Seurat object
saveRDS(seurat_obj_filtered, "results/filtered_C51_L01.rds")

# Save QC summary
qc_summary <- data.frame(
  Metric = c("Total_cells_before", "Total_cells_after", "Cells_removed", "Removal_rate_%"),
  Value = c(ncol(seurat_obj), ncol(seurat_obj_filtered), cells_removed, round(removal_rate, 1))
)
write.csv(qc_summary, "results/qc_summary.csv", row.names = FALSE)

# Clean up memory
rm(meta_data, meta_cleaned, meta_aligned)
gc()

cat("\n script 2 completed successfully!\n")
cat("   Filtered object saved to: results/filtered_C51_L01.rds\n")
cat("   QC plots saved to: figures/\n")
cat("   Ready for normalization and clustering!\n")
