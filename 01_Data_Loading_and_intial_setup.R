
# ==============================================================================
# SINGLE CELL RNA-SEQ ANALYSIS: COVID-19 vs CONTROL LUNG SAMPLES
# Script 1: Data Loading and Initial Setup
# ==============================================================================

# Load required libraries ----
library(Seurat)
library(dplyr)

# Set up directories ----
if(!dir.exists("results")) dir.create("results")
if(!dir.exists("figures")) dir.create("figures")

# SECTION 1: LOAD CONTROL SAMPLE ----
print("Loading control sample...")

# Load control sample data
ctrl_counts <- read.csv(gzfile("data/GSM5226574_C51ctr_raw_counts.csv.gz"), 
                        row.names = 1)

# Quick data exploration
print(paste("Control sample dimensions:", nrow(ctrl_counts), "genes x", ncol(ctrl_counts), "cells"))
print("First few genes and cells:")
head(ctrl_counts[1:5, 1:5])

# Create Seurat object with quality filters
# min.cells = 3: only keep genes detected in at least 3 cells
# min.features = 200: only keep cells with at least 200 detected genes
ctrl <- CreateSeuratObject(counts = ctrl_counts, 
                           project = "COVID_Lung",
                           min.cells = 3,
                           min.features = 200)

# Add sample metadata
ctrl$sample_id <- "C51"
ctrl$condition <- "Control"

print(paste("Control after filtering:", ncol(ctrl), "cells,", nrow(ctrl), "genes"))

# SECTION 2: LOAD COVID SAMPLE ----
print("Loading COVID sample...")

# Load COVID sample data
cov_counts <- read.csv(gzfile("data/GSM5226581_L01cov_raw_counts.csv.gz"), 
                       row.names = 1)

# Quick data exploration
print(paste("COVID sample dimensions:", nrow(cov_counts), "genes x", ncol(cov_counts), "cells"))

# Create Seurat object with same filters for consistency
cov <- CreateSeuratObject(counts = cov_counts, 
                          project = "COVID_Lung", 
                          min.cells = 3, 
                          min.features = 200)

# Add sample metadata
cov$sample_id <- "L01"
cov$condition <- "COVID"

print(paste("COVID after filtering:", ncol(cov), "cells,", nrow(cov), "genes"))

# SECTION 3: MERGE SAMPLES ----
print("Merging samples...")

# Merge both samples
# add.cell.ids adds prefixes to distinguish cells from different samples
combined <- merge(ctrl, y = cov, 
                  add.cell.ids = c("C51", "L01"), 
                  project = "COVID_Lung")

# Check the merge
print("Sample distribution after merging:")
print(table(combined$condition))
print(table(combined$sample_id))

# Quick summary
print(paste("Final merged object:", ncol(combined), "cells,", nrow(combined), "genes"))

# SECTION 4: SAVE PROCESSED DATA ----
print("Saving merged object...")

# Save merged object for downstream analysis
saveRDS(combined, file = "results/merged_C51_L01.rds")

print("Data loading and setup complete!")

# SECTION 5: BASIC QUALITY METRICS ----
# Calculate some basic metrics to understand the data
combined$nCount_RNA <- combined@meta.data$nCount_RNA  # Total UMI counts per cell
combined$nFeature_RNA <- combined@meta.data$nFeature_RNA  # Number of genes per cell

# Basic statistics
print("Basic quality metrics:")
print("Control cells - Median genes per cell:")
print(median(combined@meta.data[combined$condition == "Control", "nFeature_RNA"]))
print("COVID cells - Median genes per cell:")
print(median(combined@meta.data[combined$condition == "COVID", "nFeature_RNA"]))

# Clean up large objects to save memory
rm(ctrl_counts, cov_counts, ctrl, cov)
gc()  # Garbage collection

