# ===============================================================================
# ENHANCED SINGLE-CELL RNA-SEQ ANALYSIS PIPELINE
# ===============================================================================
# Purpose: Complete analysis pipeline for scRNA-seq data with UMAP visualization,
#          clustering, marker identification, and cell type annotation
# ===============================================================================

# --------------------------------
# PACKAGE LOADING & SETUP
# --------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(cowplot)
})

# Set random seed for reproducibility
set.seed(42)

# Create output directories
output_dirs <- c("results", "figures", "tables")
lapply(output_dirs, function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
})

# --------------------------------
# ANALYSIS PARAMETERS
# --------------------------------
PARAMS <- list(
  # Input/Output
  input_file = "results/filtered_C51_L01.rds",
  sample_name = "C51_L01",
  
  # Analysis parameters
  pca_dims = 30,
  neighbors_dims = 1:20,
  umap_dims = 1:20,
  cluster_resolution = 0.5,
  
  # Marker identification
  min_pct = 0.25,           # Increased for more stringent filtering
  logfc_threshold = 0.25,   # Increased for stronger markers
  top_markers_n = 5,        # More markers per cluster
  
  # Plotting
  point_size = 0.4,
  plot_width = 12,
  plot_height = 10,
  dpi = 300
)

cat("🔧 Analysis parameters loaded\n")
print(PARAMS)

# --------------------------------
# UTILITY FUNCTIONS
# --------------------------------

# Function to print section headers
print_section <- function(title) {
  cat("\n", rep("=", 60), "\n")
  cat( title, "\n")
  cat(rep("=", 60), "\n")
}

# Function to save plots with consistent formatting and size limits
save_plot <- function(plot, filename, width = PARAMS$plot_width, height = PARAMS$plot_height) {
  # Ensure dimensions don't exceed reasonable limits
  width <- min(width, 25)   # Max 25 inches wide
  height <- min(height, 25) # Max 25 inches tall
  
  full_path <- file.path("figures", filename)
  ggsave(full_path, plot, width = width, height = height, dpi = PARAMS$dpi, limitsize = FALSE)
  cat(" Plot saved:", full_path, "(", width, "x", height, "inches )\n")
}

# Function to create enhanced feature plots
create_feature_plot <- function(seurat_obj, gene, title_suffix = "") {
  FeaturePlot(
    seurat_obj, 
    features = gene, 
    reduction = "umap",
    pt.size = PARAMS$point_size,
    cols = c("lightgrey", "red"),
    order = TRUE  # Plot high-expressing cells on top
  ) +
    ggtitle(paste0(gene, title_suffix)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10)
    )
}

# --------------------------------
# DATA LOADING & VALIDATION
# --------------------------------
print_section("DATA LOADING & VALIDATION")

if (!file.exists(PARAMS$input_file)) {
  stop(" Input file not found: ", PARAMS$input_file)
}

seurat_obj <- readRDS(PARAMS$input_file)
cat(" Loaded Seurat object from:", PARAMS$input_file, "\n")

# Print basic object information
cat(" Object Summary:\n")
cat("   - Cells:", ncol(seurat_obj), "\n")
cat("   - Genes:", nrow(seurat_obj), "\n")
cat("   - Assays:", names(seurat_obj@assays), "\n")

# Check for existing cell type annotations
celltype_cols <- grep("celltype", colnames(seurat_obj@meta.data), value = TRUE)
if (length(celltype_cols) > 0) {
  cat("   - Available cell type annotations:", paste(celltype_cols, collapse = ", "), "\n")
} else {
  cat("   - No existing cell type annotations found\n")
}

# --------------------------------
# 🔬 STANDARD PREPROCESSING PIPELINE
# --------------------------------
print_section("STANDARD PREPROCESSING PIPELINE")

# Check if normalization already done
if (!"LogNormalize" %in% seurat_obj@commands) {
  cat("Running normalization...\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
} else {
  cat(" Data already normalized\n")
}

# Find variable features
cat("🔄 Finding variable features...\n")
seurat_obj <- FindVariableFeatures(
  seurat_obj, 
  selection.method = "vst", 
  nfeatures = 2000,
  verbose = FALSE
)

# Plot variable features
var_plot <- VariableFeaturePlot(seurat_obj)
top10_var <- head(VariableFeatures(seurat_obj), 10)
var_plot_labeled <- LabelPoints(
  plot = var_plot, 
  points = top10_var, 
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)
save_plot(var_plot_labeled, "variable_features.png", width = 10, height = 8)

# Scale data
cat(" Scaling data...\n")
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

# --------------------------------
# DIMENSIONALITY REDUCTION
# --------------------------------
print_section("DIMENSIONALITY REDUCTION")

# Run PCA
cat(" Running PCA...\n")
seurat_obj <- RunPCA(
  seurat_obj, 
  features = VariableFeatures(object = seurat_obj),
  npcs = PARAMS$pca_dims, 
  verbose = FALSE
)

# Plot PCA
pca_plot <- DimPlot(seurat_obj, reduction = "pca", pt.size = PARAMS$point_size) +
  ggtitle("PCA: Principal Components 1 & 2")
save_plot(pca_plot, "pca_plot.png", width = 8, height = 6)

# Elbow plot to determine optimal dimensions
elbow_plot <- ElbowPlot(seurat_obj, ndims = PARAMS$pca_dims)
save_plot(elbow_plot, "pca_elbow_plot.png", width = 8, height = 6)

# --------------------------------
# CLUSTERING ANALYSIS
# --------------------------------
print_section("CLUSTERING ANALYSIS")

# Find neighbors and clusters
cat("Finding neighbors and clusters...\n")
seurat_obj <- FindNeighbors(
  seurat_obj, 
  dims = PARAMS$neighbors_dims, 
  verbose = FALSE
)

seurat_obj <- FindClusters(
  seurat_obj, 
  resolution = PARAMS$cluster_resolution,
  verbose = FALSE
)

# Run UMAP
cat("Running UMAP...\n")
seurat_obj <- RunUMAP(
  seurat_obj, 
  dims = PARAMS$umap_dims, 
  verbose = FALSE
)

# Plot clusters
cluster_plot <- DimPlot(
  seurat_obj, 
  reduction = "umap", 
  label = TRUE, 
  repel = TRUE,
  pt.size = PARAMS$point_size,
  label.size = 5
) + 
  ggtitle(paste0("UMAP: Clusters (Resolution = ", PARAMS$cluster_resolution, ")")) +
  theme_minimal()

save_plot(cluster_plot, paste0("umap_clusters_res", PARAMS$cluster_resolution, ".png"))

# Print cluster summary
cluster_counts <- table(Idents(seurat_obj))
cat("Cluster Summary:\n")
print(cluster_counts)

# --------------------------------
# 🧬 MARKER GENE IDENTIFICATION
# --------------------------------
print_section("MARKER GENE IDENTIFICATION")

# Prepare object for marker finding
if ("RNA" %in% names(seurat_obj@assays)) {
  seurat_obj <- JoinLayers(seurat_obj)
}

# Find all markers
cat("Finding marker genes for all clusters...\n")
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = PARAMS$min_pct,
  logfc.threshold = PARAMS$logfc_threshold,
  verbose = FALSE
)

# Filter and select top markers
top_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = PARAMS$top_markers_n) %>%
  ungroup()

cat("Found", nrow(all_markers), "significant markers\n")
cat("Selected", nrow(top_markers), "top markers for visualization\n")

# Save marker tables
write.csv(all_markers, "tables/all_markers.csv", row.names = FALSE)
write.csv(top_markers, "tables/top_markers.csv", row.names = FALSE)
cat(" Marker tables saved to 'tables/' directory\n")

# --------------------------------
# FEATURE PLOT VISUALIZATION
# --------------------------------
print_section("FEATURE PLOT VISUALIZATION")

# Get genes to plot (ensure they exist in the object)
genes_to_plot <- unique(top_markers$gene[top_markers$gene %in% rownames(seurat_obj)])
cat("Plotting", length(genes_to_plot), "marker genes\n")

if (length(genes_to_plot) > 0) {
  # Limit number of genes to prevent oversized plots
  max_genes_per_plot <- 20
  
  if (length(genes_to_plot) > max_genes_per_plot) {
    cat("Too many genes (", length(genes_to_plot), "). Creating multiple plots with max", max_genes_per_plot, "genes each.\n")
    
    # Split genes into chunks
    gene_chunks <- split(genes_to_plot, ceiling(seq_along(genes_to_plot) / max_genes_per_plot))
    
    for (i in seq_along(gene_chunks)) {
      chunk_genes <- gene_chunks[[i]]
      
      # Create feature plots for this chunk
      feature_plots <- lapply(chunk_genes, function(gene) {
        create_feature_plot(seurat_obj, gene, " (Top Marker)")
      })
      
      # Calculate optimal layout
      n_genes <- length(chunk_genes)
      ncol <- min(4, n_genes)
      nrow <- ceiling(n_genes / ncol)
      
      # Combine plots with size limits
      combined_features <- wrap_plots(feature_plots, ncol = ncol)
      
      # Calculate safe dimensions (max 20 inches)
      plot_width <- min(20, ncol * 4)
      plot_height <- min(20, nrow * 3.5)
      
      save_plot(
        combined_features, 
        paste0("umap_top_markers_overlay_part", i, ".png"), 
        width = plot_width, 
        height = plot_height
      )
    }
  } else {
    # Create single plot for reasonable number of genes
    feature_plots <- lapply(genes_to_plot, function(gene) {
      create_feature_plot(seurat_obj, gene, " (Top Marker)")
    })
    
    # Calculate optimal layout
    n_genes <- length(genes_to_plot)
    ncol <- min(4, n_genes)
    nrow <- ceiling(n_genes / ncol)
    
    # Combine plots
    combined_features <- wrap_plots(feature_plots, ncol = ncol)
    
    # Calculate safe dimensions (max 20 inches)
    plot_width <- min(20, ncol * 4)
    plot_height <- min(20, nrow * 3.5)
    
    save_plot(
      combined_features, 
      "umap_top_markers_overlay.png", 
      width = plot_width, 
      height = plot_height
    )
  }
  
  # Create heatmap of top markers
  if (nrow(top_markers) > 0) {
    cat("Creating marker heatmap...\n")
    
    # Select fewer markers for cleaner heatmap (max 5 per cluster, max 50 total)
    max_markers_per_cluster <- 3
    max_total_markers <- 50
    
    heatmap_markers <- top_markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = max_markers_per_cluster) %>%
      ungroup() %>%
      slice_max(order_by = avg_log2FC, n = max_total_markers) %>%
      pull(gene)
    
    if (length(heatmap_markers) > 0) {
      # Calculate heatmap dimensions safely
      heatmap_height <- min(20, max(8, length(heatmap_markers) * 0.3 + 3))
      heatmap_width <- min(20, max(10, length(unique(Idents(seurat_obj))) * 0.8 + 4))
      
      heatmap_plot <- DoHeatmap(
        seurat_obj, 
        features = heatmap_markers,
        size = 3,
        angle = 90
      ) + 
        scale_fill_viridis_c() +
        theme(axis.text.y = element_text(size = 8))
      
      save_plot(heatmap_plot, "marker_heatmap.png", 
                width = heatmap_width, height = heatmap_height)
    } else {
      cat("No markers available for heatmap\n")
    }
  }
} else {
  cat("No valid marker genes found for plotting\n")
}

# --------------------------------
# 🏷️ CELL TYPE ANNOTATION PLOTS
# --------------------------------
print_section("CELL TYPE ANNOTATION VISUALIZATION")

# Check for existing annotations and create plots
annotation_plots <- list()

for (celltype_col in celltype_cols) {
  if (celltype_col %in% colnames(seurat_obj@meta.data)) {
    cat("🔄 Creating plot for:", celltype_col, "\n")
    
    # Count unique cell types
    n_types <- length(unique(seurat_obj@meta.data[[celltype_col]]))
    
    plot <- DimPlot(
      seurat_obj,
      group.by = celltype_col,
      label = TRUE, 
      repel = TRUE, 
      pt.size = PARAMS$point_size,
      label.size = ifelse(n_types > 20, 3, 4)
    ) + 
      ggtitle(paste0("UMAP: ", gsub("_", " ", stringr::str_to_title(celltype_col)))) +
      theme_minimal() +
      theme(legend.position = ifelse(n_types > 15, "none", "right"))
    
    annotation_plots[[celltype_col]] <- plot
    
    # Adjust plot size based on number of cell types
    plot_width <- ifelse(n_types > 15, 14, 10)
    plot_height <- ifelse(n_types > 15, 12, 8)
    
    save_plot(plot, paste0("umap_", celltype_col, ".png"), 
              width = plot_width, height = plot_height)
  }
}

# --------------------------------
# QUALITY CONTROL METRICS
# --------------------------------
print_section("QUALITY CONTROL VISUALIZATION")

# QC metrics plots
qc_metrics <- c("nFeature_RNA", "nCount_RNA")
if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
  qc_metrics <- c(qc_metrics, "percent.mt")
}

qc_plots <- lapply(qc_metrics, function(metric) {
  if (metric %in% colnames(seurat_obj@meta.data)) {
    FeaturePlot(
      seurat_obj, 
      features = metric, 
      reduction = "umap",
      pt.size = PARAMS$point_size
    ) +
      ggtitle(paste0("UMAP: ", metric)) +
      theme_minimal()
  }
})

# Remove NULL plots
qc_plots <- qc_plots[!sapply(qc_plots, is.null)]

if (length(qc_plots) > 0) {
  combined_qc <- wrap_plots(qc_plots, ncol = 2)
  save_plot(combined_qc, "umap_qc_metrics.png", width = 12, height = 8)
}

# --------------------------------
# SAVE PROCESSED OBJECT
# --------------------------------
print_section("SAVING RESULTS")

# Save processed Seurat object
output_file <- file.path("results", paste0("processed_", PARAMS$sample_name, ".rds"))
saveRDS(seurat_obj, output_file)
cat(" Processed Seurat object saved:", output_file, "\n")

# Create summary report
summary_stats <- data.frame(
  Metric = c("Total Cells", "Total Genes", "Variable Features", "PCA Dimensions", 
             "Clusters", "Significant Markers", "Top Markers"),
  Value = c(ncol(seurat_obj), nrow(seurat_obj), length(VariableFeatures(seurat_obj)),
            PARAMS$pca_dims, length(unique(Idents(seurat_obj))), 
            nrow(all_markers), nrow(top_markers))
)

write.csv(summary_stats, "tables/analysis_summary.csv", row.names = FALSE)

# --------------------------------
#  DISPLAY PLOTS
# --------------------------------
if (interactive()) {
  print_section("DISPLAYING PLOTS")
  
  # Display key plots
  print(var_plot_labeled)
  print(cluster_plot)
  
  if (length(annotation_plots) > 0) {
    lapply(annotation_plots, print)
  }
  
  if (exists("combined_features") && !is.null(combined_features)) {
    print(combined_features)
  }
}
