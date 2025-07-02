# ==============================================================================
# SINGLE CELL RNA-SEQ ANALYSIS: COVID-19 vs CONTROL LUNG SAMPLES
# Script 3: Normalization, Clustering, and Marker Gene Identification
# ==============================================================================
# Purpose:
# 1. Normalize gene expression data
# 2. Find highly variable genes
# 3. Perform dimensionality reduction (PCA, UMAP)
# 4. Cluster cells and identify marker genes
# ==============================================================================

# Load Required Libraries ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

print("Starting normalization, clustering, and marker identification...")

# SECTION 1: LOAD FILTERED DATA ----
print("Loading filtered Seurat object...")
seurat_obj <- readRDS("results/filtered_C51_L01.rds")

print(paste("Loaded object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes"))
print("Condition distribution:")
print(table(seurat_obj$condition))

# SECTION 2: DATA NORMALIZATION ----
print("Normalizing gene expression data...")

# Log-normalization: log(counts per 10,000 + 1)
# This accounts for differences in sequencing depth between cells
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

cat("✅ Data normalization complete\n")

# SECTION 3: FIND HIGHLY VARIABLE FEATURES ----
print("Identifying highly variable genes...")

# Find genes that vary the most across cells
# These are the most informative for finding cell types
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)

# Get the top 10 most variable genes
top10_variable <- head(VariableFeatures(seurat_obj), 10)
print("Top 10 most variable genes:")
print(top10_variable)

# Plot variable features
variable_plot <- VariableFeaturePlot(seurat_obj)
top10_plot <- LabelPoints(plot = variable_plot, 
                          points = top10_variable, 
                          repel = TRUE,
                          xnudge = 0,
                          ynudge = 0)
top10_plot
ggsave("figures/highly_variable_genes.png", top10_plot, 
       width = 10, height = 6, dpi = 300)

cat("✅ Found", length(VariableFeatures(seurat_obj)), "highly variable genes\n")

# SECTION 4: DATA SCALING ----
print("Scaling data for PCA...")

# Scale data: zero mean, unit variance for each gene
# This prevents highly expressed genes from dominating PCA
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

cat("✅ Data scaling complete\n")

# SECTION 5: PRINCIPAL COMPONENT ANALYSIS (PCA) ----
print("Running PCA...")

# PCA on highly variable genes
seurat_obj <- RunPCA(seurat_obj, 
                     features = VariableFeatures(seurat_obj),
                     verbose = FALSE)

# Visualize PCA results
pca_plot <- DimPlot(seurat_obj, reduction = "pca", group.by = "condition") +
  ggtitle("PCA: COVID vs Control")

# Show top genes contributing to first 2 PCs
print("Top genes contributing to PC1 and PC2:")
print(seurat_obj[["pca"]], dims = 1:2, nfeatures = 5)

# Elbow plot to determine optimal number of PCs
elbow_plot <- ElbowPlot(seurat_obj, ndims = 30) +
  ggtitle("PCA Elbow Plot - Determining Optimal Dimensions")

# Save PCA plots
ggsave("figures/pca_by_condition.png", pca_plot, width = 8, height = 6, dpi = 300)
ggsave("figures/elbow_plot.png", elbow_plot, width = 8, height = 6, dpi = 300)

cat("✅ PCA complete\n")

# SECTION 6: DETERMINE OPTIMAL DIMENSIONS ----
# Based on elbow plot, choose number of PCs to use
# Look for the "elbow" where variance explained levels off
n_dims <- 15  # Adjust based on your elbow plot

cat("Using", n_dims, "dimensions for clustering and UMAP\n")

# SECTION 7: CELL CLUSTERING ----
print("Clustering cells...")

# Build k-nearest neighbor graph
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE)

# Cluster cells using Louvain algorithm
# Resolution parameter controls cluster granularity (higher = more clusters)
resolution <- 0.5
seurat_obj <- FindClusters(seurat_obj, 
                           resolution = resolution, 
                           verbose = FALSE)

n_clusters <- length(levels(seurat_obj$seurat_clusters))
cat("✅ Found", n_clusters, "clusters at resolution", resolution, "\n")

# SECTION 8: UMAP VISUALIZATION ----
print("Running UMAP...")

# UMAP for 2D visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims, verbose = FALSE)

# Create UMAP plots
p1 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              label = TRUE, 
              repel = TRUE,
              label.size = 4) +
  ggtitle("UMAP: Cell Clusters") +
  theme(legend.position = "right")

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "condition",
              cols = c("Control" = "blue", "COVID" = "red")) +
  ggtitle("UMAP: COVID vs Control") +
  theme(legend.position = "right")

# Combine plots
combined_umap <- p1 + p2
print("Displaying UMAP plots...")
print(combined_umap)

# Save individual and combined plots
ggsave("figures/umap_clusters.png", p1, width = 8, height = 6, dpi = 300)
ggsave("figures/umap_covid_vs_control.png", p2, width = 8, height = 6, dpi = 300)
ggsave("figures/umap_combined.png", combined_umap, width = 14, height = 6, dpi = 300)

cat("✅ UMAP visualization complete\n")

# SECTION 9: CLUSTER COMPOSITION ANALYSIS ----
print("Analyzing cluster composition...")

# Create table showing condition distribution per cluster  
cluster_composition <- table(seurat_obj$seurat_clusters, seurat_obj$condition)
print("Cluster composition (COVID vs Control):")
print(cluster_composition)

# Convert to percentages
cluster_composition_pct <- prop.table(cluster_composition, margin = 1) * 100
print("Cluster composition (percentages):")
print(round(cluster_composition_pct, 1))

# Save composition data
write.csv(cluster_composition, "results/cluster_composition_counts.csv")
write.csv(cluster_composition_pct, "results/cluster_composition_percentages.csv")

# SECTION 10: MARKER GENE IDENTIFICATION ----
print("Finding marker genes for each cluster...")

# Join layers for compatibility with newer Seurat versions
seurat_obj <- JoinLayers(seurat_obj)

# Find marker genes that distinguish each cluster
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,        # Only positive markers
  min.pct = 0.25,         # Gene must be detected in at least 25% of cells
  logfc.threshold = 0.25, # Minimum log fold-change
  verbose = FALSE
)

print(paste("Found", nrow(markers), "significant marker genes"))

# Show top markers per cluster
top_markers_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_head(n = 3)

print("Top 3 marker genes per cluster:")
print(top_markers_per_cluster[, c("cluster", "gene", "avg_log2FC", "p_val_adj")])

# Save all markers
write.csv(markers, "results/cluster_markers_all.csv", row.names = FALSE)

cat("✅ Marker gene identification complete\n")

# SECTION 11: VISUALIZE TOP MARKER GENES ----
print("Creating heatmap of top marker genes...")

# Get top 5 markers per cluster for heatmap
top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Scale all genes for heatmap
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Create heatmap
p_heatmap <- DoHeatmap(
  seurat_obj,
  features = top5_markers$gene,
  group.by = "seurat_clusters",
  size = 3
) + 
  ggtitle("Top 5 Marker Genes per Cluster") +
  theme(axis.text.y = element_text(size = 6))

ggsave("figures/heatmap_top5_markers.png", p_heatmap, 
       width = 12, height = 10, dpi = 300)

cat("✅ Heatmap visualization complete\n")

# SECTION 12: FEATURE PLOTS FOR KEY MARKERS ----
print("Creating feature plots for selected marker genes...")

# Select some interesting markers to visualize
if(nrow(top_markers_per_cluster) > 0) {
  selected_genes <- top_markers_per_cluster$gene[1:min(8, nrow(top_markers_per_cluster))]
  
  feature_plots <- FeaturePlot(seurat_obj, 
                               features = selected_genes, 
                               ncol = 4,
                               reduction = "umap")
  
  ggsave("figures/feature_plots_top_markers.png", feature_plots,
         width = 16, height = 8, dpi = 300)
}

# SECTION 13: SAVE RESULTS ----
print("Saving clustered object and summary...")

# Save the final clustered object
saveRDS(seurat_obj, "results/clustered_seurat_object.rds")

# Create analysis summary
analysis_summary <- data.frame(
  Parameter = c("Total_cells", "Total_genes", "Variable_features", "PCs_used", 
                "Clustering_resolution", "Number_of_clusters", "Total_markers"),
  Value = c(ncol(seurat_obj), nrow(seurat_obj), length(VariableFeatures(seurat_obj)),
            n_dims, resolution, n_clusters, nrow(markers))
)

write.csv(analysis_summary, "results/analysis_parameters.csv", row.names = FALSE)

# Clean up memory
gc()

cat("\n Script 3 completed successfully!\n")
cat("   Clustered object saved to: results/clustered_seurat_object.rds\n")
cat("   Marker genes saved to: results/cluster_markers_all.csv\n")
cat("   Visualizations saved to: figures/\n")
cat("   Ready for cell type annotation and differential expression!\n")