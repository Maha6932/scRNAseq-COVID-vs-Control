# ============================
# Step 5: Identify Marker Genes for Each Cluster
# ============================

# Load required libraries
library(Seurat)
library(dplyr)

# Load the clustered Seurat object
seurat_obj <- readRDS("results/processed_C51_L01.rds")

# Check if clusters exist
table(Idents(seurat_obj))

# Set RNA as the default assay (important for marker detection)
DefaultAssay(seurat_obj) <- "RNA"


# View part of the RNA assay metadata (sanity check)
# Note: This will throw an error — meta.data belongs to seurat_obj@meta.data, not assay
# So this line is unnecessary or should be corrected to:
head(seurat_obj@meta.data)


# Join data layers (useful if working with multi-layered Seurat v5 object)
seurat_obj <- JoinLayers(seurat_obj)

# Try testing marker genes between two clusters (cluster 1 vs 0)
markers_test <- FindMarkers(seurat_obj, ident.1 = 1, ident.2 = 0,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
head(markers_test)

# Find all marker genes for all clusters (Wilcoxon Rank Sum test by default)
markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE,         # Keep only upregulated markers
                          min.pct = 0.25,          # At least 25% of cells in cluster express the gene
                          logfc.threshold = 0.25)  # Minimum log fold-change

# Save the full marker gene list
write.csv(markers, file = "results/C51_L01_cluster_markers.csv", row.names = FALSE)

# Extract and view top 10 markers per cluster
top10 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

# View results in a table viewer
View(top10)
