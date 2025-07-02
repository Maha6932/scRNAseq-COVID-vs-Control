# ==========================================================
# Script: 06_per_cluster_DGE_clean.R
# Purpose: Per-cluster COVID vs Control DEGs using Seurat v5
# ==========================================================

library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)

# -------------------------------
# Load clustered Seurat object
# -------------------------------
cat("\n🔹 Loading clustered Seurat object...\n")
seurat_obj <- readRDS("results/processed_C51_L01.rds")

# Join layers to enable FindMarkers in Seurat v5
cat("\n🔹 Joining layers in Seurat object for DE analysis...\n")
seurat_obj <- JoinLayers(seurat_obj)

# -------------------------------
# Check metadata
# -------------------------------
cat("\n Condition breakdown:\n")
print(table(seurat_obj$condition, useNA = "ifany"))

cat("\n Cell type breakdown:\n")
print(table(seurat_obj$celltype_intermediate, seurat_obj$condition, useNA = "ifany"))

cat("\n Cluster breakdown:\n")
print(table(seurat_obj$seurat_clusters))

# -------------------------------
# Create output directories
# -------------------------------
dir.create("results/per_cluster_dge", recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Per-cluster DEGs
# -------------------------------
cell_types <- sort(unique(na.omit(seurat_obj$celltype_intermediate)))
cat("\n Cell types to process:\n")
print(cell_types)

DefaultAssay(seurat_obj) <- "RNA"

for (ct in cell_types) {
  message("\n Processing: ", ct)
  
  # Subset Seurat object for this cell type
  subset_obj <- subset(seurat_obj, subset = celltype_intermediate == ct)
  
  if (ncol(subset_obj) == 0) {
    message(" Skipping ", ct, ": no cells found after subsetting.")
    next
  }
  
  present_conditions <- unique(subset_obj$condition)
  if (!all(c("COVID", "Control") %in% present_conditions)) {
    message(" Skipping ", ct, ": missing one of the conditions (COVID or Control).")
    next
  }
  
  Idents(subset_obj) <- "condition"
  
  dge <- tryCatch({
    FindMarkers(
      object = subset_obj,
      ident.1 = "COVID",
      ident.2 = "Control",
      min.pct = 0.2,
      logfc.threshold = 0.25
    )
  }, error = function(e) {
    message(" Error in FindMarkers for ", ct, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(dge) && nrow(dge) > 0) {
    # Save full DEGs
    out_path <- paste0("results/per_cluster_dge/DGE_", gsub("[ /]", "_", ct), ".csv")
    write.csv(dge, out_path)
    message(" Saved DEGs for ", ct, " to ", out_path)
    
    # Save top 10 DEGs (by avg_log2FC) for quick interpretation
    top10 <- dge %>%
      dplyr::arrange(desc(avg_log2FC)) %>%
      dplyr::slice_head(n = 10)
    
    top10_path <- paste0("results/per_cluster_dge/DGE_", gsub("[ /]", "_", ct), "_top10.csv")
    write.csv(top10, top10_path)
    message(" Saved top 10 DEGs for ", ct, " to ", top10_path)
    
  } else {
    message(" No DEGs found for ", ct)
  }
}

cat("\n DEG analysis complete!")

# -------------------------------
# Save session info for reproducibility
# -------------------------------
writeLines(capture.output(sessionInfo()), "results/per_cluster_dge/sessionInfo.txt")

# Directory with your DGE CSVs
dge_dir <- "results/per_cluster_dge"
plot_dir <- "figures/per_cluster_volcano"

if(!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}

dge_files <- list.files(dge_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

for(f in dge_files){
  ct <- gsub("DGE_|\\.csv", "", basename(f))
  message("\n Generating volcano plot for: ", ct)
  
  dge <- read.csv(f, row.names = 1)
  
  #skip if empty
  if(nrow(dge) == 0){
    message("Skipping ", ct, ": No DEGs")
  }
  
  p <- EnhancedVolcano(
    dge, 
    lab = rownames(dge),
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    title = paste("Volcano:", ct),
    subtitle = "COVID vs Control",
    caption = " FC cutoff = 0.5; adj p < 0.05",
    labSize = 3,
    pointSize = 2,
    col = c("grey30", "forestgreen", "royalblue", "red2")
  )
  
  out_path <- file.path(plot_dir, paste0("Volcano_", ct, ".png"))
  ggsave(out_path, p, width = 7, height = 6, dpi = 300)
  
  message("Volcano plot saved: ", out_path)
}
