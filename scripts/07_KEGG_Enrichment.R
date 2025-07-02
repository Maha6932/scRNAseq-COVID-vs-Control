# ============================================================
# Script: 07_per_cluster_KEGG_enrichment.R
# Purpose: KEGG enrichment for up/downregulated genes per cluster
# ============================================================

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(readr)
library(forcats)

deg_dir <- "results/per_cluster_dge"
kegg_dir <- "results/per_cluster_kegg"

if (!dir.exists(kegg_dir)) dir.create(kegg_dir, recursive = TRUE)

# -------------------------------
# Get DEG files (excluding _top10)
# -------------------------------
deg_files <- list.files(deg_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)
deg_files <- deg_files[!str_detect(deg_files, "_top10\\.csv$")]

for (file_path in deg_files) {
  filename <- basename(file_path)
  celltype <- str_match(filename, "^DGE_(.*)\\.csv$")[,2]
  celltype_clean <- gsub("_", " ", celltype)
  
  message("\n🔹 Processing: ", celltype_clean)
  
  # ---------------------------------------
  # Load DEG CSV safely
  # ---------------------------------------
  deg_df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  
  if (!"gene" %in% colnames(deg_df)) {
    colnames(deg_df)[1] <- "gene"
  }
  if (!"gene" %in% colnames(deg_df)) {
    stop(" Gene column not found in ", filename)
  }
  
  # ---------------------------------------
  # Get up/downregulated gene lists
  # ---------------------------------------
  up_genes <- deg_df %>%
    filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>%
    pull(gene) %>% unique()
  
  down_genes <- deg_df %>%
    filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>%
    pull(gene) %>% unique()
  
  writeLines(up_genes, file.path(kegg_dir, paste0("up_genes_", celltype, ".txt")))
  writeLines(down_genes, file.path(kegg_dir, paste0("down_genes_", celltype, ".txt")))
  
  message("Up genes: ", length(up_genes), " | Down genes: ", length(down_genes))
  
  # ---------------------------------------
  # KEGG: Upregulated
  # ---------------------------------------
  if (length(up_genes) >= 5) {
    up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (!is.null(up_entrez) && nrow(up_entrez) > 0) {
      ekegg_up <- enrichKEGG(gene = up_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
      ekegg_up_df <- as.data.frame(ekegg_up)
      
      if (nrow(ekegg_up_df) > 0) {
        write.csv(ekegg_up_df, file.path(kegg_dir, paste0("KEGG_up_", celltype, ".csv")), row.names = FALSE)
        
        p_up <- barplot(ekegg_up, showCategory = 10) + 
          ggtitle(paste("KEGG Upregulated:", celltype_clean))
        ggsave(file.path(kegg_dir, paste0("KEGG_up_", celltype, ".png")), p_up, width = 8, height = 6)
      } else {
        message("No KEGG enrichment found for upregulated genes in ", celltype_clean)
      }
    }
  }
  
  # ---------------------------------------
  # KEGG: Downregulated
  # ---------------------------------------
  if (length(down_genes) >= 5) {
    down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (!is.null(down_entrez) && nrow(down_entrez) > 0) {
      ekegg_down <- enrichKEGG(gene = down_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
      ekegg_down_df <- as.data.frame(ekegg_down)
      
      if (nrow(ekegg_down_df) > 0) {
        write.csv(ekegg_down_df, file.path(kegg_dir, paste0("KEGG_down_", celltype, ".csv")), row.names = FALSE)
        
        p_down <- barplot(ekegg_down, showCategory = 10) + 
          ggtitle(paste("KEGG Downregulated:", celltype_clean))
        ggsave(file.path(kegg_dir, paste0("KEGG_down_", celltype, ".png")), p_down, width = 8, height = 6)
      } else {
        message("No KEGG enrichment found for downregulated genes in ", celltype_clean)
      }
    }
  }
}

message("\n Per-cluster KEGG enrichment completed for all cell types.\n")

# ============================================================
# Combine all KEGG CSV results into one table safely
# ============================================================

output_file <- file.path(kegg_dir, "combined_KEGG_results.csv")
kegg_files <- list.files(kegg_dir, pattern = "^KEGG_(up|down)_.*\\.csv$", full.names = TRUE)

combined_df <- data.frame()

for (file in kegg_files) {
  filename <- basename(file)
  direction <- ifelse(str_detect(filename, "^KEGG_up_"), "Up", "Down")
  celltype <- str_match(filename, "^KEGG_(up|down)_(.*)\\.csv$")[,3]
  celltype_clean <- gsub("_", " ", celltype)
  
  kegg_df <- read_csv(file, show_col_types = FALSE)
  
  # Skip empty files safely
  if (nrow(kegg_df) == 0) {
    message(" Skipping empty file: ", filename)
    next
  }
  
  # Force critical columns to character safely
  kegg_df <- kegg_df %>%
    mutate(
      geneID = as.character(geneID),
      ID = as.character(ID),
      Description = as.character(Description)
    ) %>%
    mutate(CellType = celltype_clean, Direction = direction)
  
  combined_df <- bind_rows(combined_df, kegg_df)
  
  message(" Added: ", filename, " | Pathways: ", nrow(kegg_df))
}

View(combined_df)
write_csv(combined_df, output_file)

message("\n combined KEGG results saved to: ", output_file)

# Read combined KEGG file

kegg_df <- combined_df
# Check columns
print(colnames(kegg_df))

View(kegg_df)

# Select top pathways across all cell types for clean plotting
top_kegg <- kegg_df %>%
  group_by(Description, Direction) %>%
  summarise(Total_Count = sum(Count), min_p = min(p.adjust), .groups = "drop") %>%
  arrange(min_p) %>%
  slice_head(n = 20)

# Filter original df for these top pathways only
plot_df <- kegg_df %>%
  filter(Description %in% top_kegg$Description)
plot_df <- plot_df %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text = x))))
sum(is.na(plot_df$GeneRatio))

# Dot plot
p <- ggplot(plot_df, aes(
  x = fct_reorder(CellType, GeneRatio, .desc = TRUE),
  y = fct_reorder(Description, GeneRatio),
  size = GeneRatio,
  color = p.adjust
)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(
    low = "red",
    high = "blue",
    trans = "reverse",
    name = "Adj p-value"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "GeneRatio"
  ) +
  facet_wrap(~Direction, scales = "free_y") +
  labs(
    x = "Cell Type",
    y = "KEGG Pathway",
    title = "KEGG Pathway Enrichment Across Cell Types",
    subtitle = "Dot size = GeneRatio, color = adjusted p-value"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p)

# Save
ggsave("results/per_cluster_kegg/kegg_dotplot.png", p, width = 12, height = 8, dpi = 300)
message("KEGG dot plot saved: results/per_cluster_kegg/kegg_dotplot.png")
