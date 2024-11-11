library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)



perform_differential_expression <- function(seurat_list, output_dir) {
  # Merge Seurat objects
  combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
  # Join layers
  combined_seurat <- JoinLayers(combined_seurat)
  
  # Add condition metadata if identifiable by project names
  combined_seurat$Condition <- ifelse(grepl("NML", combined_seurat@meta.data$orig.ident), "Normal", "Scleroderma")
  
  # Check normalization
  if (!"data" %in% slotNames(combined_seurat[["RNA"]])) {
    combined_seurat <- NormalizeData(combined_seurat)
  }
  
  # Run PCA for visualization
  combined_seurat <- RunPCA(combined_seurat)
  pca_plot <- DimPlot(combined_seurat, reduction = "pca", group.by = "Condition")
  ggsave(file.path(output_dir, "PCA_plot.png"), plot = pca_plot)
  
  # Perform differential expression analysis
  de_results <- FindMarkers(combined_seurat, ident.1 = "Normal", ident.2 = "Scleroderma", group.by = "Condition")
  write.csv(de_results, file.path(output_dir, "DE_results_Normal_vs_Scleroderma.csv"))
  
  # Create volcano plot
  de_results$gene <- rownames(de_results)
  de_results <- de_results %>%
    mutate(Significance = case_when(
      avg_log2FC > 1 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < -1 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  
  volcano_plot <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = "white")) +  
    labs(title = "Volcano Plot for DE Genes", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    geom_text(data = head(de_results[order(-abs(de_results$avg_log2FC)), ], 10), aes(label = gene), vjust = -0.5, size = 3)
  
  ggsave(file.path(output_dir, "volcano_plot.png"), plot = volcano_plot, width = 10, height = 7)
}


# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

seurat_files <- args[1:(length(args) - 1)]

# Print the Seurat files to verify the paths
cat("Seurat files :\n")
print(seurat_files)

# Extract the output directory from the arguments
output_dir <- args[length(args)]
cat("output dir :\n")
print(output_dir)

# Load each Seurat object and store in a list
seurat_list <- lapply(seurat_files, readRDS)

# Run the differential expression function
perform_differential_expression(seurat_list, output_dir)



