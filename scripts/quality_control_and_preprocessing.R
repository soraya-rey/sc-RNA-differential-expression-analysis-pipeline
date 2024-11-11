library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

# Step 1 : quality control of the reads and preprocessing

perform_quality_control <- function(seurat_object, output_dir, project_name) {
  # Calculate mitochondrial content
  cat("Calculating mitochondrial content...\n")
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate QC plots
  cat("Generating QC plots...\n")
  vln_plot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  ggsave(filename = file.path(output_dir, paste0(project_name, "_QC_violin.png")), plot = vln_plot)
  
  scatter_plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  ggsave(filename = file.path(output_dir, paste0(project_name, "_scatter_nCount_vs_mt.png")), plot = scatter_plot1)
  
  scatter_plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(filename = file.path(output_dir, paste0(project_name, "_scatter_nCount_vs_nFeature.png")), plot = scatter_plot2)
  
  # Filter low-quality cells
  cat("Filtering low-quality cells...\n")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  
  return(seurat_object)
}


perform_preprocessing <- function(seurat_object, output_dir, sample_name) {
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  saveRDS(seurat_object, file.path(output_dir, paste0(sample_name, "_seurat_object.rds")))
  return(seurat_object)
}

load_and_process_data <- function(data_dirs, output_dir) {
  seurat_list <- list()
  for (dir_name in names(data_dirs)) {
    seurat_obj <- CreateSeuratObject(counts = Read10X(data.dir = data_dirs[[dir_name]]), project = dir_name, min.cells = 3, min.features = 200)
    seurat_obj <- perform_quality_control(seurat_obj, output_dir, dir_name)
    seurat_obj <- perform_preprocessing(seurat_obj, output_dir, dir_name)
    seurat_list[[dir_name]] <- seurat_obj
  }
  return(seurat_list)
}


# Main execution
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]    
output_dir <- args[2]

cat("Data directory:", data_dir, "\n")


# Get all sample directories
sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
sample_dict <- setNames(sample_dirs, basename(sample_dirs))

# Print sample_dict to check its contents
cat("Sample directories and their names:\n")
print(sample_dict)

# Process the data for all samples
seurat_objects <- load_and_process_data(sample_dict, output_dir)