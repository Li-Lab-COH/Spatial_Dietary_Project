# Author: Jonathan Anzules
# Project: Spatial Dietary
# Objective: Spatial visualization of these genes: "Ulk1", "Atg2a", "Lc3", "Becn1"

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Define directories
input_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/data/Rogelio/Analysis_031625/BANKSY_Normalized_QC_Filtered"
out_folder <- "/Users/janzules/Roselab/Spatial/dietary_project/figures/Kenichi_figures/Ar_gene"

# Define autophagy genes
# autophagy_genes <- c("Ulk1", "Atg2a", "Map1lc3a", "Becn1")
# Define the gene of interest
ar_gene <- "Ar"

# List of files
file_list <- list(
  RT_Files = c("F07833_5_RT.rds", "F07835_3_RT.rds", "F07837_31_RT.rds", 
               "F07836_11_RT.rds", "F07838_32_RT.rds"),
  
  CRRT_Files = c("F07833_8_CRRT.rds", "F07834_8_CRRT.rds", 
                 "F07835_7_CRRT.rds", "F07836_10_CRRT.rds"),
  
  SKO_CRRT_Files = c("F07838_36_CR_RT.rds", "F07838_38_CR_RT.rds", 
                     "F07837_37_CR_RT.rds")
)

# Loop through each category in file_list
for (category in names(file_list)) {
  
  # Define subfolder path
  # category_folder <- file.path(out_folder, category)
  
  # Ar folder
  category_folder <- file.path(out_folder, category, "Ar_Gene")
  
  # Create folder if it does not exist
  if (!dir.exists(category_folder)) {
    dir.create(category_folder, recursive = TRUE)
    message("Created directory: ", category_folder)
  }
  
  # Loop through each file in the category
  for (file in file_list[[category]]) {
    
    # Remove ".rds" extension
    file_name <- gsub("\\.rds$", "", file)
    
    # Remove the F##### pattern (assuming it always starts with 'F' followed by 5 digits)
    file_name <- gsub("^F\\d+_", "", file_name)
    
    # Print extracted name for debugging
    message("Processing file: ", file, " -> Cleaned name: ", file_name)
    
    # Load the Seurat object
    seurat_obj <- readRDS(file.path(input_dir, file))
    
    # Set the default assay to Spatial.008um
    DefaultAssay(seurat_obj) <- "Spatial.008um"
    
    #-------------Ar gene----------------------------
    # Generate spatial expression plot for Ar gene
    ar_gene_plot <- SpatialFeaturePlot(seurat_obj, features = ar_gene,
                                       image.alpha = 0.2,
                                       pt.size.factor = 10,
                                       crop = TRUE) +
      ggtitle(paste("Ar Gene Expression for", file_name)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    # Define output file path
    ar_plot_path <- file.path(ar_folder, paste0(file_name, "_Ar_expression.tiff"))
    
    # Save plot as high-resolution TIFF file
    ggsave(ar_plot_path, plot = ar_gene_plot, dpi = 300, width = 8, height = 8, units = "in", device = "tiff")
    
    message("Saved Ar gene plot for: ", file_name)
    
    #--------------------Autophagy Genes----------------
    # # Compute autophagy score
    # genes_present <- autophagy_genes[autophagy_genes %in% rownames(seurat_obj[["Spatial.008um"]])]
    # if (length(genes_present) > 0) {
    #   seurat_obj$autophagy_score <- colMeans(as.matrix(GetAssayData(seurat_obj, assay = "Spatial.008um", slot = "data")[genes_present, , drop = FALSE]))
    # } else {
    #   message("Skipping ", file_name, " - No autophagy genes found.")
    #   next  # Skip if no genes are found
    # }
    # 
    # # Generate spatial expression plot for autophagy genes
    # autophagy_gene_plot <- SpatialFeaturePlot(seurat_obj, features = autophagy_genes,
    #                                           image.alpha = 0.2,
    #                                           pt.size.factor = 10,
    #                                           crop = TRUE) +
    #   plot_annotation(title = paste("Autophagy Genes for", file_name))
    # 
    # # Generate the combined autophagy score plot
    # autophagy_score_plot <- SpatialFeaturePlot(seurat_obj, features = "autophagy_score",
    #                                            image.alpha = 0.2,
    #                                            pt.size.factor = 10,
    #                                            crop = TRUE) +
    #   ggtitle(paste("Autophagy Score for", file_name, "\n(Mean Expression of Ulk1, Atg2a, Map1lc3a, Becn1)")) +
    #   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    # 
    # # Define output file paths
    # gene_plot_path <- file.path(category_folder, paste0(file_name, "_autophagy_genes.tiff"))
    # score_plot_path <- file.path(category_folder, paste0(file_name, "_autophagy_score.tiff"))
    # 
    # # Save plots as high-resolution TIFF files
    # ggsave(gene_plot_path, plot = autophagy_gene_plot, dpi = 300, width = 8, height = 8, units = "in", device = "tiff")
    # ggsave(score_plot_path, plot = autophagy_score_plot, dpi = 300, width = 8, height = 8, units = "in", device = "tiff")
    
    message("Saved plots for: ", file_name)
  }
}































library(Seurat)
library(ggplot2)
library(patchwork)  # Needed for adding a title

input_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/data/Rogelio/Analysis_031625/BANKSY_Normalized_QC_Filtered"
out_folder <- "/Users/janzules/Roselab/Spatial/dietary_project/figures/Kenichi_figures"



file_list <- list(
  RT_Files = c("F07833_5_RT.rds", "F07835_3_RT.rds", "F07837_31_RT.rds", 
               "F07836_11_RT.rds", "F07838_32_RT.rds"),
  
  CRRT_Files = c("F07833_8_CRRT.rds", "F07834_8_CRRT.rds", 
                 "F07835_7_CRRT.rds", "F07836_10_CRRT.rds"),
  
  SKO_CRRT_Files = c("F07838_36_CR_RT.rds", "F07838_38_CR_RT.rds", 
                  "F07837_37_CR_RT.rds")
)

autophagy_genes <- c("Ulk1", "Atg2a", "Map1lc3a", "Becn1")
# Loop through the file list
for (category in names(file_list)) {
  
  # Define subfolder path
  category_folder <- file.path(out_folder, category)
  
  # Create folder if it does not exist
  if (!dir.exists(category_folder)) {
    dir.create(category_folder, recursive = TRUE)
    message("Created directory: ", category_folder)
  }
  
  # Loop through files in the category
  for (file in file_list[[category]]) {
    
    # Remove ".rds" extension
    file_name <- gsub("\\.rds$", "", file)
    
    # Remove the F##### pattern
    file_name <- gsub("^F\\d+_", "", file_name)
    
    # Print extracted name (for debugging)
    message("Processed file: ", file, " -> Cleaned name: ", file_name)
    
  }
}








#Lc3 = Map1lc3a


seurat_obj <- readRDS(paste0(input_dir,'/', file_list$RT_Files[1]))

DefaultAssay(seurat_obj) <- "Spatial.008um"

seurat_obj$autophagy_score <- colMeans(as.matrix(GetAssayData(seurat_obj, assay = "Spatial.008um", slot = "data")[autophagy_genes, , drop = FALSE]))

SpatialFeaturePlot(seurat_obj, features = autophagy_genes,
                   image.alpha = 0.2,
                   pt.size.factor = 10,
                   crop = TRUE)+
  plot_annotation(title = "Autophagy Genes for [obj name]")

# Generate the plot and store it in a variable
autophagy_plot <- SpatialFeaturePlot(seurat_obj, features = "autophagy_score",
                                     image.alpha = 0.2,
                                     pt.size.factor = 10,
                                     crop = TRUE) +
  ggtitle("Autophagy Score for [obj name]\n(Mean Expression of Ulk1, Atg2a, Map1lc3a, Becn1)") +  # Add title
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))  # Center and style title

# Display the plot
print(autophagy_plot)
