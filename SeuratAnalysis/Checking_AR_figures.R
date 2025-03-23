# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)

# Define directories
input_dir <- "C:/Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/Rogelio/Analysis/BANKSY_Normalized/Seurat_Objects/"
output_dir <- "C:/Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/Figures/Visualizing_sections/Ar_Gene"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created directory: ", output_dir)
}

# Define the gene of interest
ar_gene <- "Ar"

# Get list of .rds files from the input directory
file_list <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Process each file
for (file_path in file_list) {
  
  # Extract file name without path
  file_name <- basename(file_path)
  
  # Remove ".rds" extension
  sample_name <- gsub("\\.rds$", "", file_name)
  
  # Print extracted name for debugging
  message("Processing file: ", file_name, " -> Sample name: ", sample_name)
  
  # Load the Seurat object
  seurat_obj <- readRDS(file_path)
  
  # Set the default assay to Spatial.008um
  DefaultAssay(seurat_obj) <- "Spatial.008um"
  
  # Generate spatial expression plot for Ar gene
  ar_gene_plot <- SpatialFeaturePlot(seurat_obj, features = ar_gene,
                                     image.alpha = 0.2,
                                     pt.size.factor = 10,
                                     crop = TRUE) +
    ggtitle(paste("Ar Gene Expression for", sample_name)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  # Define output file path
  ar_plot_path <- file.path(output_dir, paste0(sample_name, "_Ar_expression.tiff"))
  
  # Save plot as high-resolution TIFF file
  ggsave(ar_plot_path, plot = ar_gene_plot, dpi = 300, width = 8, height = 8, units = "in", device = "tiff")
  
  message("Saved Ar gene plot for: ", sample_name)
}

message("Processing complete. All plots saved in: ", output_dir)





#--------------------------------- Solo view ----------------------------------
fig_directory_Ar <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/"

DefaultAssay(prostate_ST) <- "Spatial.008um"

ar_plot <- SpatialFeaturePlot(prostate_ST, features = "Ar",interactive = FALSE,
                   image.alpha = 0.2,
                   pt.size.factor = 3
                   )


ggsave(
  filename = file.path(fig_directory_Ar, "F07834_28_LFRT_AR.png"),
  plot = ar_plot,
  width = 8,
  height = 6,
  dpi = 300
)
# SpatialFeaturePlot(seurat_obj, features = ar_gene,
#                    image.alpha = 0.2,
#                    pt.size.factor = 10,
#                    crop = TRUE) 


SpatialFeaturePlot(prostate_ST, features = "Adgre1",interactive = FALSE,
                   image.alpha = 0.2,
                   pt.size.factor = 6
)
