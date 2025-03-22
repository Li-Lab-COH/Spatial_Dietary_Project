# Load required libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(pryr)
library(future)

# Prevent unwanted parallelization
plan("sequential")

# # Define paths
input_dir <- "~/1Work/RoseLab/Spatial/Dietary_Project/data/Rogelio/Analysis_final/BANKSY_Normalized_QC_Filtered_minUMI_50_minGene50_MT5/"
output_base <- "~/1Work/RoseLab/Spatial/Dietary_Project/Figures/Rogelio_final"
# input_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/data/Rogelio/Analysis_031625/BANKSY_Normalized_QC_Filtered"
# output_base <- "/Users/janzules/Roselab/Spatial/dietary_project/figures/from_Rogelio_objects/cell_type_search"

# List .rds files (ignoring directories such as "Plots" and "Cluster_Markers")
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Define grouping names (resolution folders)
groupings <- c("BANKSY_snn_res.0.8", "BANKSY_snn_res.0.5")

# Define marker sets for each section as a nested list
marker_sets <- list(
  Immune_Cell_Markers = list(
    NK_Cell = c("Gzmb", "Ifng", "Prf1"),
    T_cells = c("Cd3e", "Cd4", "Cd8a", "Foxp3", "Tbx21", "Gzmk", "Il2ra", "Gzmb", "Prf1"),
    TCR_Markers = c("Trac", "Trbc1", "Trbc2", "Trdc"),
    B_cells = c("Cd19", "Cd79a", "Ms4a1", "Ighm", "Pax5", "Cd22"),
    Macrophage = c("Adgre1", "Cd68", "Itgam", "Mrc1", "Tnf", "Il1b"),
    Dendritic_Cell = c("Itgax", "H2-Ab1", "Cd74", "Xcr1", "Cd86", "Batf3")
  ),
  Other_Immune_Cells = list(
    MDSCs = c("Ly6g", "Arg1", "Nos2", "S100a8", "S100a9"),
    Regulatory_T_Cells = c("Foxp3", "Il2ra", "Ctla4", "Tigit", "Ikzf2"),
    Exhausted_T_Cells = c("Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tox"),
    TAM_anti = c("Nos2", "Cd86", "Il12b"),
    TAM_pro = c("Mrc1", "Arg1"),
    N1_Neutrophils = c("Cxcr3", "Il12b", "Nos2", "Tnf", "Il1b", "Ccl3", "Ccl4", "Ccl5", "Mpo"),
    N2_Neutrophils = c("S100a8", "S100a9", "Arg1", "Tgfb1", "Vegfa"),
    Mast_Cells = c("Kit", "Fcer1a", "Mcpt1", "Cpa3"),
    Fibroblast_Associated_Immune_Cells = c("Fap", "Pdgfrb", "Col1a1", "Acta2")
  ),
  MyC_CaP_and_SKO = list(
    Mycs = c("Myc", "Mycn", "Mycl", "Ar"),
    Cell_growth = c("Ccna2", "Ccnd1", "Ccne1", "Cdk1", "Cdc25a", "Rcc1", "Odc1", "Tert", "Tk1"),
    Survival = c("Cdkn2a", "Trp53", "Bax"),
    Apoptosis = c("Bcl2", "Fas"),
    Metabolism = c("Ldha", "Cad", "Dhfr", "Eif2s1", "Eif4e", "Rida"),
    Downregulated_differentiation = c("Cebpa", "Thbs1", "H2-K1", "H2-D1", "Col1a1", "Col1a2", "Col6a3", "Col3a1"),
    Prostate_specific_markers = c("Epcam", "Nkx3-1", "Krt8", "Krt18"),
    Proliferation_markers = c("Mki67", "E2f1", "E2f2", "E2f3"),
    Myc_metabolic_transporters = c("Slc7a5"),
    Probasin_expression = c("Pbsn")
  ),
  Epithelial_Cell_Markers = list(
    Luminal = c("Tgm4", "Msmb", "Ppp1r1b"),
    LumP = c("Ppp1r1b", "Cldn10"),
    Periurethral = c("Ly6d", "Aqp3", "Ppp1r1b")
  )
)

# Define the new, crisper plotting function
generate_violin_plots <- function(marker_list, cell_type, pointSize = 0.3) {
  plots <- lapply(marker_list, function(gene) {
    # Try to fetch data safely
    expr_data <- tryCatch(
      FetchData(seurat_obj, vars = gene),
      error = function(e) {
        message(paste("Warning: Gene", gene, "not found in Seurat object. Skipping..."))
        return(NULL)  # Skip plotting if gene isn't found
      }
    )
    
    if (is.null(expr_data)) {
      return(ggplot() + ggtitle(paste(cell_type, "-", gene, "(Not Found)")))
    }
    
    meta_data <- seurat_obj@meta.data
    plot_data <- data.frame(Expression = expr_data[, 1], Group = meta_data[[grouping_name]])
    
    # Remove zero values
    plot_data <- plot_data[plot_data$Expression > 0, ]
    
    # If all values are zero, return a blank plot
    if (nrow(plot_data) == 0) {
      return(ggplot() + ggtitle(paste(cell_type, "-", gene, "(No Non-Zero Values)")))
    }
    
    # Suppress warnings only for the ggplot call
    suppressWarnings(
      suppressMessages(
        ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
          geom_violin(scale = "width") +
          geom_jitter(width = 0.2, size = pointSize, alpha = 0.5) +
          coord_cartesian(ylim = c(min(plot_data$Expression), NA)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle(paste(cell_type, "-", gene))
      )
    )
  })
  
  wrap_plots(plots)  # Arrange plots together
}

# Loop over each .rds file (each Seurat object)
for(rds_file in rds_files) {
  object_name <- tools::file_path_sans_ext(basename(rds_file))
  cat("\nProcessing object:", object_name, "\n")
  
  # Load the Seurat object
  seurat_obj <- readRDS(rds_file)
  cat("Memory after loading:", pryr::mem_used(), "\n")
  
  # Loop over each grouping resolution
  for(grouping in groupings) {
    grouping_name <- grouping  # set the grouping for use in the plotting function
    grouping_folder <- file.path(output_base, object_name, grouping_name)
    if (!dir.exists(grouping_folder)) dir.create(grouping_folder, recursive = TRUE)
    
    # Loop over each section (e.g., Immune_Cell_Markers, Other_Immune_Cells, etc.)
    for(section in names(marker_sets)) {
      section_folder <- file.path(grouping_folder, section)
      if (!dir.exists(section_folder)) dir.create(section_folder)
      
      # Loop over each marker set / cell type within the section
      for(cell_type in names(marker_sets[[section]])) {
        markers <- marker_sets[[section]][[cell_type]]
        
        # Use the new plotting function to generate the combined plot
        combined_plot <- generate_violin_plots(markers, cell_type, pointSize = 0.3)
        
        # Construct the output file name (e.g., F07834_9_KDRT_Mycs.png)
        output_file <- file.path(section_folder, paste0(object_name, "_", cell_type, ".png"))
        
        # Save the combined plot
        ggsave(filename = output_file, plot = combined_plot, width = 10, height = 8, dpi = 600)
        cat("  Saved:", output_file, "\n")
      }
    }
  }
  
  # Properly remove the Seurat object and clean up memory
  cat("Memory before clearing object:", pryr::mem_used(), "\n")
  rm(seurat_obj)
  gc()
  cat("Memory after clearing object:", pryr::mem_used(), "\n")
}
