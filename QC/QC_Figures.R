####################################
# Generating QC metrics for all objects
# 

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)

# Read the QC metrics CSV
qc_metrics <- read.csv("~/Roselab/Spatial/dietary_project/data/Raw_objects/qc_metrics_all_edited.csv", stringsAsFactors = FALSE)

# Define output folder for final figures and create it if it doesn't exist
qc_figs_folder <- "~/Roselab/Spatial/dietary_project/figures/QC_Figures/"
if (!dir.exists(qc_figs_folder)) {
  dir.create(qc_figs_folder, recursive = TRUE)
}

# Group by Tgen_ID and TMA (each group corresponds to one composite figure)
grouped_qc <- qc_metrics %>% group_by(Tgen_ID, TMA) %>% group_split()

for (group_df in grouped_qc) {
  
  # Get unique group identifiers (assuming one unique Tgen_ID and TMA per group)
  Tgen_ID <- unique(group_df$Tgen_ID)
  TMA <- unique(group_df$TMA)
  group_name <- paste(Tgen_ID, TMA, sep = "_")
  
  # Create a list to store spatial composite plots (one per row in the group)
  spatial_plots_list <- list()
  
  # Loop through each row in the group (should be 4 rows per group)
  for (i in seq_len(nrow(group_df))) {
    row_data <- group_df[i, ]
    label <- row_data$label
    file_path <- row_data$file_path
    
    # Load the Seurat object from the saved .rds file
    seurat_obj <- readRDS(file_path)
    
    # Create the first SpatialFeaturePlot (with dots and title)
    p1 <- SpatialFeaturePlot(
      seurat_obj,
      features = "nCount_Spatial.008um",
      pt.size.factor = 5,
      image.alpha = 0.6,
      max.cutoff = 400
    ) +
      scale_fill_viridis_c(option = "magma") +
      ggtitle(label) +
      theme_void() +
      theme(legend.position = "right")
    
    # Create the second SpatialFeaturePlot (only image, no dots)
    p2 <- SpatialFeaturePlot(
      seurat_obj,
      features = "nCount_Spatial.008um",
      pt.size.factor = 0,
      image.alpha = 1,
      max.cutoff = 500
    ) +
      theme(legend.position = "none")
    
    # Combine the two plots side by side using patchwork
    combined_spatial <- p1 + p2 + plot_layout(ncol = 2)
    spatial_plots_list[[i]] <- combined_spatial
  }
  
  # Arrange the 4 composite spatial plots in a 2x2 grid
  spatial_grid <- wrap_plots(spatial_plots_list, ncol = 2)
  
  # Create bar plots for QC metrics from group_df
  bp_fraction_under_100 <- ggplot(group_df, aes(x = label, y = fraction_under_100, fill = label)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Fraction Under 100")
  
  bp_fraction_over_100 <- ggplot(group_df, aes(x = label, y = fraction_over_100, fill = label)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Fraction Over 100")
  
  bp_mean_umi <- ggplot(group_df, aes(x = label, y = mean_umi, fill = label)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Mean UMI")
  
  bp_median_umi <- ggplot(group_df, aes(x = label, y = median_umi, fill = label)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Median UMI")
  
  bp_mean_percent_mt <- ggplot(group_df, aes(x = label, y = mean_percent_mt, fill = label)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Mean Percent MT")
  
  # Arrange the bar plots in a grid (here we use 3 columns; adjust layout as desired)
  qc_bar_grid <- wrap_plots(
    bp_mean_percent_mt,
    bp_fraction_under_100,
    bp_fraction_over_100,
    bp_median_umi,
    bp_mean_umi,
    ncol = 3
  )
  
  # Combine the spatial grid (top) with the QC bar plots grid (bottom)
  final_plot <- spatial_grid / qc_bar_grid +
    plot_annotation(title = group_name)
  
  # Define output file name based on group (Tgen_ID + TMA)
  out_file <- file.path(qc_figs_folder, paste0(group_name, ".png"))
  
  # Save the composite figure as a high-resolution PNG
  ggsave(filename = out_file, plot = final_plot, width = 16, height = 12, dpi = 300)
}

