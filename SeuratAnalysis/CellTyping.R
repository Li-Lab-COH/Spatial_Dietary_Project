library(Seurat)
library(dplyr)
library(tibble)
library(pryr)  # for mem_used()
library(future)

# Prevent unwanted parallelization
plan("sequential")

# Define your tumor clusters list for the samples of interest.
# NOTE: If you have duplicate names in your list, you may want to combine the clusters (e.g. union).
tumor_clusters_list <- list(
  "F07833_5_RT"    = c(1, 6, 8, 9),
  "F07834_28_LFRT" = c(3, 16, 7),  # combined clusters if duplicates occur
  "F07833_8_CRRT"  = c(1, 13),
  "F07833_9_KDRT"  = c(2, 10),
  "F07834_5_RT"    = c(3, 6, 9),
  "F07834_8_CRRT"  = c(3, 4),
  "F07834_9_KDRT"  = c(3, 4, 0),
  "F07835_7_CRRT"  = c(2, 5, 7, 9),
  "F07836_11_RT"   = c(8, 9),
  "F07837_31_RT"   = c(0, 2, 6, 8, 9),
  "F07837_37_CR_RT"= c(8, 9, 11),
  "F07838_32_RT"   = c(2, 4, 10, 11),
  "F07838_36_CR_RT"= c(4, 5),
  "F07838_38_CR_RT"= c(2, 7)
)
# tumor_clusters_list <- tumor_clusters_list[0:2]
# Define the expected cell types for the output table.
expected_cell_types <- c("B cell", "Cd4", "Cd8", "Dendritic cell", "NK cell", 
                         "N1", "N2", "M2", "Suppressor", "Treg", "Tumor", "Other")

# Specify the input directory containing your .rds files.
input_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/data/Rogelio/Analysis_031625/BANKSY_Normalized_QC_Filtered"

# Initialize a master counts data frame.
master_counts <- data.frame(Cell_Type = expected_cell_types)

# Loop over each sample listed in the tumor_clusters_list.
for (sample_name in names(tumor_clusters_list)) {
  
  # Construct the file path. Assumes files are named as "sample_name.rds".
  file_path <- file.path(input_dir, paste0(sample_name, ".rds"))
  
  if (!file.exists(file_path)) {
    cat("File not found for sample", sample_name, "\n")
    next
  }
  
  cat("Processing sample:", sample_name, "\n")
  seurat_obj <- readRDS(file_path)
  
  ### 1. Tumor Cell Identification ###
  DefaultAssay(seurat_obj) <- "BANKSY"
  
  # Use the clusters from the tumor_clusters_list (converted to character).
  target_clusters <- as.character(tumor_clusters_list[[sample_name]])
  cells_to_check <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$BANKSY_snn_res.0.5 %in% target_clusters, ])
  
  tumor_gene_list <- c("Myc", "Ar", "Ccna2", "Ccnd1", "Cdk1", "Odc1", "Slc7a5", "Ldha", "Cad")
  expr_data <- FetchData(seurat_obj, vars = tumor_gene_list, slot = "data")
  colnames(expr_data) <- gsub("^spatial008um_", "", colnames(expr_data))
  
  expr_data_subset <- expr_data[cells_to_check, ]
  
  tumor_cells <- expr_data_subset %>%
    filter(
      (Myc > 0 | Ar > 0) &  # at least one of these markers
        rowSums(.[, c("Ccna2", "Ccnd1", "Cdk1", "Odc1", "Slc7a5", "Ldha", "Cad")] > 0) >= 1
    ) %>%
    rownames()
  
  # Label all cells as "Other" then mark tumor cells.
  seurat_obj@meta.data$cell_types <- "Other"
  seurat_obj@meta.data[tumor_cells, "cell_types"] <- "Tumor"
  
  ### 2. Immune Cell Classification ###
  immune_gene_list <- c(
    "Trac", "Trbc1", "Trbc2", "Cd3e", "Cd4", "Cd8a", "Foxp3",      # T cell markers
    "Cd19", "Cd22", "Ighm",                                         # B cell markers
    "Itgax", "Batf3", "H2-Ab1",                                     # Dendritic cell markers
    "Ncr1",                                                         # NK cell marker
    "Cxcr3", "Il12b", "Nos2", "Tnf", "Ccl5",                        # N1 neutrophil markers
    "Tgfb1", "Vegfa", "Arg1", "Cd274", "Mmp9",                      # Pro-tumor markers for neutrophils
    "Retnlg", "Cxcr2",                                              # N2 markers
    "Cd163", "Mrc1", "Il10", "Ccl22"                                 # M2 markers
  )
  
  expr_data_ext <- FetchData(seurat_obj, vars = immune_gene_list, slot = "data")
  colnames(expr_data_ext) <- gsub("^spatial008um_", "", colnames(expr_data_ext))
  
  non_tumor_ids <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$cell_types == "Other"]
  immune_expr <- expr_data_ext[non_tumor_ids, ] %>% rownames_to_column("cell")
  
  immune_expr <- immune_expr %>%
    rowwise() %>%
    mutate(
      immune_type = case_when(
        (Foxp3 > 0) ~ "Treg",
        (Cd4 > 0) ~ "Cd4",
        (Cd8a > 0) ~ "Cd8",
        (Cd19 > 0) ~ "B cell",
        (Itgax > 0) ~ "Dendritic cell",
        (Ncr1 > 0.5) ~ "NK cell",
        (Cxcr3 > 1 & Il12b > 1 & Nos2 > 1.5 & Tnf > 1 & Ccl5 > 1) ~ "N1",
        (( (Tgfb1 > 2) + (Vegfa > 1.75) + (Arg1 > 1.75) + (Cd274 > 1.75) + (Mmp9 > 1.5) ) >= 2 &
           (Retnlg > 1 | Cxcr2 > 1.75)) ~ "N2",
        (( (Tgfb1 > 2) + (Vegfa > 1.75) + (Arg1 > 1.75) + (Cd274 > 1.75) + (Mmp9 > 1.5) ) >= 2 &
           ((Cd163 > 1.75) + (Mrc1 > 2) + (Il10 > 1.75) + (Ccl22 > 1)) >= 2) ~ "M2",
        (( (Tgfb1 > 2) + (Vegfa > 1.75) + (Arg1 > 1.75) + (Cd274 > 1.75) + (Mmp9 > 1.5) ) >= 3) ~ "Suppressor",
        TRUE ~ "Other"
      )
    ) %>%
    ungroup()
  
  immune_assignments <- immune_expr %>% select(cell, immune_type)
  seurat_obj@meta.data$cell_types[match(immune_assignments$cell, rownames(seurat_obj@meta.data))] <- immune_assignments$immune_type
  
  Idents(seurat_obj) <- seurat_obj@meta.data$cell_types
  
  ### 3. Generate and Store Cell Count Data ###
  cell_counts <- as.data.frame(table(seurat_obj@meta.data$cell_types))
  colnames(cell_counts) <- c("Cell_Type", "Count")
  
  # Ensure all expected cell types are included
  all_cell_counts <- data.frame(Cell_Type = expected_cell_types)
  all_cell_counts <- merge(all_cell_counts, cell_counts, by = "Cell_Type", all.x = TRUE)
  all_cell_counts$Count[is.na(all_cell_counts$Count)] <- 0
  
  # Rename the 'Count' column with the current sample name
  sample_counts <- all_cell_counts %>% select(Cell_Type, Count)
  colnames(sample_counts)[2] <- sample_name
  
  # Merge into the master_counts data frame
  master_counts <- merge(master_counts, sample_counts, by = "Cell_Type", all = TRUE)
  
  ### 4. Clean Up Memory ###
  cat("Memory before clearing object for", sample_name, ":", pryr::mem_used(), "\n")
  
  # Explicitly remove all large objects
  rm(seurat_obj, expr_data, expr_data_subset, expr_data_ext, immune_expr, immune_assignments, cell_counts, all_cell_counts)
  
  # Aggressive garbage collection
  gc(reset = TRUE)
  invisible(gc())
  
  cat("Memory after clearing object for", sample_name, ":", pryr::mem_used(), "\n\n")
}

# Print the final master cell counts data frame.
print(master_counts)
