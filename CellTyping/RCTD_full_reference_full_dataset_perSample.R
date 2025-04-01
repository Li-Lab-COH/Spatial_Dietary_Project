##########################################################
##########################################################
# Title: RCTD_full_reference_perSample
# Author: Jonathan Anzules
# Purpose: Generates spatial plots from the RCTD run.
#          This script uses the fully annotated reference
#          (as opposed to the condensed) and outputs:
#          - cell counts,
#          - Spatial plots,
#          - Test parameters,
#          - Uncertain cutoff values.
#
# Important: This script is intended to run one sample per job.
##########################################################

#############################
# Setup
#############################

library(Seurat)
library(spacexr)
library(ggplot2)
library(pryr)
set.seed(1337)

############################################
# Setting up test parameters - user defined
############################################

test_name <- "test_2"
gene_cutoff <- 0.000125
fc_cutoff <- 0.5
confidence_threshold <- 5
doublet_threshold <- 20

############################################
# Setting up functions
############################################

# Show memory output in a human readable format
pretty_mem <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB")
  power <- ifelse(bytes > 0, floor(log(bytes, 1024)), 0)
  power <- min(power, length(units) - 1)
  converted <- bytes / (1024^power)
  sprintf("%.2f %s", converted, units[power + 1])
}

log_block <- function(...) {
  msg <- paste(rep("=", 88), collapse = "")
  message(msg)
  message("==> ", ...)
  message(msg)
}


ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

############################################
# Setting up directories and aggregated output files
############################################

results_dir <- "/home/janzules/Spatial/dietary_project/data/RCTD_param_tests/CellTypeFullAnnotationTest"
test_dir <- file.path(results_dir, test_name)
ensure_dir(test_dir)

# Folder for spatial plots (all samples go here)
spatial_plot_dir <- file.path(test_dir, "spatial_plots")
ensure_dir(spatial_plot_dir)

# Base folder for confidence histograms (each sample gets its own subfolder)
hist_base_dir <- file.path(test_dir, "confidence_histograms")
ensure_dir(hist_base_dir)

# Define paths for aggregated output files; these will be appended if they exist.
cell_counts_file <- file.path(test_dir, paste0("Cell_counts_", test_name, ".csv"))
cutoff_file <- file.path(test_dir, paste0("myccap_cutoff_", test_name, ".csv"))

# Save test parameters for record
param_lines <- c(
  paste0("test_name: ", test_name),
  paste0("gene_cutoff: ", gene_cutoff),
  paste0("fc_cutoff: ", fc_cutoff),
  paste0("confidence_threshold: ", confidence_threshold),
  paste0("doublet_threshold: ", doublet_threshold)
)

writeLines(param_lines, file.path(test_dir, "rctd_params.txt"))

############################################
# Preparing important variables for aggregated data
############################################

# Define expected cell types (using underscores) with the proper commas!
expected_cell_types <- c(
  "MyC_CaP",
  "uncertain_MyC_CaP",
  "Dendritic_Cells",
  "ISG_high_Macrophages",
  "Proliferating_Macrophages",
  "M2_Macrophages",
  "Macrophages",
  "T_Cells",
  "NK_Cells",
  "Tregs",
  "LumP",
  "Mesenchymal",
  "Endothelial_Cells",
  "Fibroblasts"
)

# Initialize master data frames for cell counts and MyC-CaP cutoff values
label_df <- data.frame(row.names = expected_cell_types)
cutoff_df <- data.frame(sample = character(), myccap_cutoff = numeric(), stringsAsFactors = FALSE)

# Prepare color palette and label order for final spatial plotting using "myccap_strict"
my_labels <- c(
  "MyC_CaP",
  "uncertain_MyC_CaP",
  "Dendritic_Cells",
  "ISG_high_Macrophages",
  "Proliferating_Macrophages",
  "M2_Macrophages",
  "Macrophages",
  "T_Cells",
  "NK_Cells",
  "Tregs",
  "LumP",
  "Mesenchymal",
  "Endothelial_Cells",
  "Fibroblasts",
  "Unknown"
)

my_colors <- c(
  `MyC_CaP`                 = "black",
  `uncertain_MyC_CaP`       = "grey80",
  Dendritic_Cells           = "steelblue1",
  ISG_high_Macrophages      = "lightskyblue3",
  Proliferating_Macrophages = "dodgerblue3",
  M2_Macrophages            = "deepskyblue4",
  Macrophages               = "blue4",
  T_Cells                   = "firebrick1",
  NK_Cells                  = "darkred",
  Tregs                     = "violetred3",
  LumP                      = "olivedrab",
  Mesenchymal               = "darkolivegreen4",
  Endothelial_Cells         = "yellowgreen",
  Fibroblasts               = "chartreuse4",
  Unknown                   = "white"
)
# Order the color vector by my_labels:
my_colors <- my_colors[my_labels]

############################################
# Creating Reference (outside the loop)
############################################

MyC_CaP_ref <- readRDS("/home/janzules/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_unal_reference.rds")
Idents(MyC_CaP_ref) <- "celltype"
# (Optional: remove unwanted cell types if needed)
counts_ref <- MyC_CaP_ref[["RNA"]]$counts
cluster <- as.factor(MyC_CaP_ref$celltype)
nUMI <- MyC_CaP_ref$nCount_RNA
levels(cluster) <- gsub(" ", "_", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts_ref, cluster, nUMI)

############################################
# Job Array: Process one sample per job using SLURM_ARRAY_TASK_ID
############################################

args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])  # Should be 0 to (n_samples - 1)
file_list <- readLines("/home/janzules/Spatial/dietary_project/code/addresses/Myc_CaP_raw.txt")
file <- file_list[array_id + 1]
sample_name <- gsub("\\.rds$", "", basename(file))
message("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
message("Processing file: ", basename(file), " -> Sample name: ", sample_name)
message("Current memory usage: ", pretty_mem(pryr::mem_used()))
message("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# Load the Seurat object for this sample
prostate_ST <- readRDS(file)
prostate_ST <- subset(
  prostate_ST,
  subset = nCount_Spatial.008um >= 100 &
           nFeature_Spatial.008um >= 100 &
           percent.mt <= 15
)

#############################
# Setup output directories for this sample
#############################
hist_folder <- file.path(hist_base_dir, sample_name)
ensure_dir(hist_folder)
message("Created histogram directory: ", hist_folder)

#############################
# Run PCA on full dataset (using Spatial.008um assay)
#############################
log_block("Processing full dataset")
DefaultAssay(prostate_ST) <- "Spatial.008um"

# Scale the data — required before PCA
prostate_ST <- NormalizeData(prostate_ST)
prostate_ST <- FindVariableFeatures(prostate_ST)
prostate_ST <- ScaleData(prostate_ST, verbose = TRUE)

prostate_ST <- RunPCA(
  object = prostate_ST,
  assay = "Spatial.008um",
  reduction.name = "pca.prostate.full",
  verbose = TRUE
)

log_block("Completed PCA. Current memory usage: ", pretty_mem(pryr::mem_used()))

prostate_ST <- FindNeighbors(prostate_ST, assay = "Spatial.008um",
              reduction = "pca.prostate.full", dims = 1:15)

log_block("Completed Find Neighbors. Current memory usage: ", pretty_mem(pryr::mem_used()))


#############################
# Create Query Object and Run RCTD for this sample
#############################
counts_hd <- prostate_ST[["Spatial.008um"]]$counts
prostate_cells_hd <- colnames(prostate_ST[["Spatial.008um"]])

coords <- GetTissueCoordinates(prostate_ST)
coords <- coords[intersect(rownames(coords), prostate_cells_hd), 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# log_memory_and_cpu(label = paste("RCTD:", sample_name), interval_sec = 120, repetitions = 240)

log_block("Creating RCTD object... Current memory usage: ", pretty_mem(pryr::mem_used()))

RCTD <- create.RCTD(query, reference, max_cores = 12, UMI_min = 100,
                    gene_cutoff = gene_cutoff,
                    fc_cutoff = fc_cutoff,
                    CONFIDENCE_THRESHOLD = confidence_threshold,
                    DOUBLET_THRESHOLD = doublet_threshold,
                    CELL_MIN_INSTANCE = 24)

log_block("Running RCTD... Current memory usage: ", pretty_mem(pryr::mem_used()))
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
log_block("RCTD run complete! Current memory usage: ", pretty_mem(pryr::mem_used()))


#############################
# Project RCTD Labels from the Sketched Cells to the Full Dataset
#############################

# Add RCTD results to metadata
prostate_ST <- AddMetaData(prostate_ST, metadata = RCTD@results$results_df)
prostate_ST$first_type <- as.character(prostate_ST$first_type)
prostate_ST$first_type[is.na(prostate_ST$first_type)] <- "Unknown"


#############################
# Labeling Uncertain MyC-CaP Bins Based on Confidence for this sample
#############################
myccap_scores <- prostate_ST$full_first_type.score[prostate_ST$full_first_type == "MyC_CaP"]
myccap_cutoff <- quantile(myccap_scores, probs = 0.15, na.rm = TRUE)
log_block("MyC-CaP cutoff for ", sample_name, ": ", myccap_cutoff)
cutoff_df <- rbind(cutoff_df, data.frame(sample = sample_name, myccap_cutoff = myccap_cutoff, stringsAsFactors = FALSE))
prostate_ST$myccap_strict <- as.character(prostate_ST$first_type)
prostate_ST$myccap_strict[
  prostate_ST$first_type == "MyC_CaP" & prostate_ST$full_first_type.score < myccap_cutoff
] <- "uncertain_MyC_CaP"
print(table(prostate_ST$myccap_strict))


#############################
# Saving annotated object
#############################
saveRDS(prostate_ST, file = file.path("/home/janzules/Spatial/dietary_project/data/RCTD_annotated_n_PCA_full", paste0(sample_name, "RCTD_annotated.rds")))


#############################
# Aggregate Cell Count Data for this sample
#############################
label_counts <- table(prostate_ST$myccap_strict)
label_vector <- as.numeric(label_counts[match(expected_cell_types, names(label_counts))])
names(label_vector) <- expected_cell_types
label_vector[is.na(label_vector)] <- 0
label_df[[sample_name]] <- label_vector


#############################
# Generate and Save Histograms for QC (for each cell type except uncertain_MyC-CaP)
#############################
cell_types_hist <- setdiff(unique(prostate_ST$full_first_type), "uncertain_MyC_CaP")
for (ct in cell_types_hist) {
  scores <- prostate_ST$full_first_type.score[prostate_ST$full_first_type == ct]
  file_name <- paste0(gsub("[^A-Za-z0-9_]", "_", ct), ".png")
  file_path <- file.path(hist_folder, file_name)
  png(filename = file_path, width = 800, height = 600, res = 150)
  hist(scores, main = paste("Confidence Scores for", ct),
       xlab = "full_first_type.score", breaks = 50, col = "skyblue", border = "white")
  dev.off()
}

#############################
# Generate High-Resolution Spatial Plot using myccap_strict Labels
#############################
prostate_ST$myccap_strict <- factor(prostate_ST$myccap_strict, levels = my_labels)
p1_strict <- SpatialDimPlot(
  object = prostate_ST,
  group.by = "myccap_strict",
  label = FALSE,
  repel = TRUE,
  image.alpha = 0.9,
  pt.size.factor = 3,
  crop = TRUE,
  cols = my_colors
)
spatial_file <- file.path(spatial_plot_dir, paste0(sample_name, "_spatial_celltype.png"))
ggsave(filename = spatial_file, plot = p1_strict, width = 8, height = 6, dpi = 500)

#############################
# Save the Processed RDS Object for this sample
#############################
# rds_file <- file.path(labeled_rds_dir, paste0(sample_name, ".rds"))
# saveRDS(prostate_ST, file = rds_file)
message("Finished processing sample: ", sample_name)

# Clean up memory: remove heavy slots and run gc
message("Memory before clearing object: ", pretty_mem(pryr::mem_used()))
prostate_ST@reductions <- list()
prostate_ST@assays <- list()
prostate_ST@images <- list()
prostate_ST@meta.data <- data.frame()
prostate_ST@graphs <- list()
prostate_ST@neighbors <- list()
gc(full = TRUE)
rm(prostate_ST)
message("Memory after clearing object: ", pretty_mem(pryr::mem_used()))
message("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
message("Job for sample ", sample_name, " completed successfully.")

############################################
# Append Aggregated Data to Files (if desired)
############################################
if (file.exists(cell_counts_file)) {
  old_label_df <- read.csv(cell_counts_file, row.names = 1)
  new_label_df <- cbind(old_label_df, label_df)
  write.csv(new_label_df, cell_counts_file, row.names = TRUE)
} else {
  write.csv(label_df, cell_counts_file, row.names = TRUE)
}

if (file.exists(cutoff_file)) {
  old_cutoff_df <- read.csv(cutoff_file, stringsAsFactors = FALSE)
  new_cutoff_df <- rbind(old_cutoff_df, cutoff_df)
  write.csv(new_cutoff_df, cutoff_file, row.names = FALSE)
} else {
  write.csv(cutoff_df, cutoff_file, row.names = FALSE)
}

message("Saved aggregated cell counts to: ", cell_counts_file)
message("Saved MyC-CaP cutoff values to: ", cutoff_file)
message("Job for sample ", sample_name, " completed successfully.")
