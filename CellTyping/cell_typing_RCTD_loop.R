##########################
# Author: Jonathan M. Anzules
# e-mail: jonanzule@gmail.com
# purpose: Testing various parameters for cell type annotation, exporting, cell counts, visuals, and parameters used


# Load required libraries
library(Seurat)
library(devtools)
library(spacexr)
library(ggplot2)
set.seed(1337)

################################
# Setting up Testing parameters
################################
args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])  # This is your SLURM_ARRAY_TASK_ID

param_df <- read.csv("/home/janzules/Spatial/dietary_project/data/cell_typing_reference/rctd_parameter_grid.csv",
               row.names = 1, check.names = FALSE)
# Assign parameter values
test_name <- colnames(param_df)[array_id + 1]
gene_cutoff <- as.numeric(param_df["gene_cutoff", array_id +1])
fc_cutoff <- as.numeric(param_df["fc_cutoff", array_id +1])
confidence_threshold <- as.numeric(param_df["confidence_threshold", array_id +1])
doublet_threshold <- as.numeric(param_df["doublet_threshold", array_id +1])


# Example: Print or use the parameters
cat("Running test", test_name, "\n")
cat("gene_cutoff:", gene_cutoff, "\n")
cat("fc_cutoff:", fc_cutoff, "\n")
cat("confidence_threshold:", confidence_threshold, "\n")
cat("doublet_threshold:", doublet_threshold, "\n")



# Create a named character vector of parameters
param_lines <- c(
  paste0("test_name: ", colnames(param_df)[array_id + 1]),
  paste0("gene_cutoff: ", gene_cutoff),
  paste0("fc_cutoff: ", fc_cutoff),
  paste0("confidence_threshold: ", confidence_threshold),
  paste0("doublet_threshold: ", doublet_threshold)
)

################################
# Setting up directories
################################
# Major folder forl all the saved
results_dir <- "/home/janzules/Spatial/dietary_project/data/RCTD_param_tests"

# Test specific folders
test_dir <- file.path(results_dir, test_name)
spatial_plot_dir <- file.path(test_dir, "spatial_plots")
hist_base_dir <- file.path(test_dir, "confidence_histograms")

# Function to create directories
ensure_dir <- function(dir_path){
  if (!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
}

# Create directories
ensure_dir(test_dir)
ensure_dir(spatial_plot_dir)
ensure_dir(hist_base_dir)

# saving parameters used for this test
writeLines(param_lines, paste0(test_dir, "rctd_params", ".txt"))


################################
# Preparing important variables
################################

# Define expected cell types (using underscores)
expected_cell_types <- c(
  "MyC-CaP",
  "uncertain_MyC-CaP",
  "Mesenchymal",
  "Macrophages",
  "Endothelial_Cells",
  "NK_n_Tcells",
  "Dendritic_Cells",
  "Fibroblasts"
)

# Initialize master data frames for cell counts and MyC-CaP cutoff values
label_df <- data.frame(row.names = expected_cell_types)
cutoff_df <- data.frame(sample = character(), myccap_cutoff = numeric(), stringsAsFactors = FALSE)

# Prepare color palette and label order for final spatial plotting using the "myccap_strict" label:
my_labels <- c(
  "MyC-CaP",
  "uncertain_MyC-CaP",
  "Macrophages",
  "NK_n_Tcells",
  "Dendritic_Cells",
  "Endothelial_Cells",
  "Mesenchymal",
  "Fibroblasts",
  "Unknown"  # optional for unmatched cells
)


my_colors <- c(
  `MyC-CaP`            = "black",
  `uncertain_MyC-CaP`  = "grey80",
  Macrophages          = "blue2",
  NK_n_Tcells          = "firebrick1",
  Dendritic_Cells      = "steelblue1",
  Endothelial_Cells    = "darkolivegreen4",
  Mesenchymal          = "gold",
  Fibroblasts          = "darkgoldenrod2",
  Unknown              = "white"
)

# Ensure the color vector is in the same order as `my_labels`
my_colors <- my_colors[my_labels]

################################
# Creating Reference
################################
# Load the scRNA-seq reference object once
# TODO: move this to the where all the reference building is happening
MyC_CaP_ref <- readRDS("/home/janzules/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_condensed.rd")
# Prepare the reference object (only once, outside loop)
Idents(MyC_CaP_ref) <- "celltype"
# Remove unwanted cell types (e.g., "Resident Macrophages")
# cells_to_keep <- rownames(MyC_CaP_ref@meta.data)[!(MyC_CaP_ref$celltype %in% c("Resident Macrophages"))]
# MyC_CaP_ref <- subset(MyC_CaP_ref, cells = cells_to_keep)
# MyC_CaP_ref$celltype <- droplevels(MyC_CaP_ref$celltype)
counts_ref <- MyC_CaP_ref[["RNA"]]$counts
cluster <- as.factor(MyC_CaP_ref$celltype)
nUMI <- MyC_CaP_ref$nCount_RNA
levels(cluster) <- gsub(" ", "_", levels(cluster))
cluster <- droplevels(cluster)
# Create RCTD reference object once:
reference <- Reference(counts_ref, cluster, nUMI)

# Get list of RDS files to process
file_list <- list.files(path = "~/Roselab/Spatial/dietary_project/data/Rogelio/Analysis_final", 
                        pattern = "\\.rds$", full.names = TRUE)

#############################
# MAIN LOOP: Process each RDS file
#############################

for (file in file_list) {
  
  # Extract sample name from file path
  sample_name <- gsub("\\.rds$", "", basename(file))
  message("Processing file: ", basename(file), " -> Sample name: ", sample_name)
  
  # Load the Seurat object for this sample
  prostate_ST <- readRDS(file)
  prostate_ST <- subset(prostate_ST, subset = nCount_Spatial.008um >= 100)
  
  #############################
  # Setup output directories for this sample
  #############################
  hist_folder <- file.path(hist_base_dir, sample_name)
  ensure_dir(hist_folder)
  
  message("Created histogram directory: ", hist_folder)
  
  #############################
  # Run PCA on full dataset (using Spatial.008um assay)
  #############################
  message("Starting PCA on full dataset")

  DefaultAssay(prostate_ST) <- "Spatial.008um"
  prostate_ST <- RunPCA(
    object = prostate_ST,
    assay = "Spatial.008um",
    reduction.name = "pca.prostate.full",
    verbose = TRUE
  )
  message("PCA of full ST object complete...")
  #############################
  # Sketching (subsample using LeverageScore)
  #############################
  message("Starting sketch processing...")
 
  DefaultAssay(prostate_ST) <- "Spatial.008um"
  prostate_ST <- SketchData(
    object = prostate_ST,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  
  DefaultAssay(prostate_ST) <- "sketch"
  prostate_ST <- ScaleData(prostate_ST)
  prostate_ST <- RunPCA(prostate_ST, assay = "sketch",
                        reduction.name = "pca.prostate.sketch",
                        verbose = TRUE)
  # ElbowPlot(prostate_ST, reduction = "pca.prostate.sketch")
  
  prostate_ST <- FindNeighbors(prostate_ST, reduction = "pca.prostate.sketch",
                               dims = 1:15)
  prostate_ST <- RunUMAP(prostate_ST, reduction = "pca.prostate.sketch",
                         reduction.name = "umap.prostate.sketch",
                         return.model = TRUE,
                         dims = 1:15,
                         verbose = TRUE)
  
  message("Sketch processing complete...")
  #############################
  # Create Query Object and Run RCTD for this sample
  #############################

  counts_hd <- prostate_ST[["sketch"]]$counts
  prostate_cells_hd <- colnames(prostate_ST[["sketch"]])
  coords <- GetTissueCoordinates(prostate_ST)[prostate_cells_hd, 1:2]
  
  query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
  
  RCTD <- create.RCTD(query, reference, max_cores = 32, UMI_min = 100,
                        gene_cutoff = gene_cutoff
                        fc_cutoff = fc_cutoff
                        CONFIDENCE_THRESHOLD = confidence_threshold
                        DOUBLET_THRESHOLD = doublet_threshold
              )
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # Add RCTD results to metadata
  prostate_ST <- AddMetaData(prostate_ST, metadata = RCTD@results$results_df)
  
  #############################
  # Project RCTD labels from the sketched cells to the full dataset
  #############################
  message("Query complete, projecting...")
  
  prostate_ST$first_type <- as.character(prostate_ST$first_type)
  prostate_ST$first_type[is.na(prostate_ST$first_type)] <- "Unknown"
  
  prostate_ST <- ProjectData(
    object = prostate_ST,
    assay = "Spatial.008um",
    full.reduction = "pca.prostate.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca.prostate.sketch",
    umap.model = "umap.prostate.sketch",
    dims = 1:15,
    refdata = list(full_first_type = "first_type")
  )
  
  #############################
  # Labeling uncertain MyC-CaP bins based on confidence for this sample
  #############################
  myccap_scores <- prostate_ST$full_first_type.score[prostate_ST$full_first_type == "MyC-CaP"]
  myccap_cutoff <- quantile(myccap_scores, probs = 0.15, na.rm = TRUE)
  message("MyC-CaP cutoff for ", sample_name, ": ", myccap_cutoff)
  
  # Save cutoff value to master cutoff_df
  cutoff_df <- rbind(cutoff_df, data.frame(sample = sample_name, myccap_cutoff = myccap_cutoff, stringsAsFactors = FALSE))
  
  # Create uncertain label: convert first_type to character first
  prostate_ST$myccap_strict <- as.character(prostate_ST$first_type)
  prostate_ST$myccap_strict[
    prostate_ST$first_type == "MyC-CaP" & prostate_ST$full_first_type.score < myccap_cutoff
  ] <- "uncertain_MyC-CaP"
  
  # Check label counts for debugging
  print(table(prostate_ST$myccap_strict))
  
  #############################
  # Aggregate Cell Count Data for this sample
  #############################
  label_counts <- table(prostate_ST$myccap_strict)
  # Align label_counts with expected_cell_types
  label_vector <- as.numeric(label_counts[expected_cell_types])
  names(label_vector) <- expected_cell_types
  # Replace any NA with 0
  label_vector[is.na(label_vector)] <- 0
  
  # Add this column to the master label_df (column name = sample_name)
  label_df[[sample_name]] <- label_vector
  
  #############################
  # Generate and Save Histograms for QC (for each cell type except uncertain_MyC-CaP)
  #############################
  cell_types_hist <- setdiff(unique(prostate_ST$full_first_type), "uncertain_MyC-CaP")
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
  # Generate High-Resolution Spatial Plot using myccap_strict labels
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
  
  spatial_file <- file.path(spatial_plot_dir, paste0(sample_name, ".png"))
  ggsave(filename = spatial_file, plot = p1_strict, width = 8, height = 6, dpi = 500)
  
  
  message("Finished processing sample: ", sample_name)
  
  # Properly remove the Seurat object and clean up memory
  cat("Memory before clearing object:", pryr::mem_used(), "\n")
  prostate_ST@reductions <- list()
  prostate_ST@assays <- list()
  prostate_ST@images <- list()
  prostate_ST@meta.data <- data.frame()
  prostate_ST@graphs <- list()
  prostate_ST@neighbors <- list()
  gc(full = TRUE)
  cat("Memory after clearing object:", pryr::mem_used(), "\n")
}


#############################
# After the loop: Save the aggregated data frames with date-time stamps
#############################
cell_counts_file <- file.path(test_dir, paste0("Cell_counts_", test_name, ".csv"))
cutoff_file <- file.path(test_dir, paste0("myccap_cutoff_", test_name, ".csv"))

#  TODO: Make sure the right locations are being used here
write.csv(label_df, file = cell_counts_file, row.names = TRUE)
write.csv(cutoff_df, file = cutoff_file, row.names = FALSE)

message("Saved aggregated cell counts to: ", cell_counts_file)
message("Saved MyC-CaP cutoff values to: ", cutoff_file)
