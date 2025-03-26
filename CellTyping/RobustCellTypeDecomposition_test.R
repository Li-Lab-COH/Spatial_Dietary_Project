##########################
# Author: Jonathan M. Anzules
# e-mail: jonanzule@gmail.com
# purpose: testing out robust cell type deconvolution and then will restructure
# as a loop for all the objects.

library(Seurat)
library(devtools)
library(spacexr)
library(ggplot2)
set.seed(1337)

#-------------------------------- loading data -------------------------------
MyC_CaP_ref <- readRDS("~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_unal_reference.rds")


#----------------------------- Initializing variables -----------------------------
# labeled RDS
labeled_rds_dir <- "~/Roselab/Spatial/dietary_project/data/CellTyping_objects/"
figure_dir <- "~/Roselab/Spatial/dietary_project/figures/RCTD/"

# Define expected cell types
expected_cell_types <- c(
  "Dendritic_Cells", "Endothelial_Cells", "Fibroblasts", "ISG-high_Macrophages",
  "LumP", "M2_Macrophages", "Macrophages", "Mesenchymal", "MyC-CaP", "uncertain_MyC-CaP",
  "Proliferating_Macrophages", "Unknown"
)

# Initialize master label_df
label_df <- data.frame(row.names = expected_cell_types)

#------------------------ Within loop variables ------------------------------
# # This is how we starting the loop
# for (file_path in file_list) {
#   
#   # Extract file name without path
#   file_name <- basename(file_path)
#   
#   # Remove ".rds" extension
#   sample_name <- gsub("\\.rds$", "", file_name)
#   
#   # Print extracted name for debugging
#   message("Processing file: ", file_name, " -> Sample name: ", sample_name)
#   
#   # Load the Seurat object
#   prostate_ST <- readRDS(file_path)

# The next two variables should be dynamically prepared in the loop
file <- "~/Roselab/Spatial/dietary_project/data/Rogelio/Analysis_final/F07833_28_LFRT.rds"
sample_name <- gsub("\\.rds$", "", basename(file))

prostate_ST <- readRDS(file)

hist_folder <- file.path(figure_dir, "confidence_histograms", sample_name)
# Create output directory if it does not exist
if (!dir.exists(hist_folder)) {
  dir.create(hist_folder, recursive = TRUE)
  message("Created directory: ", hist_folder)
}

spatial_plot <- file.path(figure_dir, "spatial_plot", sample_name)
if (!dir.exists(spatial_plot)) {
  dir.create(spatial_plot, recursive = TRUE)
  message("Created directory: ", spatial_plot)
}

# BELOW EVERYTHING SHOULD BE IN THE LOOP
#----------------------- Running PCA on full dataset -------------------------
# Run PCA command on the entire object, needs to be done for all objects
# dims will be kept at a constant 1:15
DefaultAssay(prostate_ST) <- "Spatial.008um"
prostate_ST <- RunPCA(
  object = prostate_ST,
  assay = "Spatial.008um",
  reduction.name = "pca.prostate.full",
  verbose = T
)




#--------------------------- Sketching -----------------------------------
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
                      verbose = T)
ElbowPlot(prostate_ST, reduction = "pca.prostate.sketch")

prostate_ST <- FindNeighbors(prostate_ST, reduction = "pca.prostate.sketch",
                             dims = 1:15)
prostate_ST <- RunUMAP(prostate_ST, reduction = "pca.prostate.sketch",
                       reduction.name = "umap.prostate.sketch",
                       return.model = T,
                       dims = 1:15,
                       verbose = T)


#------------------------ Prepping Reference data ---------------------------
Idents(MyC_CaP_ref) <-  "celltype"


# Removing unlabeled and low cell numbers
cells_to_keep <- rownames(MyC_CaP_ref@meta.data)[
  !(MyC_CaP_ref$celltype %in% c("Resident Macrophages"))
]

MyC_CaP_ref <- subset(MyC_CaP_ref, cells = cells_to_keep)
MyC_CaP_ref$celltype <- droplevels(MyC_CaP_ref$celltype)
table(MyC_CaP_ref$celltype)

counts <- MyC_CaP_ref[["RNA"]]$counts
cluster <- as.factor(MyC_CaP_ref$celltype)
nUMI <- MyC_CaP_ref$nCount_RNA
levels(cluster) <- gsub(" ", "_", levels(cluster))
droplevels(cluster)

# creating RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- prostate_ST[["sketch"]]$counts
prostate_cells_hd <- colnames(prostate_ST[["sketch"]])
coords <- GetTissueCoordinates(prostate_ST)[prostate_cells_hd, 1:2]

#------------------------ Querying and Projecting -------------------------
# Create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
#return results
prostate_ST <- AddMetaData(prostate_ST, metadata = RCTD@results$results_df)


#Projecting RCTD labels from sketched to the rest of the data
prostate_ST$first_type <- as.character(prostate_ST$first_type)
prostate_ST$first_type[is.na(prostate_ST$first_type)] <- "Unknown"

# Projecting data
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


#---------------------- Labeling uncertain MyC-CaP ------------------------
# Labeling uncertain MyCaps

# Subset scores only for MyC-CaP bins
myccap_scores <- prostate_ST$full_first_type.score[prostate_ST$full_first_type == "MyC-CaP"]

# Dynamically calculate the 25th percentile
myccap_cutoff <- quantile(myccap_scores, probs = 0.15, na.rm = TRUE)

# Print or use it
print(myccap_cutoff)

# Default to keeping existing label
prostate_ST$myccap_strict <- as.character(prostate_ST$first_type)

# Relabel low-confidence MyC-CaP bins as "uncertain_MyC-CaP"
prostate_ST$myccap_strict[
  prostate_ST$first_type == "MyC-CaP" &
    prostate_ST$full_first_type.score < myccap_cutoff
] <- "uncertain_MyC-CaP"


saveRDS(prostate_ST, file_path(labeled_rds_dir, sample_name, ".rds"))
#------------------------------- Getting QC metrics --------------------------

label_counts <- table(prostate_ST$myccap_strict)
# Convert to a named vector for easy alignment
label_vector <- as.numeric(label_counts[expected_cell_types])  # will insert NAs where types are missing
names(label_vector) <- expected_cell_types

# Add this as a column to your master dataframe
label_df[[sample_name]] <- label_vector
# save df at this variable at the end of the script: labeled_rds_dir


# PROBABILITY SCORES
# Get unique identities
cell_types <- unique(prostate_ST$full_first_type)

# Loop through and save histograms
for (ct in cell_types) {
  # Subset scores for this cell type
  scores <- prostate_ST$full_first_type.score[prostate_ST$full_first_type == ct]
  
  # Create a safe filename (replace spaces or special characters)
  file_name <- paste0(gsub("[^A-Za-z0-9_]", "_", ct), ".png")
  file_path <- file.path(hist_folder, file_name)
  
  # Open PNG device
  png(filename = file_path, width = 800, height = 600, res = 150)
  
  # Plot histogram
  hist(
    scores,
    main = paste("Confidence Scores for", ct),
    xlab = "full_first_type.score",
    breaks = 50,
    col = "skyblue",
    border = "white"
  )
  
  # Close PNG device
  dev.off()
}


#--------------------------- fixing the labeling -----------------------------
prostate_ST <- readRDS("~/Roselab/Spatial/dietary_project/data/CellTyping_objects/F07833_28_LFRT.rds")
# rm(prostate_ST)


# Extract the current labels
# my_labels <- levels(factor(prostate_ST$first_type))


table(prostate_ST$myccap_strict)

my_labels <- c(
  "MyC-CaP",
  "uncertain_MyC-CaP",
  "T_Cells",
  "NK_Cells",
  "Dendritic_Cells", 
  "ISG-high_Macrophages",
  "Macrophages",
  "Proliferating_Macrophages",
  "Tregs",
  "M2_Macrophages",
  "Mesenchymal", 
  "LumP",
  "Endothelial_Cells", 
  "Fibroblasts", 
  "Unknown"
)

my_colors <- c(
  `MyC-CaP`                  = "black",
  `uncertain_MyC-CaP`        = "grey80",
  T_Cells                    = "firebrick1",
  NK_Cells                   = "darkred",
  Dendritic_Cells            = "indianred3",
  `ISG-high_Macrophages`     = "steelblue3",
  Macrophages                = "blue4",
  Proliferating_Macrophages  = "steelblue2",
  Tregs                      = "violetred3",
  M2_Macrophages             = "deepskyblue4",
  Mesenchymal                = "darkolivegreen4",
  LumP                       = "darkgoldenrod3",
  Endothelial_Cells          = "goldenrod2",
  Fibroblasts                = "darkgoldenrod2",
  Unknown                    = "white"
)

# my_labels_ordered <- my_labels  # or reorder if needed
my_colors <- my_colors[my_labels]  # Ensures correct order and matching names

prostate_ST$myccap_strict <- factor(
  prostate_ST$myccap_strict,
  levels = my_labels
)

fig_directory <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/"

# Generate the plot and assign to variable
p1_strict <- SpatialDimPlot(
  object         = prostate_ST,
  group.by       = "myccap_strict",
  label          = FALSE,
  repel          = TRUE,
  image.alpha    = 1,
  pt.size.factor = 3.5,
  crop           = TRUE,
  cols           = my_colors
)
# print(p1_strict)


# Save as PNG
ggsave(
  filename = file.path(fig_directory, "spatial_celltype3328lfrt.png"),
  plot = p1_strict,
  width = 8,
  height = 6,
  dpi = 600
)

print("stop")

















#----------------- first visualizing attempt -------------------------




# For clarity, re-level them in a desired order
my_labels_ordered <- c(
  "CAFs",
  "CD8+-T-cells",
  "Endothelial-cells",
  "Macrophages",
  "MyC-CaP",       # We'll give this a near-white color
  "Neutrophils",
  "Pericytes",
  "Proliferating-T-cells",
  "uncertain",
  "Unknown"
)

# Ensure that factor levels match your actual labels
my_labels_ordered <- intersect(my_labels_ordered, my_labels) 

# Apply to the metadata
prostate_ST$first_type <- factor(
  prostate_ST$first_type,
  levels = my_labels_ordered
)

# Create a color vector, one color per label in `my_labels_ordered`
my_colors <- c(
  CAFs                = "darkgreen",
  `CD8+-T-cells`      = "tomato",
  `Endothelial-cells` = "orange",
  Macrophages         = "blue",
  `MyC-CaP`           = "black",   # near-white to de-emphasize
  Neutrophils         = "red",
  Pericytes           = "yellow",
  `Proliferating-T-cells` = "deeppink3",
  uncertain           = "pink",
  Unknown             = "grey90"
)[my_labels_ordered]  # subset in correct order








print("stop")

#------------------------------- Garbage? ---------------------------------
# Loading object with RCTD results
# prostate_ST <- readRDS("~/Roselab/Spatial/dietary_project/data/CellTyping_objects/F07834_28_LFRT.rds")

DefaultAssay(prostate_ST) <- "Spatial.008um"
Idents(prostate_ST) <- "full_first_type"
table(Idents(prostate_ST))
fig_directory <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/"





# Generate the plot and assign to variable
p1 <- SpatialDimPlot(
  prostate_ST,
  group.by = "full_first_type",
  label = FALSE,
  repel = TRUE,
  pt.size.factor = 2.5,
  crop = TRUE
)

# Save as PNG
ggsave(
  filename = file.path(fig_directory, "spatial_celltype_map.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

png(filename = file.path(fig_directory, "rctd_assignment_confidence_hist_noline.png"),
    width = 800, height = 600, res = 150)

abline(v = 0.95, col = "red", lwd = 2, lty = 2)

dev.off()







sum(prostate_ST$full_first_type.score < 0.95)

table(prostate_ST$first_type_strict)
table(prostate_ST$full_first_type)

SpatialDimPlot(
  prostate_ST,
  group.by = "first_type_strict",
  label = FALSE,
  repel = TRUE,
  image.alpha = 0.5,
  pt.size.factor = 3,
  crop = TRUE
)

