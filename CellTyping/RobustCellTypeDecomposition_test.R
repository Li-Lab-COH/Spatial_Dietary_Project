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
prostate_ST <- readRDS("~/Roselab/Spatial/dietary_project/data/Rogelio/Analysis_final/F07834_28_LFRT.rds")
MyC_CaP_ref <- readRDS("~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/MyC_CaP_Reference.rds")




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
  !(MyC_CaP_ref$celltype %in% c("Unlabeled", "B cells", "M2 Macrophages", "NK cells"))
]

MyC_CaP_ref <- subset(MyC_CaP_ref, cells = cells_to_keep)
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

# TODO: have this be the first think you do
# Run PCA command on the entire object
DefaultAssay(prostate_ST) <- "Spatial.008um"
prostate_ST <- RunPCA(
  object = prostate_ST,
  assay = "Spatial.008um",
  reduction.name = "pca.prostate.full",
  verbose = T
)
# TODO: change all dims to 1:15
ElbowPlot(prostate_ST, reduction = "pca.prostate.full", ndims = 50)


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

saveRDS(prostate_ST, "~/Roselab/Spatial/dietary_project/data/CellTyping_objects/F07834_28_LFRT.rds")
#------------------------------- Visualizing ---------------------------------
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
hist(prostate_ST$full_first_type.score, breaks = 100,
     main = "Confidence of RCTD Assignments",
     xlab = "first_type_prob", col = "gray")
abline(v = 0.95, col = "red", lwd = 2, lty = 2)

dev.off()


# After a 0.95 cutoff
prostate_ST$first_type_strict <- ifelse(
  prostate_ST$full_first_type == "MyC-CaP" & prostate_ST$full_first_type.score < 0.95,
  "uncertain",
  prostate_ST$full_first_type
)

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

#--------------------------- fixing the labeling -----------------------------
# prostate_ST <- readRDS("~/Roselab/Spatial/dietary_project/data/CellTyping_objects/F07834_28_LFRT.rds")
rm(prostate_ST)


# Extract the current labels
my_labels <- levels(factor(prostate_ST$first_type_strict))

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
prostate_ST$first_type_strict <- factor(
  prostate_ST$first_type_strict,
  levels = my_labels_ordered
)

# Create a color vector, one color per label in `my_labels_ordered`
my_colors <- c(
  CAFs                = "darkgreen",
  `CD8+-T-cells`      = "tomato",
  `Endothelial-cells` = "orange",
  Macrophages         = "blue",
  `MyC-CaP`           = "grey50",   # near-white to de-emphasize
  Neutrophils         = "red",
  Pericytes           = "yellow",
  `Proliferating-T-cells` = "deeppink3",
  uncertain           = "pink",
  Unknown             = "black"
)[my_labels_ordered]  # subset in correct order

fig_directory <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/"

# Generate the plot and assign to variable
p1_strict <- SpatialDimPlot(
  object         = prostate_ST,
  group.by       = "first_type_strict",
  label          = FALSE,
  repel          = TRUE,
  image.alpha    = 0.9,
  pt.size.factor = 3,
  crop           = TRUE,
  cols           = my_colors
)
print(p1_strict)


# Save as PNG
ggsave(
  filename = file.path(fig_directory, "spatial_celltype_map_strict.png"),
  plot = p1_strict,
  width = 8,
  height = 6,
  dpi = 300
)



