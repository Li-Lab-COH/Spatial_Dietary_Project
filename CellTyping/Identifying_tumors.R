library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
install.packages("enrichR")
library(enrichR)
set.seed(1337)

# Set your working directory (optional)
data_dir <- "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/"

# Load the data
sc_data <- Read10X(data.dir = data_dir)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "Unal_MyC_CaP_WT", min.cells = 3, min.features = 200)

# Calculate % mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)
# Get limits for nFeature_RNA
nf <- seurat_obj@meta.data$nFeature_RNA
nf_low <- quantile(nf, 0.025)
nf_high <- quantile(nf, 0.975)

# Apply all filters together
seurat_obj <- subset(seurat_obj, subset = 
                       nCount_RNA >= 1000 &
                       percent.mt <= 10 &
                       nFeature_RNA >= nf_low &
                       nFeature_RNA <= nf_high
)



# Run standard workflow for clustering
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, dims = 1:10)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

# DotPlot(seurat_obj, features = c("Wfdc12", "Fgb", "HSpb1", "Lcn11", "Cldn3")) + RotatedAxis()



# seurat_obj <- AddModuleScore(seurat_obj, features = list(c("Wfdc12", "Fgb", "Hspb1", "Lcn11", "Cldn3")), name = "tumor_signature")
# VlnPlot(seurat_obj, features = "tumor_signature1", group.by = "seurat_clusters")



# head(FindMarkers(seurat_obj, ident.1 = 12), 20)

tumor_clusters <- c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10)


Idents(seurat_obj) <- "seurat_clusters"
seurat_obj$celltype <- "Other"

seurat_obj$celltype[seurat_obj$seurat_clusters %in% tumor_clusters] <- "MyC-CaP"
Idents(seurat_obj) <- "celltype"

#--------------------------Cell Idents --------------------------------------#
Idents(seurat_obj) <- "seurat_clusters"

# Find markers for all clusters
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 20 markers per cluster
top20_markers_all_clusters <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


tumor_clusters <- c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10)

non_tumor_clusters <- c(6, 11, 12, 13, 14, 15)

top20_non_tumor_markers <- top20_markers_all_clusters %>%
  filter(cluster %in% non_tumor_clusters)

write.csv(top20_non_tumor_markers, "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/top20NonTumor.csv")


#-----------------------subsetting non tumors-------------------------------#

# Define tumor clusters (already done)
tumor_clusters <- c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10)

# Subset to non-tumor cells
non_tumor_obj <- subset(seurat_obj, idents = setdiff(0:15, tumor_clusters))

sum(table(Idents(non_tumor_obj)))

# Normalize and scale
non_tumor_obj <- NormalizeData(non_tumor_obj)
non_tumor_obj <- FindVariableFeatures(non_tumor_obj)
non_tumor_obj <- ScaleData(non_tumor_obj)
non_tumor_obj <- RunPCA(non_tumor_obj)

# Elbow plot (optional)
# ElbowPlot(non_tumor_obj)

# Choose dims and cluster
non_tumor_obj <- FindNeighbors(non_tumor_obj, dims = 1:7)
non_tumor_obj <- FindClusters(non_tumor_obj, dims = 1:7, resolution = 1)  
non_tumor_obj <- RunUMAP(non_tumor_obj, dims = 1:7)

length(unique(Idents(non_tumor_obj)))
DimPlot(non_tumor_obj, reduction = "umap", label = TRUE, pt.size = 0.5)


# Find markers for all clusters
all_markers_nontummor <- FindAllMarkers(non_tumor_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 20 markers per cluster
top20_markers_non_tumors <- all_markers_nontummor %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(top20_markers_non_tumors, "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/top20NonTumor.csv")

VlnPlot(non_tumor_obj, features = c("Xcr1", "Flt3", "Timd4"), group.by = "seurat_clusters", pt.size = 0)
#Tregs
VlnPlot(non_tumor_obj, features = "Foxp3", group.by = "seurat_clusters", pt.size = 0)+
  ylim(0.1,NA)
FeaturePlot(non_tumor_obj, features = "Foxp3")

VlnPlot(non_tumor_obj, features = "Cd8a", group.by = "seurat_clusters", pt.size = 0)+
  ylim(0.1,NA)

#nk cells
FeaturePlot(non_tumor_obj, features = c("Nkg7", "Klrk1", "Klrd1", "Gzmb"))
VlnPlot(non_tumor_obj, features = c("Nkg7", "Klrk1", "Gzmb"), group.by = "seurat_clusters", pt.size = 0)
FeaturePlot(non_tumor_obj, features = c("Nkg7", "Klrk1", "Klrd1"), reduction = "umap")
FeaturePlot(non_tumor_obj, features = c("Gzmb", "Prf1", "Cd244"), reduction = "umap")

#--------------------Final marker list---------------------------------------
# Map cluster numbers to labels
cluster_annotations <- c(
  "2" = "Macrophages",
  "5" = "Macrophages",
  "3" = "Pericytes",
  "7" = "CD8+ T cells",
  "9" = "Proliferating T cells",
  "10" = "Dendritic cells",
  "11" = "Endothelial cells",
  "4" = "Neutrophils",
  "12" = "CAFs",
  "13" = "CAFs",
  "14" = "M2 Macrophages"
)

# Set identities to original clusters if not already
Idents(non_tumor_obj) <- "seurat_clusters"

# Create a new metadata column based on annotations
non_tumor_obj$celltype <- plyr::mapvalues(
  x = as.character(Idents(non_tumor_obj)),
  from = names(cluster_annotations),
  to = cluster_annotations,
  warn_missing = FALSE
)

non_tumor_obj$celltype[non_tumor_obj$celltype %in% c("0", "1", "6", "8")] <- "Unlabeled"

# Set updated identities
Idents(non_tumor_obj) <- "celltype"

# Confirm it worked
table(Idents(non_tumor_obj))

##### labeling Tregs

# Check expression
foxp3_expr <- FetchData(non_tumor_obj, vars = "Foxp3")

# Define Tregs as those with Foxp3 expression above a cutoff (e.g., > 0)
treg_cells <- rownames(foxp3_expr)[foxp3_expr$Foxp3 > 0]

# Update cell type for those cells
non_tumor_obj$celltype[treg_cells] <- "Tregs"

# Set updated identities
Idents(non_tumor_obj) <- "celltype"

# Confirm it worked
table(Idents(non_tumor_obj))


##### B Cells

# Pull expression
b_markers <- FetchData(non_tumor_obj, vars = c("Cd79a", "Cd19"))
FeaturePlot(non_tumor_obj, features = c("Cd79a", "Cd19"))


# Define B cells: expression of at least one B marker > 0
bcell_cells <- rownames(b_markers)[rowSums(b_markers > 0) >= 1]

# Label them
non_tumor_obj$celltype[bcell_cells] <- "B cells"

Idents(non_tumor_obj) <- "celltype"
table(Idents(non_tumor_obj))

# NK cells
# Fetch NK and T cell marker data
nk_markers <- FetchData(non_tumor_obj, vars = c("Nkg7", "Klrk1", "Prf1", "Klrd1", "Klrb1c", "Cd244", "Ifng"))
t_cell_markers <- FetchData(non_tumor_obj, vars = c("Trac", "Trbc1", "Cd3d"))

# Define NK-like cells (express at least 3 NK markers)
nk_like_cells <- rownames(nk_markers)[rowSums(nk_markers > 0.5) >= 3]

# Remove cells that express TCR (any T cell marker > 0.5)
is_t_cell <- rowSums(t_cell_markers > 0.5) > 0
non_t_cells <- rownames(t_cell_markers)[!is_t_cell]

# Final NK cell list: NK-like but not TCR-positive
nk_cells <- intersect(nk_like_cells, non_t_cells)

# Label them
non_tumor_obj$celltype[nk_cells] <- "NK cells"
Idents(non_tumor_obj) <- "celltype"

table(Idents(non_tumor_obj))


#-------------------Transfering labels---------------------------------------

label_df <- data.frame(
  barcode = colnames(non_tumor_obj),
  celltype = Idents(non_tumor_obj)
)

table(Idents(seurat_obj))
table(Idents(non_tumor_obj))

# Ensure that the barcodes and celltype columns in label_df are characters
label_df$barcode <- as.character(label_df$barcode)
label_df$celltype <- as.character(label_df$celltype)

# Find the matching barcodes between seurat_obj and label_df
matching_barcodes <- intersect(colnames(seurat_obj), label_df$barcode)

# (Optional) If seurat_obj$celltype already exists as a factor, convert it to character.
seurat_obj$celltype <- as.character(seurat_obj$celltype)

# Update only the matching barcodes with the new labels from label_df
seurat_obj$celltype[matching_barcodes] <- label_df$celltype[match(matching_barcodes, label_df$barcode)]

# Check the result
table(seurat_obj$celltype, useNA = "ifany")


non_tumor_obj$celltype











#match barcodes
matching_barcodes <- intersect(colnames(seurat_obj), label_df$barcode)

#label according to barcodes
seurat_obj$celltype[matching_barcodes] <- label_df$celltype[match(matching_barcodes, label_df$barcode)]

table(seurat_obj$celltype)



saveRDS(seurat_obj, file = "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/seurat_obj_annotated.rds")

seurat_obj_test <- readRDS("~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/seurat_obj_annotated.rds")


table(seurat_obj_test$celltype)




