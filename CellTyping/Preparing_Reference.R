library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)
library(enrichR)
set.seed(1337)


#------------------------------- prepping data -------------------------------
# data_dir <- "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/"
data_dir <- "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/countMatrix/"
fig_directory <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/Unal_reference/"
# out_dir <- "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/"
# Load the data
sc_data <- Read10X(data.dir = data_dir)
# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "Unal_MyC_CaP_WT", min.cells = 3, min.features = 200)


#------------------------- Working with Seurat object -------------------------
# Calculate % mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
#         ncol = 3, pt.size = 0)
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

# fig_directory <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/Unal_reference/"

# Plot clusters
# p1_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
# DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
# 
# ggsave(
#   filename = file.path(fig_directory, "uMap.png"),
#   plot = p1_umap,
#   width = 8,
#   height = 6,
#   dpi = 300
# )
# DotPlot(seurat_obj, features = c("Wfdc12", "Fgb", "HSpb1", "Lcn11", "Cldn3")) + RotatedAxis()



# seurat_obj <- AddModuleScore(seurat_obj, features = list(c("Wfdc12", "Fgb", "Hspb1", "Lcn11", "Cldn3")), name = "tumor_signature")
# VlnPlot(seurat_obj, features = "tumor_signature1", group.by = "seurat_clusters")



# head(FindMarkers(seurat_obj, ident.1 = 12), 20)

tumor_clusters <- c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10)


Idents(seurat_obj) <- "seurat_clusters"
seurat_obj$celltype <- "Other"

seurat_obj$celltype[seurat_obj$seurat_clusters %in% tumor_clusters] <- "MyC-CaP"
Idents(seurat_obj) <- "celltype"

table(Idents(seurat_obj))

#--------------- Finding the 20 markers for MyC-CaP cells ----------------------

Idents(seurat_obj) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
  )


top20_markers_all <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


tumor_clusters <- c(0, 1, 2, 3, 4, 5, 7, 8, 9, 10)

top20_tumors <- top20_markers_all %>%
  filter(cluster %in% tumor_clusters)

output_file <- "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/top20_tumor_markers.csv"
write.csv(top20_tumors, file = output_file, row.names = FALSE)





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
# write.csv(top20_non_tumor_markers, "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/top20NonTumor.csv")



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
View(top20_markers_non_tumors)

# write.csv(top20_markers_non_tumors, "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/top20NonTumor.csv")




#------------------------------ SingleR Non-Tumor------------------------------
ref <- ImmGenData()

# Extract normalized gene expression matrix from your Seurat object
expr_mat <- GetAssayData(non_tumor_obj, slot = "data")
table(non_tumor_obj$seurat_clusters)
# Get cluster identities
clusters <- Idents(non_tumor_obj)

# Run SingleR to annotate each cluster using ImmGenData
singleR_preds <- SingleR(test = expr_mat,
                         ref = ref,
                         labels = ref$label.fine,
                         clusters = clusters)

# See predicted cell type for each cluster
singleR_preds$labels

# View in table format: cluster ID vs predicted label
table(singleR_preds$labels)
View(singleR_preds$scores)

rownames(singleR_preds)
# Get top 3 labels for each cluster
# This doesn't take in the cluster label lapply does, here it uses the first 
# column as row name
# top_matches <- apply(singleR_preds$scores, 1, function(x) {
#   sorted <- sort(x, decreasing = TRUE)
#   data.frame(
#     Top1 = names(sorted)[1],
#     Score1 = sorted[1],
#     Top2 = names(sorted)[2],
#     Score2 = sorted[2],
#     Top3 = names(sorted)[3],
#     Score3 = sorted[3]
#   )
# })
# top_matches_df <- do.call(rbind, top_matches)
# top_matches_df$Cluster <- rownames(top_matches_df)

# top_matches <- apply(singleR_preds, 1, function())

scores <- singleR_preds$scores
rownames(scores) <- rownames(singleR_preds)


top_matches <- lapply(rownames(scores), function(cell_id){
  x <- scores[cell_id, ]
  sorted <- sort(x, decreasing = TRUE)
  data.frame(
    Top1 = names(sorted)[1], 
    Score1 = sorted[1],
    Top2 = names(sorted)[2],
    Score2 = sorted[2],
    Top3 = names(sorted)[3],
    Score3 = sorted[3],
    row.names = cell_id,
    stringsAsFactors = FALSE
  )
})

top_matches_df <- do.call(rbind, top_matches)
top_matches_df$Cluster <- rownames(top_matches_df)

# Reorder columns
top_matches_df <- top_matches_df[, c("Cluster", "Top1", "Score1", "Top2", "Score2", "Top3", "Score3")]

# View table
print(top_matches_df)

# non_tumor_df <- file.path(out_dir, "Non_tumor_matches.csv")
# write.csv(top_matches_df, non_tumor_df)
# 
# out_data_dir <- "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/non_tumor_cells.rds"
# saveRDS(non_tumor_obj, out_data_dir)

#----------------------- Transfering SingleR lables ---------------------------

cluster_annotations <- c(
  "0" = "ISG-high Macrophages",
  "1" = "M2 Macrophages",
  "2" = "Macrophages",
  "3" = "Mesenchymal",
  "4" = "LumP",
  "5" = "Macrophages",
  "6" = "Macrophages",
  "7" = "NK + T Cells",
  "8" = "Ambiguous/Low-Quality",
  "9" = "Proliferating Macrophages",
  "10" = "Dendritic Cells",
  "11" = "Endothelial Cells",
  "12" = "Fibroblasts",
  "13" = "Fibroblasts",
  "14" = "Resident Macrophages"
)
Idents(non_tumor_obj) <- "seurat_clusters"

non_tumor_obj$celltype <- plyr::mapvalues(
  x = as.character(Idents(non_tumor_obj)),
  from = names(cluster_annotations),
  to = cluster_annotations,
  warn_missing = FALSE
)

Idents(non_tumor_obj) <- "celltype"

table(Idents(non_tumor_obj))

other_cells_fig_dir <- "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/UMAP_non_tumors.png"
non_tumor_imgage <- DimPlot(non_tumor_obj, group.by = "celltype", label = TRUE, repel = TRUE)

ggsave(
  filename = other_cells_fig_dir,
  plot = non_tumor_imgage,
  width = 8,
  height = 6,
  dpi = 300
)

# Testing macrophage markers

# Define canonical macrophage markers
macrophage_markers <- c("Lyz2", "Cd68", "C1qa", "Adgre1")

# Violin plot across all clusters
VlnPlot(non_tumor_obj, features = macrophage_markers, group.by = "seurat_clusters", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Investigating cluster 4
VlnPlot(non_tumor_obj, features = c("nCount_RNA", "nFeature_RNA"), group.by = "seurat_clusters")
VlnPlot(non_tumor_obj, features = c("Krt8", "Krt7", "Epcam", "Cldn3", "Rab25")) #epithelial markers
#luminal markers
VlnPlot(non_tumor_obj, features = c("Ppp1r1b", "Msmb", "Cldn10", "Ly6d", "Aqp3"), group.by = "seurat_clusters")

#
# ----------------------------- Subsetting Nk + T cells ----------------------
Idents(non_tumor_obj) <- "seurat_clusters"
# Subset cluster 7 from the main object
cluster7_obj <- subset(non_tumor_obj, idents = "7")
cluster7_obj <- NormalizeData(cluster7_obj)
cluster7_obj <- FindVariableFeatures(cluster7_obj)
cluster7_obj <- ScaleData(cluster7_obj)


cluster7_obj <- RunPCA(cluster7_obj)
ElbowPlot(cluster7_obj)
cluster7_obj <- RunUMAP(cluster7_obj, dims = 1:20)
cluster7_obj <- FindNeighbors(cluster7_obj, dims = 1:20)
cluster7_obj <- FindClusters(cluster7_obj, resolution = 0.8)  # You can tweak resolution

DimPlot(cluster7_obj, label = TRUE, group.by = "seurat_clusters")

FeaturePlot(cluster7_obj, features = c("Cd3d", "Cd3e", "Cd2"))    # T cells
FeaturePlot(cluster7_obj, features = c("Nkg7", "Klrk1", "Gzmb"))  # NK cells

# Assign new labels based on cluster ID
Idents(cluster7_obj) <- "seurat_clusters"
cluster7_obj$celltype <- "Unknown"
cluster7_obj$celltype[Idents(cluster7_obj) == "1"] <- "T Cells"
cluster7_obj$celltype[Idents(cluster7_obj) == "0"] <- "NK Cells"

Idents(cluster7_obj) <- "celltype"
table(Idents(cluster7_obj))

# Transfering labels to non_tumor_obj
t_nk_labels <- data.frame(
  barcode = colnames(cluster7_obj),
  tnk_identity = Idents(cluster7_obj),
  stringsAsFactors = FALSE
)
non_tumor_obj$celltype <- as.character(non_tumor_obj$celltype)

#Assign new labels where barcodes match
non_tumor_obj$celltype[t_nk_labels$barcode] <- as.character(t_nk_labels$tnk_identity)
non_tumor_obj$celltype <- factor(non_tumor_obj$celltype)


table(non_tumor_obj$celltype)

#-------------------------------- Finding Foxp3 --------------------------
# Check expression
foxp3_expr <- FetchData(non_tumor_obj, vars = "Foxp3")

# Define Tregs as those with Foxp3 expression above a cutoff (e.g., > 0)
treg_cells <- rownames(foxp3_expr)[foxp3_expr$Foxp3 > 0.4]
length(treg_cells)
# Convert to character so you can safely assign new label
non_tumor_obj$celltype <- as.character(non_tumor_obj$celltype)

# Now assign "Tregs" to Foxp3+ cells
non_tumor_obj$celltype[treg_cells] <- "Tregs"

# (Optional) Convert back to factor for clean plotting
non_tumor_obj$celltype <- factor(non_tumor_obj$celltype)


# Saving fully annotated non_tumor_cells
Idents(non_tumor_obj) <- "celltype"
table(Idents(non_tumor_obj))

df <- as.data.frame(table(Idents(non_tumor_obj)))
colnames(df) <- c("CellType", "Count")



non_tumor_cells <- ggplot(df, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic() +
  labs(title = "Cell Type Distribution", x = "Cell Type", y = "Number of Cells") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave(
  filename = "~/Roselab/Spatial/dietary_project/figures/Cell_Typing_RCTD/Unal_reference/non_tumor_reference_count.png",
  plot = non_tumor_cells,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
# saveRDS(non_tumor_obj, "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/non_tumor_cells_obj.rds")
#-------------------Transfering All labels---------------------------------------
table(Idents(non_tumor_obj))
table(Idents(seurat_obj))
dim(non_tumor_obj)

# Step 1: Create label dataframe from non_tumor_obj
label_df <- data.frame(
  barcode = colnames(non_tumor_obj),
  celltype = as.character(Idents(non_tumor_obj)),  # Ensures character from the start
  stringsAsFactors = FALSE
)

# Step 2: Match barcodes
matching_barcodes <- intersect(colnames(seurat_obj), label_df$barcode)

# Step 3: Prep seurat_obj$celltype to avoid factor issues
seurat_obj$celltype <- as.character(seurat_obj$celltype)

# Step 4: Transfer labels using match() to align positions
seurat_obj$celltype[matching_barcodes] <- label_df$celltype[
  match(matching_barcodes, label_df$barcode)
]

# Step 5: (Optional) Convert to factor for plotting
seurat_obj$celltype <- factor(seurat_obj$celltype)

# Step 6: Confirm it worked
table(seurat_obj$celltype, useNA = "ifany")

#-----------Removing ambiguous population and saving --------------------------

# Subset the object to keep only cells that are NOT labeled as "Ambiguous/Low-Quality"
tumor_reference_clean <- subset(seurat_obj, subset = celltype != "Ambiguous/Low-Quality")
Idents(tumor_reference_clean) <- "celltype"

table(Idents(tumor_reference_clean))

saveRDS(tumor_reference_clean, "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_unal_reference.rds")






#------------------------------ Older stuff -----------------------------------

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
# cluster_annotations <- c(
#   "2" = "Macrophages",
#   "5" = "Macrophages",
#   "3" = "Pericytes",
#   "7" = "CD8+ T cells",
#   "9" = "Proliferating T cells",
#   "10" = "Dendritic cells",
#   "11" = "Endothelial cells",
#   "4" = "Neutrophils",
#   "12" = "CAFs",
#   "13" = "CAFs",
#   "14" = "M2 Macrophages"
# )

# Reannotating
BiocManager::install("SingleR")
BiocManager::install("celldex")  

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













#match barcodes
matching_barcodes <- intersect(colnames(seurat_obj), label_df$barcode)

#label according to barcodes
seurat_obj$celltype[matching_barcodes] <- label_df$celltype[match(matching_barcodes, label_df$barcode)]

table(seurat_obj$celltype)



saveRDS(seurat_obj, file = "~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/seurat_obj_annotated.rds")

seurat_obj_test <- readRDS("~/1Work/RoseLab/Spatial/Dietary_Project/data/Unal_2024_Myc_CaP/seurat_obj_annotated.rds")


table(seurat_obj_test$celltype)




