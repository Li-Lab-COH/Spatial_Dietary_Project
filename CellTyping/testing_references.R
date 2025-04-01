################################################################################
# Title: Per sample analysis: Initial set up (Single parameter test)
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)       # for ImmGenData
library(SummarizedExperiment)
library(data.table)
library(Matrix)

# Set seed for reproducibility
set.seed(1337)

# Function to create directory if it does not exist
ensure_dir <- function(dir_path) {
  dir_path <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

#----------------------- Preparing Data and Folders -----------------------
file <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/data/Raw_objects/F07835_27_LFRT.rds"
sample_name <- gsub("\\.rds$", "", basename(file))
fig_loc <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/figures/per_sample_analysis/"
fig_loc_sample <- file.path(fig_loc, sample_name)
ensure_dir(fig_loc_sample)

#----------------------- Read and Preprocess Spatial Data -----------------------
prostate_ST <- readRDS(file)
DefaultAssay(prostate_ST) <- "Spatial.008um"

# Quality filtering
prostate_ST <- subset(
  prostate_ST,
  subset = nCount_Spatial.008um >= 100 &
    nFeature_Spatial.008um >= 100 &
    percent.mt <= 15
)

# Normalization and scaling
prostate_ST <- NormalizeData(prostate_ST)
prostate_ST <- FindVariableFeatures(prostate_ST)
prostate_ST <- ScaleData(prostate_ST)

# Run PCA using the Spatial assay; note the custom reduction name
prostate_ST <- RunPCA(
  object = prostate_ST,
  assay = "Spatial.008um",
  reduction.name = "pca.prostate.full",
  verbose = TRUE
)

#----------------------- Set Clustering Parameters -----------------------
k <- 20
res <- 0.2

# Create folder for this parameter combination
comb_folder <- file.path(fig_loc_sample, paste0("neighbors_", k, "_resolution_", res))
ensure_dir(comb_folder)

# Work on a fresh copy of the object
st_obj <- prostate_ST

# Run the neighbor search and clustering using the specified PCA reduction
st_obj <- FindNeighbors(
  st_obj,
  reduction = "pca.prostate.full",  # use the custom reduction name
  dims = 1:15,
  k.param = k,
  prune.SNN = 1/15,
  verbose = FALSE
)

st_obj <- FindClusters(
  st_obj,
  resolution = res,
  verbose = FALSE
)

# Run UMAP using the same PCA reduction
st_obj <- RunUMAP(
  st_obj,
  dims = 1:15,
  reduction = "pca.prostate.full",
  verbose = FALSE
)

#----------------------- UMAP Plot -----------------------
p_umap <- DimPlot(
  st_obj,
  reduction = "umap",
  label = TRUE
) + ggtitle(paste("UMAP: k =", k, "resolution =", res))

ggsave(
  filename = file.path(comb_folder, paste0("UMAP_clusters_k", k, "_res", res, ".png")),
  plot = p_umap,
  width = 8,
  height = 6
)

#----------------------- Violin Plot -----------------------
# Make sure identities are set to clusters
Idents(st_obj) <- "seurat_clusters"

p_vln <- VlnPlot(
  st_obj,
  features = "nCount_Spatial.008um",
  group.by = "seurat_clusters",
  pt.size = 0
) + ggtitle(paste("Violin plot (nCount) k =", k, "resolution =", res))

ggsave(
  filename = file.path(comb_folder, paste0("VlnPlot_nCount_Spatial_k", k, "_res", res, ".png")),
  plot = p_vln,
  width = 8,
  height = 6
)

#----------------------- Marker Identification -----------------------
markers <- FindAllMarkers(
  st_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top20_markers_all <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(
  top20_markers_all,
  file = file.path(comb_folder, paste0("top20_markers_k", k, "_res", res, ".csv")),
  row.names = FALSE
)

#----------------------- SingleR Annotation -----------------------
# Extract normalized expression matrix and cluster identities
expr_mat <- GetAssayData(st_obj, slot = "data")
clusters <- Idents(st_obj)

# (A) SingleR Annotation using ImmGen reference
# Run SingleR using ImmGen reference (immune cells)
ref_immgen <- ImmGenData()
preds_immgen <- SingleR(
  test = expr_mat,
  ref = ref_immgen,
  labels = ref_immgen$label.fine,
  clusters = clusters
)

# Extract scores and add rownames
scores <- preds_immgen$scores
rownames(scores) <- rownames(preds_immgen)

# Top 7 prediction summary
top_matches <- lapply(rownames(scores), function(cluster_id) {
  x <- scores[cluster_id, ]
  sorted <- sort(x, decreasing = TRUE)
  data.frame(
    Top1 = names(sorted)[1], Score1 = sorted[1],
    Top2 = names(sorted)[2], Score2 = sorted[2],
    Top3 = names(sorted)[3], Score3 = sorted[3],
    Top4 = names(sorted)[4], Score4 = sorted[4],
    Top5 = names(sorted)[5], Score5 = sorted[5],
    Top6 = names(sorted)[6], Score6 = sorted[6],
    Top7 = names(sorted)[7], Score7 = sorted[7],
    row.names = cluster_id,
    stringsAsFactors = FALSE
  )
})

top_matches_df <- do.call(rbind, top_matches)
top_matches_df$Cluster <- rownames(top_matches_df)

top_matches_df <- top_matches_df[, c("Cluster",
                                     "Top1", "Score1", "Top2", "Score2", "Top3", "Score3",
                                     "Top4", "Score4", "Top5", "Score5", "Top6", "Score6",
                                     "Top7", "Score7")]

# Save
write.csv(
  top_matches_df,
  file = file.path(comb_folder, paste0("singleR_top7_ImmGen_k", k, "_res", res, ".csv")),
  row.names = FALSE
)


# (B) SingleR Annotation using Karhtaus reference
ref_kar <- readRDS("C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Karthaus_reference_cleaned.rds")

# Extract log-normalized expression data from counts slot
ref_expr <- GetAssayData(ref_kar, slot = "counts")

# Pull cell type labels
ref_labels <- ref_kar$celltype

# Keep only shared genes between reference and test
shared_genes <- intersect(rownames(expr_mat), rownames(ref_expr))
expr_mat_sub <- expr_mat[shared_genes, , drop = FALSE]
ref_expr_sub <- ref_expr[shared_genes, , drop = FALSE]

# Run SingleR
preds_kar <- SingleR(
  test = expr_mat_sub,
  ref = ref_expr_sub,
  labels = ref_labels,
  clusters = clusters
)

# Extract and label scores
scores <- preds_kar$scores
rownames(scores) <- rownames(preds_kar)  # cluster IDs

# Extract top 7 matches per cluster
top_matches <- lapply(rownames(scores), function(cluster_id) {
  x <- scores[cluster_id, ]
  sorted <- sort(x, decreasing = TRUE)
  data.frame(
    Top1 = names(sorted)[1], Score1 = sorted[1],
    Top2 = names(sorted)[2], Score2 = sorted[2],
    Top3 = names(sorted)[3], Score3 = sorted[3],
    Top4 = names(sorted)[4], Score4 = sorted[4],
    Top5 = names(sorted)[5], Score5 = sorted[5],
    Top6 = names(sorted)[6], Score6 = sorted[6],
    Top7 = names(sorted)[7], Score7 = sorted[7],
    row.names = cluster_id,
    stringsAsFactors = FALSE
  )
})


top_matches_df <- do.call(rbind, top_matches)
top_matches_df$Cluster <- rownames(top_matches_df)

top_matches_df <- top_matches_df[, c("Cluster",
                                     "Top1", "Score1", "Top2", "Score2", "Top3", "Score3",
                                     "Top4", "Score4", "Top5", "Score5", "Top6", "Score6",
                                     "Top7", "Score7")]

# Save to CSV
write.csv(
  top_matches_df,
  file = file.path(comb_folder, paste0("singleR_top7_Karthaus_k", k, "_res", res, ".csv")),
  row.names = FALSE
)


# (C) SingleR Annotation using Unal reference
# Load Unal reference (MyC-CaP + others)
ref_unal <- readRDS("C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_condensed.rds")


ref_expr <- GetAssayData(ref_unal, slot = "data")  # this is the normalized data
ref_labels <- ref_unal$celltype

# Match genes between reference and test
shared_genes <- intersect(rownames(expr_mat), rownames(ref_expr))
expr_mat_sub <- expr_mat[shared_genes, , drop = FALSE]
ref_expr_sub <- ref_expr[shared_genes, , drop = FALSE]

# Run SingleR
preds_unal <- SingleR(
  test = expr_mat_sub,
  ref = ref_expr_sub,
  labels = ref_labels,
  clusters = clusters
)

# Extract score matrix and set cluster names
scores <- preds_unal$scores
rownames(scores) <- rownames(preds_unal)

# Build Top 7 predictions table
top_matches <- lapply(rownames(scores), function(cluster_id) {
  x <- scores[cluster_id, ]
  sorted <- sort(x, decreasing = TRUE)
  data.frame(
    Top1 = names(sorted)[1], Score1 = sorted[1],
    Top2 = names(sorted)[2], Score2 = sorted[2],
    Top3 = names(sorted)[3], Score3 = sorted[3],
    Top4 = names(sorted)[4], Score4 = sorted[4],
    Top5 = names(sorted)[5], Score5 = sorted[5],
    Top6 = names(sorted)[6], Score6 = sorted[6],
    Top7 = names(sorted)[7], Score7 = sorted[7],
    row.names = cluster_id,
    stringsAsFactors = FALSE
  )
})

top_matches_df <- do.call(rbind, top_matches)
top_matches_df$Cluster <- rownames(top_matches_df)

top_matches_df <- top_matches_df[, c("Cluster",
                                     "Top1", "Score1", "Top2", "Score2", "Top3", "Score3",
                                     "Top4", "Score4", "Top5", "Score5", "Top6", "Score6",
                                     "Top7", "Score7")]

# Save to CSV
write.csv(
  top_matches_df,
  file = file.path(comb_folder, paste0("singleR_top7_Unal_k", k, "_res", res, ".csv")),
  row.names = FALSE
)

message("Completed test run for k = ", k, " and resolution = ", res)





# Example: the folder for a specific clustering result
comb_folder <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/figures/per_sample_analysis/F07835_27_LFRT/20_0.2_solo_run/"

# List all CSVs you want to combine
csv_files <- list.files(comb_folder, pattern = "\\.csv$", full.names = TRUE)

# Create workbook
wb <- createWorkbook()

# Loop through each file and add as a sheet
for (file_path in csv_files) {
  sheet_name <- tools::file_path_sans_ext(basename(file_path))
  sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet names must be ≤ 31 characters
  
  data <- read.csv(file_path)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}

# Define Excel output file path
excel_output_path <- file.path(comb_folder, "combined_output.xlsx")
saveWorkbook(wb, file = excel_output_path, overwrite = TRUE)






