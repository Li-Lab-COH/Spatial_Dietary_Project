library(Seurat)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

source("../CustomFunctions/Annotating_subclusters.R")

prostate_ST <- readRDS("~/1Work/RoseLab/Spatial/dietary_droject/data/RCTD_annotated_n_PCA_full/F07835_27_LFRTRCTD_annotated.rds")



prostate_ST <- FindNeighbors(prostate_ST, dims = 1:15, k.param = 20, prune.SNN = 1/15,
                             reduction = "pca.prostate.full")

prostate_ST <- FindClusters(prostate_ST, 
                            resolution = 0.8)

prostate_ST <- RunUMAP(prostate_ST,
                       dims = 1:15,
                       reduction = "pca.prostate.full",
                       verbose = FALSE
                       )


Idents(prostate_ST) <- "first_type"

unique(Idents(prostate_ST))


#-------------------------- level 2 idents ----------------------------------

prostate_ST$level_3 <- prostate_ST$first_type
prostate_ST$first_type <- NULL

Idents(prostate_ST) <- "level_3"

current_idents_level_3 <- Idents(prostate_ST)
head(current_idents_level_3)

# Create a named vector for mapping
level_2_map <- c(
  "MyC_CaP" = "MyC_CaP",
  "Endothelial_Cells" = "Epithelial",
  "Mesenchymal" = "Epithelial",
  "LumP" = "Epithelial",
  "Fibroblasts" = "Epithelial",
  "ISG_high_Macrophages" = "Macrophages",
  "M2_Macrophages" = "Macrophages",
  "Proliferating_Macrophages" = "Macrophages",
  "Macrophages" = "Macrophages",
  "Resident_Macrophages" = "Macrophages",
  "Dendritic_Cells" = "Dendritic_Cells",
  "NK_Cells" = "NK_Cells",
  "Tregs" = "T_Cells",
  "T_Cells" = "T_Cells"
)


# 3. Create a named vector with cell barcodes as names

level_2_vector <- level_2_map[as.character(current_idents_level_3)]
names(level_2_vector) <- names(current_idents_level_3)

prostate_ST$level_2 <- level_2_vector


#-------------------------- level 1 idents ----------------------------------
Idents(prostate_ST) <- "level_2"

current_idents_level_2 <- prostate_ST$level_2


level_1_map <- c(
  "MyC_CaP" = "MyC_CaP",
  "Epithelial" = "Epithelial",
  "Macrophages" = "Immune_Cells",
  "Dendritic_Cells" = "Immune_Cells",
  "NK_Cells" = "Immune_Cells",
  "T_Cells" = "Immune_Cells"
)


level_1_vector <- level_1_map[as.character(current_idents_level_2)]
names(level_1_vector) <- names(current_idents_level_2)

prostate_ST$level_1 <- level_1_vector



#----------------------- level 2 plots -------------------------------------

Idents(prostate_ST) <- "level_2"

custom_colors <- c(
  "MyC_CaP" = "grey93",
  "Epithelial" = "#1f78b4",       # blue
  "Macrophages" = "#33a02c",     # green
  "Dendritic_Cells" = "#ff7f00", # orange
  "NK_Cells" = "#20b2aa",        # purple
  "T_Cells" = "#e31a1c"          # red
)

DimPlot(
  prostate_ST,
  reduction = "umap",
  # label = TRUE,
  cols = custom_colors
)

