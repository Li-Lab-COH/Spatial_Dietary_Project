library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(EBImage)
library(Rfast2)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


#----------------------------Loading Data--------------------------------------
# Preparing Data
# tma1_address <- "../../data/SpacerangerCountOutput/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07837_22WJCYLT3/outs/"
tma1_address <- "../../data/Rose_Li_VisiumHD/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07837_22WJCYLT3/outs"
# tma2_address <- "../../data/SpacerangerCountOutput/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07838_22WJCYLT3/outs/"

tma1_subset <- read.csv("../../data/ExportedSectionLabels/F07837_RL_SKO_TMA1_Slide1.csv")
# tma2_subset <- read.csv("../../data/ExportedSectionLabels/F07838_RL_SKO_TMA2_Slide1.csv")

tma1 <- Load10X_Spatial(tma1_address, bin.size = c(8))
# tma2 <- Load10X_Spatial(tma2_address, bin.size = 8)


#------------------------Subsetting data--------------------------------


# Preparing barcodes
barcodes_31_RT <- subset(tma1_subset, F07837 == "31_RT")
barcodes_31_RT <- as.character(barcodes_31_RT$Barcode)

# barcodes_35_LGI_RT <- subset(tma1_subset, F07837 == "35_LGI_RT")
# barcodes_35_LGI_RT <- as.character(barcodes_35_LGI_RT$Barcode)
# 
# barcodes_33_LGI_RT <- subset(tma1_subset, F07837 == "33_LGI_RT")
# barcodes_33_LGI_RT <- as.character(barcodes_33_LGI_RT$Barcode)
# 
# barcodes_37_CR_RT <- subset(tma1_subset, F07837 == "37_CR_RT")
# barcodes_37_CR_RT <- as.character(barcodes_37_CR_RT$Barcode)


# Subsetting
ms_31RT <- subset(tma1, cells = barcodes_31_RT)

dim(tma1)
dim(ms_31RT)
#------------------------Visualizing--------------------------------

vln.plot <- VlnPlot(ms_31RT, features = "nCount_Spatial.008um", pt.size = 0) + 
  theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(ms_31RT, features = "nCount_Spatial.008um") + 
  theme(legend.position = "right")

vln.plot | count.plot

SpatialFeaturePlot(ms_31RT, features = "nCount_Spatial.008um",
                   min.cutoff = "q10",
                   max.cutoff = "q90",
                   image.alpha = 0.2,
                   pt.size.factor = 5,
                   crop = TRUE)




ms_31RT <- ScaleData(ms_31RT,
                     features = VariableFeatures(ms_31RT),
                     verbose = TRUE)





#------------------------Processing--------------------------------

DefaultAssay(ms_31RT) <- "Spatial.008um"

# perform clustering workflow
ms_31RT <- FindVariableFeatures(ms_31RT)
ms_31RT <- NormalizeData(ms_31RT)
ms_31RT <- ScaleData(ms_31RT)
ms_31RT <- RunPCA(ms_31RT)
ElbowPlot(ms_31RT)
ms_31RT <- FindNeighbors(ms_31RT, dims = 1:15)
ms_31RT <- FindClusters(ms_31RT, resolution = 1)
ms_31RT <- RunUMAP(ms_31RT, reduction = "pca", reduction.name = "umap", 
                   return.model = T, dims = 1:15)
ms_31RT <- RunTSNE(ms_31RT, reduction = "pca", reduction.name = "tsne", 
                   return.model = T, dims = 1:15)


DimPlot(ms_31RT, reduction = "umap", label = F) + 
  ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

SpatialDimPlot(ms_31RT, label = T, repel = T, label.size = 4,
               image.alpha = 0.2,
               pt.size.factor = 5)

FindSpatiallyVariableFeatures(
  tma1 = ms_31RT,  # Your Seurat tma1
  assay = "Spatial.008um",   # Specify the correct assay (check `DefaultAssay(F07837_topright)`)
  selection.method = "moransi",  # OR "markvariogram"
  features = VariableFeatures(ms_31RT),            # Analyze all genes or specify a subset
  nfeatures = 2000,           # Number of top spatially variable genes to select
  image = "slice1.008um",     # Spatial image name (check `F07837_topright@images`)
  verbose = TRUE
)



install.packages('Rfast2')




as.sparse(ms_31RT@assays$Spatial.008um@counts)










