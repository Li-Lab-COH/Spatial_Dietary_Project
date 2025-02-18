library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(EBImage)
library(Rfast2)
library(SeuratWrappers)
library(Banksy)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

BiocManager::install("Banksy")

#----------------------------Loading Data--------------------------------------
# Preparing Data
tma1_address <- "../../data/SpacerangerCountOutput/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07837_22WJCYLT3/outs/"
# tma1_address <- "../../data/Rose_Li_VisiumHD/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07837_22WJCYLT3/outs"

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




# ms_31RT <- ScaleData(ms_31RT,
#                      features = VariableFeatures(ms_31RT),
#                      verbose = TRUE)





#-----------------------------Processing----------------------------------------

DefaultAssay(ms_31RT) <- "Spatial.008um"

# perform clustering workflow
ms_31RT <- NormalizeData(ms_31RT)
ms_31RT <- FindVariableFeatures(ms_31RT)
ms_31RT <- ScaleData(ms_31RT)
ms_31RT <- RunPCA(ms_31RT)
ElbowPlot(ms_31RT)
ms_31RT <- FindNeighbors(ms_31RT, dims = 1:15)
ms_31RT <- FindClusters(ms_31RT, resolution = 1)
ms_31RT <- RunUMAP(ms_31RT, reduction = "pca", reduction.name = "umap", 
                   return.model = T, dims = 1:15)
ms_31RT <- RunTSNE(ms_31RT, reduction = "pca", reduction.name = "tsne", 
                   return.model = T, dims = 1:15)



#################
# Visualizing

p1 <- DimPlot(ms_31RT, reduction = "umap", label = F) + 
  ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p2 <- SpatialDimPlot(ms_31RT, label = T, repel = T, label.size = 4,
               image.alpha = 0.2,
               pt.size.factor = 5)

p1 | p2

SpatialDimPlot(ms_31RT, label = T, repel = T,#, interactive = TRUE, 
               pt.size.factor = 5,
               image.alpha = 0.5,
               alpha = 5)

# Figuring out cluster of interest

DimPlot(subset(ms_31RT, idents = c("13", "9")), reduction = "umap", label = FALSE)

SpatialDimPlot(subset(ms_31RT, idents = c("13", "9")), label = TRUE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 0.5, alpha = 5)

#-------------------------------Gene Expression-------------------------------
#Trying this without subsetting

ms_31RT <- BuildClusterTree(ms_31RT, assay = "Spatial.008um",
                            reduction = "pca", reorder = T)

markers <- FindAllMarkers(ms_31RT, assay = "Spatial.008um", only.pos = TRUE)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

ms_31RT <- ScaleData(ms_31RT, assay = "Spatial.008um", features = top5$gene)

p <- DoHeatmap(ms_31RT, assay = "Spatial.008um", features = top5$gene, size = 2.5) + 
  theme(axis.text = element_text(size = 5.5)) + 
  NoLegend()

p

#TODO: Export table of top expressed genes in each cluster

#-------------------------Spatial tissues-------------------------------

DefaultAssay(ms_31RT) <- "Spatial.008um"
ms_31RT <- RunBanksy(ms_31RT,
                    lambda = 0.2, verbose = TRUE,
                    assay = "Spatial.008um", slot = "data", features = "variable",
                    k_geom = 30
)

DefaultAssay(ms_31RT) <- "BANKSY"
ms_31RT <- RunPCA(ms_31RT, assay = "BANKSY", reduction.name = "pca.banksy", 
                  features = rownames(ms_31RT), npcs = 30)

ms_31RT <- FindNeighbors(ms_31RT, reduction = "pca.banksy", dims = 1:30)
ms_31RT <- FindClusters(ms_31RT, cluster.name = "banksy_cluster", resolution = 0.5)
Idents(ms_31RT) <- "banksy_cluster"
SpatialDimPlot(ms_31RT, group.by = "banksy_cluster", label = T, 
                    repel = T, label.size = 4,
                    image.alpha = 0.8,
                    pt.size.factor = 4)

# DimPlot(ms_31RT, reduction = "pca.banksy", group.by = "banksy_cluster", label = TRUE)
# DimPlot(ms_31RT, reduction = "umap", group.by = "banksy_cluster", label = TRUE)


banksy_cells <- CellsByIdentities(ms_31RT)
SpatialDimPlot(ms_31RT, 
               cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], 
               cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
               combine = T,
               image.alpha = 0.8,
               pt.size.factor = 4) + NoLegend()

############################################
save_path <- "../../figures/ms_31RT/SpatialDimPlot_ms_31RT.png"

# Generate the plot and store it in a variable
spatial_plot <- SpatialDimPlot(ms_31RT, 
                               cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], 
                               cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, 
                               combine = T,
                               image.alpha = 0.8,
                               pt.size.factor = 4) + NoLegend()

# Save as high-quality PNG (300 DPI for publication)
ggsave(filename = save_path, plot = spatial_plot, width = 8, height = 6, dpi = 300)


save_path2 <- "../../figures/ms_31RT/Clusters_4n2n3n14.png"
spatial_4n14 <- SpatialDimPlot(subset(ms_31RT, idents = c("4", "14", "2", "3")), label = FALSE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 2, alpha = 5)
# Save as high-quality PNG (300 DPI for publication)
ggsave(filename = save_path2, plot = spatial_4n14, width = 8, height = 6, dpi = 300)


#########################
# More on cluster 4

# Extract top 10 genes
top10_genes <- c("Troap.m0", "Aspm.m0", "Hmmr.m0", "Kif18b.m0", "Prr11.m0", 
                 "Cenpf.m0", "Aurkb.m0", "Ube2c.m0", "Cdc20.m0", "Cenpe.m0")

# Generate violin plot
VlnPlot(ms_31RT, features = top10_genes, group.by = "seurat_clusters", pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#---------------------------Marker Genes ------------------------------

markers <- FindAllMarkers(
  object = ms_31RT, 
  assay = "BANKSY", 
  only.pos = TRUE, # Only return upregulated genes
  min.pct = 0.25,  # Gene must be expressed in at least 25% of cells
  logfc.threshold = 0.25 # Genes must have at least a 0.25 log-fold change
)

top_markers <- markers %>%
  group_by(cluster) %>%      # Group by cluster
  slice_max(n = 15, order_by = avg_log2FC) %>%  # Select top 5 genes per cluster by log fold change
  ungroup()
subset(top_markers, cluster == 4 )
print(top_markers)

View(top_markers)
VlnPlot(ms_31RT, features = c("Ppp1r1b.m0", "Gene2", "Gene3"), group.by = "banksy_cluster")

dim(markers)
summary(markers$p_val)
# View the top marker genes

# TODO: Need to customize the file output here
write.csv(markers, "../../data/ClusterGenes/ms_31RT_top15.csv", row.names = FALSE)






