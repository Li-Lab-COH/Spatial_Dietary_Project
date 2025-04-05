library(Seurat)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

source("../CustomFunctions/Annotating_subclusters.R")
outdir <- "../../Results/RnJ/Updating_clustering04.04.25/Figures/"
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

custom_colors_level_2 <- c(
  "MyC_CaP" = "grey93",
  "Epithelial" = "grey50",       # blue
  "Macrophages" = "#33a02c",     # green
  "Dendritic_Cells" = "#ff7f00", # orange
  "NK_Cells" = "#20b2aa",        # purple
  "T_Cells" = "#e31a1c"          # red
)

level_2_umap <- DimPlot(
  prostate_ST,
  reduction = "umap",
  cols = custom_colors_level_2
)

Idents(prostate_ST) <- "seurat_clusters"

cluster_umap <- DimPlot(
  prostate_ST,
  reduction = "umap",
  label = TRUE
)

Idents(prostate_ST) <- "level_1"
# unique(Idents(prostate_ST))

custom_colors_level_1 <- c(
  "MyC_CaP" = "grey93",
  "Epithelial" = "grey50",       # blue
  "Immune_Cells" = "#33a02c"     # green
)


level_1_umap <- DimPlot(
  prostate_ST,
  reduction = "umap",
  cols = custom_colors_level_1
)

level_2_umap | cluster_umap | level_1_umap
combined_umap <- level_2_umap | cluster_umap | level_1_umap

ggsave(
  filename = file.path(outdir, "UMAP_level1_level2_clusters.png"),
  plot = combined_umap,
  width = 15, height = 6, units = "in",
  dpi = 300
)


#---------------------------- Immune plots ------------------------------------
library(Seurat)
library(patchwork)

features_to_plot <- c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Foxp3")
# features_to_plot <- c("Adgre1", "Cd68", "Csf1r", "Cd14", "Cd163")
# For storing both types of plots
marker_plots <- list()
ident_plots <- list()

for (gene in features_to_plot) {
  
  umap_coords <- Embeddings(prostate_ST, reduction = "umap")
  xlims <- range(umap_coords[, 1])
  ylims <- range(umap_coords[, 2])
  
  
  # Get expression for this gene
  expr_vals <- GetAssayData(prostate_ST, assay = "Spatial.008um", layer = "data")[gene, ]
  nonzero_cells <- names(expr_vals[expr_vals > 0])
  
  if (length(nonzero_cells) == 0) {
    message("No non-zero cells for ", gene)
    next
  }
  
  p_marker <- FeaturePlot(
    subset(prostate_ST, cells = nonzero_cells),
    features = gene,
    cols = c("grey95", "blue"),
    order = TRUE,
    min.cutoff = 0
  ) +
    xlim(xlims) +
    ylim(ylims) +
    ggtitle(gene)
  
  p_ident <- DimPlot(
    subset(prostate_ST, cells = nonzero_cells),
    group.by = "level_2",
    label = FALSE,
    pt.size = 1
  ) +
    xlim(xlims) +
    ylim(ylims) +
    ggtitle(paste(gene, "- level_2"))
  
  
  # Store them
  marker_plots[[gene]] <- p_marker
  ident_plots[[gene]] <- p_ident
}

# # Combine into a grid: for each gene, show marker + level_2 plot side-by-side
# combined_plots <- mapply(function(m, i) m | i, marker_plots, ident_plots, SIMPLIFY = FALSE)
# wrap_plots(combined_plots, ncol = 1)

feature_plot_markers <- patchwork::wrap_plots(marker_plots, ncol = 3)
feature_plot_idents <- patchwork::wrap_plots(ident_plots, ncol = 3)




# ggsave(
#   filename = file.path(outdir, "Macrophage_markers.png"),
#   plot = feature_plot_markers,
#   width = 15, height = 15, units = "in",
#   dpi = 300
# )
# 
# ggsave(
#   filename = file.path(outdir, "Macrophage_idents.png"),
#   plot = feature_plot_idents,
#   width = 15, height = 15, units = "in",
#   dpi = 300
# )


ggsave(
  filename = file.path(outdir, "Tcell_markers.png"),
  plot = feature_plot_markers,
  width = 15, height = 15, units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "Tcell_idents.png"),
  plot = feature_plot_idents,
  width = 15, height = 15, units = "in",
  dpi = 300
)








# T cells
FeaturePlot(prostate_ST, features = c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Foxp3"), ncol = 3)

# Core Macrophage marker
FeaturePlot(prostate_ST, features = c("Adgre1", "Cd68", "Csf1r", "Cd14", "Cd163"), ncol = 3)
# FeaturePlot(prostate_ST, features = c("Nos2", "Il1b", "Tnf", "Mrc1", "Arg1", "Retnla"), ncol = 3)
p_m1 <- FeaturePlot(
  prostate_ST,
  features = c("Nos2", "Il1b", "Tnf"),
  ncol = 3,
  cols = c("grey95", "firebrick"),
  order = TRUE,
  min.cutoff = "q5"
)
p_m2 <- FeaturePlot(
  prostate_ST,
  features = c("Mrc1", "Arg1", "Retnla"),
  ncol = 3,
  cols = c("grey95", "forestgreen"),
  order = TRUE,
  min.cutoff = "q5"
)

p_m1 / p_m2

# NK Cells
FeaturePlot(prostate_ST, features = c("Nkg7", "Klrb1c", "Prf1", "Gzma", "Gzmb", "Ifng"), ncol = 3)
# Dendritic Cells
FeaturePlot(prostate_ST, features = c("Itgax", "H2-Aa", "Zbtb46", "Flt3"), ncol = 2)
# Neutrophils
FeaturePlot(prostate_ST, features = c("S100a8", "S100a9", "Ly6g", "Mpo", "Elane"), ncol = 5)


FeaturePlot(
  prostate_ST,
  features = c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Foxp3"),
  ncol = 3,
  min.cutoff = 1,  # or "q5"
  order = TRUE,
  cols = c("grey95", "blue")
)

FeaturePlot(
  prostate_ST,
  features = c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Foxp3"),
  ncol = 3,
  order = TRUE,
  cols = c("lightgrey", "blue"),
  min.cutoff = "q5"
)

FeaturePlot(prostate_ST, features = c("Nos2", "Il1b", "Tnf", "Mrc1", "Arg1", "Retnla"), ncol = 3)


FeaturePlot(prostate_ST, features = c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Foxp3"), ncol = 3)



FeaturePlot(
  subset(prostate_ST, subset = level_1 == "Immune_Cells"),
  features = c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Foxp3"),
  ncol = 3,
  order = TRUE,
  cols = c("lightgrey", "blue"),
  min.cutoff = "q5"
)
