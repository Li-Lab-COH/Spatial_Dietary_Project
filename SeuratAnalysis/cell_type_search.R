# Load required library
library(Seurat)
library(patchwork)  # For arranging multiple plots

# Define Seurat object
seurat_obj <- F07835_3_RT  # Change if using another dataset
grouping_name <- "BANKSY_snn_res.0.8"


# Function to generate violin plots dynamically
generate_violin_plots <- function(marker_list, cell_type, pointSize = 0.3) {
  plots <- lapply(marker_list, function(gene) {
    VlnPlot(seurat_obj, features = gene, group.by = grouping_name, pt.size = pointSize) + 
      ylim(0.1, NA) + ggtitle(paste(cell_type, "-", gene))
  })
  wrap_plots(plots)  # Arrange plots together
}

#----------------------------- Immune Cell markers-------------------------------
# Natural Killer cells
# nk_markers <- c("Ncr1", "Klrb1c", "Klra8", "klra", "Gzmb", "Ifng", "Prf1")
nk_markers <- c("Gzmb", "Ifng", "Prf1")
generate_violin_plots(nk_markers, "NK Cell")

#t_cell_markers
t_cell_markers <- c("Cd3e", "Cd4", "Cd8a", "Foxp3", "Tbx21", "Gzmk")
generate_violin_plots(t_cell_markers, "T cells")

# tcr markers
tcr_markers <- c("Trac", "Trbc1", "Trbc2", "Trdc")
generate_violin_plots(tcr_markers, "TCR Markers")

# Macrophages
macrophage_markers <- c("Adgre1", "Cd68", "Itgam", "Mrc1")
generate_violin_plots(macrophage_markers, "Macrophage")

# Dendritic Cells (DCs)
dc_markers <- c("Itgax", "H2-Ab1", "Cd74", "Xcr1")
generate_violin_plots(dc_markers, "Dendritic Cell")

# CD4+ T Cells
cd4_markers <- c("Cd3e", "Cd4", "Foxp3", "Il2ra")
generate_violin_plots(cd4_markers, "CD4+ T Cell")

# B Cells 
b_markers <- c("Cd19", "Ms4a1", "Cd79a")
generate_violin_plots(b_markers, "B Cell")

# CD8+ T Cells
cd8_markers <- c("Cd3e", "Cd8a", "Gzmb", "Prf1")
generate_violin_plots(cd8_markers, "CD8+ T Cell")


#------------------------------- Searching MyC-CaP variants ------------------------------
grep("Myc", rownames(F07833_5_RT[["Spatial.008um"]]), value = TRUE, ignore.case = TRUE)
# based on CHI V. DANG* 1999 - mini review of myc genes
# Myc
mycs <- c("Myc", "Mycn", "Mycl")
generate_violin_plots(mycs, "Mycs", 1)

# AR
VlnPlot(seurat_obj, features = "Ar", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

# Cell Growth and Proliferation
myc_growth <- c("Ccna2", "Ccnd1", "Ccne1", "Cdk1", "Cdc25a", "Rcc1", "Odc1", "Tert", "Tk1")
generate_violin_plots(myc_growth, "Cell growth")

# survival
myc_survival <- c("Cdkn2a", "Trp53", "Bax")
generate_violin_plots(myc_survival, "Survival")

# Apoptosis
myc_apoptosis <- c("Bcl2", "Fas")
generate_violin_plots(myc_apoptosis, "Apoptosis")

# Metabolism
myc_metabolism <- c("Ldha", "Cad", "Dhfr", "Eif2s1", "Eif4e", "Rida")
generate_violin_plots(myc_metabolism, "Metabolism")

# Differentiation and Adhesion
myc_differentiation <- c("Cebpa", "Thbs1", "H2-K1", "H2-D1", "Col1a1", "Col1a2", "Col6a3", "Col3a1")
generate_violin_plots(myc_differentiation, "Downregulated differentiation")


#------------------------Searching for markers------------------------------

b_cell_markers <- c("Cd19", "Cd79a", "Ms4a1", "Ighm", "Pax5", "Cd22")
macrophage_markers <- c("Adgre1", "Cd68", "Itgam", "Tnf", "Il1b")
dendritic_cell_markers <- c("Itgax", "Xcr1", "Cd86", "Batf3")
neutrophil_markers <- c("S100a8", "S100a9", "Ly6g", "Cxcr2")

for (gene in tcr_markers) {
  print("")
  print(paste0("current gene search:", gene) )
  found <- grep(gene, rownames(seurat_obj[["Spatial.008um"]]), value = TRUE, ignore.case = TRUE)
  if (length(found) >= 1){
    print(paste0("Genes found: ", paste(found, collapse=", ")))
    # print(found)
  }else{
    print("Nothing found :(")
  }
}



#











  
grep("ccne1", rownames(seurat_obj[["Spatial.008um"]]), value = TRUE, ignore.case = TRUE)


downstream_myc <- c()

VlnPlot(seurat_obj, features = "Ar", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)+
  ggtitle("Ar")



SpatialDimPlot(subset(seurat_obj, idents = c("0")), label = FALSE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 2, alpha = 5)










VlnPlot(seurat_obj, features = "Myc", group.by = grouping_name, pt.size = 0.000000000001) +
  ylim(0.1, NA)  # Set lower bound to 0.5

#--------------------------------CD4 T cells-----------------------------------
#CD4+ T cell marker
VlnPlot(seurat_obj, features = "Cd4", group.by = grouping_name, pt.size = 1) +
  ylim(0.1, NA)  # Set lower bound to 0.5

# Cytotoxic granzyme A
VlnPlot(seurat_obj, features = "Gzma", group.by = grouping_name, pt.size = 1) +
  ylim(0.1, NA)  # Set lower bound to 0.5

# Perforin
VlnPlot(seurat_obj, features = "Prf1", group.by = grouping_name, pt.size = 1) +
  ylim(0.1, NA)  # Set lower bound to 0.5

SpatialDimPlot(subset(seurat_obj, idents = c("3", "4")), label = FALSE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 2, alpha = 5)


#--------------------------Macrophages T cells----------------------------------
plots <- VlnPlot(seurat_obj, features = c("Adgre1", "Cd68", "Itgam"), group.by = grouping_name, pt.size = 0.02, combine = FALSE)

# Apply ylim to each plot
plots <- lapply(plots, function(p) p + ylim(0.5, NA))

# Arrange and display all modified violin plots
patchwork::wrap_plots(plots)



SpatialFeaturePlot(seurat_obj, features = "Adgre1", alpha = c(0.1, 1),
                   image.alpha = 0.3,
                   pt.size.factor = 5) + 
  ggtitle(paste("Spatial Expression of", "Adgre1")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

SpatialFeaturePlot(seurat_obj, features = "Cd68", alpha = c(0.1, 1),
                   image.alpha = 0.3,
                   pt.size.factor = 5) + 
  ggtitle(paste("Spatial Expression of", "Cd68")) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))  

#-------------------------------NK Cells--------------------------------------
VlnPlot(seurat_obj, features = "Cdk4", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5


VlnPlot(seurat_obj, features = "Myc", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5


#-----------------------------MyC-Cap-----------------------------------------
grep("Myc", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
grep("myc", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)


grep("myc|bHLHe39|ENSMUSG00000022346", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

myc_targets <- c("Ccnd1", "Cdk4", "E2f1", "Mki67", "Hk2", "Ldha", "Pkm", "Glut1")
grep(paste(myc_targets, collapse = "|"), rownames(seurat_obj), value = TRUE, ignore.case = TRUE)

VlnPlot(seurat_obj, features = "Ar", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5


VlnPlot(seurat_obj, features = "Myc", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

VlnPlot(seurat_obj, features = "Cdk4", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

VlnPlot(seurat_obj, features = "E2f1", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

VlnPlot(seurat_obj, features = "Mki67", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

VlnPlot(seurat_obj, features = "Hk2", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

VlnPlot(seurat_obj, features = "Ldha", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

VlnPlot(seurat_obj, features = "Pkm", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5

# From online
VlnPlot(seurat_obj, features = "Ck5", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5


Idents(seurat_obj) <- grouping_name

# seurat_obj - search for MyC-CaP
VlnPlot(seurat_obj, features = "Krt8", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5
VlnPlot(seurat_obj, features = "Psca", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5
VlnPlot(seurat_obj, features = "Pkm", group.by = grouping_name, pt.size = 0.3) +
  ylim(0.1, NA)  # Set lower bound to 0.5


# SpatialDimPlot(subset(seurat_obj, idents = c("14", "16"), label = FALSE, repel = TRUE, 
#                pt.size.factor = 5, image.alpha = 0.5, alpha = 5))

# SpatialDimPlot(subset(seurat_obj, idents = c("14", "16"), label = FALSE, repel = TRUE, 
#                       pt.size.factor = 5, image.alpha = 0.05, alpha = 5))

SpatialDimPlot(subset(seurat_obj, idents = c("0")), label = FALSE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 2, alpha = 5)

# SpatialDimPlot(subset(seurat_obj, idents = c("5", "4", "0")), label = FALSE, repel = TRUE, 
#                pt.size.factor = 5, image.alpha = 2, alpha = 5)
