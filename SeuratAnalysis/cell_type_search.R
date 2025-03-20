# Load required library
library(Seurat)
library(patchwork)  # For arranging multiple plots
library(ggplot2)

# Define Seurat object
rm(F07835_27_LFRT)
seurat_obj <- F07834_28_LFRT  # Change if using another dataset
# seurat_obj <- F07833_5_RT
grouping_name <- "BANKSY_snn_res.0.8"
grouping_name <- "BANKSY_snn_res.0.5"

# Function to generate violin plots dynamically
# generate_violin_plots <- function(marker_list, cell_type, pointSize = 0.3) {
#   plots <- lapply(marker_list, function(gene) {
#     VlnPlot(seurat_obj, features = gene, group.by = grouping_name, pt.size = pointSize) + 
#       ylim(0.1, NA) + ggtitle(paste(cell_type, "-", gene))
#   })
#   wrap_plots(plots)  # Arrange plots together
# }
# generate_violin_plots <- function(marker_list, cell_type, pointSize = 0.3) {
#   plots <- lapply(marker_list, function(gene) {
#     VlnPlot(seurat_obj, features = gene, group.by = grouping_name, pt.size = 0.3) +
#       coord_cartesian(ylim = c(1, NA)) +
#       ggtitle(paste(cell_type, "-", gene))
#   })
#   wrap_plots(plots)  # Arrange plots together
# }

generate_violin_plots <- function(marker_list, cell_type, pointSize = 0.3) {
  plots <- lapply(marker_list, function(gene) {
    # Try to fetch data safely
    expr_data <- tryCatch(
      FetchData(seurat_obj, vars = gene),
      error = function(e) {
        message(paste("Warning: Gene", gene, "not found in Seurat object. Skipping..."))
        return(NULL)  # Skip plotting if gene isn't found
      }
    )
    
    if (is.null(expr_data)) return(ggplot() + ggtitle(paste(cell_type, "-", gene, "(Not Found)")))
    
    meta_data <- seurat_obj@meta.data
    plot_data <- data.frame(Expression = expr_data[, 1], Group = meta_data[[grouping_name]])
    
    # Remove zero values
    plot_data <- plot_data[plot_data$Expression > 0, ]
    
    # If all values are zero, return a blank plot
    if (nrow(plot_data) == 0) {
      return(ggplot() + ggtitle(paste(cell_type, "-", gene, "(No Non-Zero Values)")))
    }
    
    # Suppress warnings only for the ggplot call
    suppressWarnings(
      suppressMessages(
        ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
          geom_violin(scale = "width") +
          geom_jitter(width = 0.2, size = pointSize, alpha = 0.5) +
          coord_cartesian(ylim = c(min(plot_data$Expression), NA)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle(paste(cell_type, "-", gene))
      )
    )
  })
  
  wrap_plots(plots)  # Arrange plots together
}




#----------------------------- Immune Cell markers-------------------------------
# Natural Killer cells
# nk_markers <- c("Ncr1", "Klrb1c", "Klra8", "klra", "Gzmb", "Ifng", "Prf1")
nk_markers <- c("Gzmb", "Ifng", "Prf1", "Ncr1", "Itga2")
generate_violin_plots(nk_markers, "NK Cell")

#t_cell_markers - cd4, treg, and cd8
t_cell_markers <- c("Cd3e", "Cd4", "Cd8a", "Foxp3", "Tbx21", "Gzmk", "Il2ra", "Gzmb", "Prf1")
generate_violin_plots(t_cell_markers, "T cells")

# tcr markers
tcr_markers <- c("Trac", "Trbc1", "Trbc2", "Trdc")
generate_violin_plots(tcr_markers, "TCR Markers")

# be cell markers
b_cell_markers <- c("Cd19", "Cd79a", "Ms4a1", "Ighm", "Pax5", "Cd22")
generate_violin_plots(b_cell_markers, "B cells")

# Macrophages
macrophage_markers <- c("Adgre1", "Cd68", "Itgam", "Mrc1", "Tnf", "Il1b", "Csf1r")
generate_violin_plots(macrophage_markers, "Macrophage")

# Dendritic Cells (DCs)
dc_markers <- c("Itgax", "H2-Ab1", "Cd74", "Xcr1", "Cd86", "Batf3")
generate_violin_plots(dc_markers, "Dendritic Cell")


# -------------------------- Other immune cells ------------------------------
# Myeloid-derived suppressor cells
mdsc_markers <- c("Ly6g", "Arg1", "Nos2", "S100a8", "S100a9")
generate_violin_plots(mdsc_markers, "MDSCs")

# Tregs
treg_markers <- c("Foxp3", "Il2ra", "Ctla4", "Tigit", "Ikzf2")
generate_violin_plots(treg_markers, "Regulatory T Cells")

# Exhausted T cells
exhausted_t_cell_markers <- c("Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tox")
generate_violin_plots(exhausted_t_cell_markers, "Exhausted T Cells")

# Tumor associated macrophages (TAM)
# anti-tumor
tam_markers_anti_tumor <- c("Nos2", "Cd86", "Il12b") 
generate_violin_plots(tam_markers_anti_tumor, "Tumor-Associated Macrophages")

# pro tumorigenic TAM
pro_tumor_TAM <- c("Mrc1", "Arg1")
generate_violin_plots(pro_tumor_TAM, "Pro tumorigenic TAM")

# Neutrophils and Tumor-associated neutrophils (TANs)
# N1 Neutrophils - Anti Tumor
n1_markers <- c("Cxcr3", "Il12b", "Nos2", "Tnf", "Il1b", "Ccl3", "Ccl4", "Ccl5", "Mpo")
generate_violin_plots(n1_markers, "N1 Neutrophils")

# N2 Neutrophils (Pro-Tumor)
n2_markers <- c("S100a8", "S100a9", "Arg1", "Tgfb1", "Vegfa")
generate_violin_plots(n2_markers, "N2 Neutrophils")

# Mast cells
mast_cell_markers <- c("Kit", "Fcer1a", "Mcpt1", "Cpa3")
generate_violin_plots(mast_cell_markers, "Mast Cells")

# Fibroblast
fibroblast_markers <- c("Fap", "Pdgfrb", "Col1a1", "Acta2")
generate_violin_plots(fibroblast_markers, "Fibroblast-Associated Immune Cells")

#--------------------------- Searching MyC-CaP Effects--------------------------
# based on CHI V. DANG* 1999 - mini review of myc genes
# Myc
mycs <- c("Myc", "Mycn", "Mycl", "Ar")
generate_violin_plots(mycs, "Mycs", 0.3)

# AR
# VlnPlot(seurat_obj, features = "Ar", group.by = grouping_name, pt.size = 0.3) +
#   ylim(0.1, NA)  # Set lower bound to 0.5

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

# -------------------Additional for SKO cells ---------------------------------
# Prostate-Specific & Epithelial Markers
myc_prostate_markers <- c("Epcam", "Nkx3-1", "Krt8", "Krt18")
generate_violin_plots(myc_prostate_markers, "Prostate-specific markers")

# Additional Proliferation Markers
myc_proliferation <- c("Mki67", "E2f1", "E2f2", "E2f3")
generate_violin_plots(myc_proliferation, "Proliferation markers")

# Myc-Driven Metabolism and Transporters
myc_metabolic_transporters <- c("Slc7a5")
generate_violin_plots(myc_metabolic_transporters, "Myc metabolic transporters")

# Optional: If your SKO line was derived from a probasin-driven model
myc_probasin <- c("Pbsn")
generate_violin_plots(myc_probasin, "Probasin expression")

#--------------------- Epithelial cell markers---------------------------------
# Based on the Crowley paper - Crowley, Laura, et al. Elife 9 (2020): e59465.

# Lum A, D, L, P
luminal_markers <- c("Tgm4", "Msmb", "Ppp1r1b")
generate_violin_plots(luminal_markers, "Lumin A/D/L/P")

# lumP
lump_markers <- c("Ppp1r1b", "Cldn10")
generate_violin_plots(lump_markers, "Proximal Luminal Cells (LumP)")

# Periurethral epithelial population
pru_markers <- c("Ly6d", "Aqp3", "Ppp1r1b")
generate_violin_plots(pru_markers, "Periurethral Epithelial Cells (PrU)")


#------------------------------- Pro-tumorogenic --------------------------------
pro_tumor <- c("Cd274", "Arg1", "Mmp9")
generate_violin_plots(pro_tumor, "pro_tumor")

N2_markers <- c("Csf3", "Retnlg", "Tgm2", "Cxcr2")
generate_violin_plots(N2_markers, "N2")

M2_markers <- c("Cd163", "Mrc1", "Il10", "Ccl22")
generate_violin_plots(M2_markers, "M2")
#------------------------Searching for markers------------------------------
adgre <- c("Lc3a")

for (gene in adgre) {
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


# Visualizing 
Idents(seurat_obj) <- seurat_obj$BANKSY_snn_res.0.8
SpatialDimPlot(subset(seurat_obj, idents = c("0", "4", "6")), label = FALSE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 2, alpha = 5)


Idents(seurat_obj) <- seurat_obj$BANKSY_snn_res.0.5
SpatialDimPlot(subset(seurat_obj, idents = c("0", "1")), label = FALSE, repel = TRUE, 
               pt.size.factor = 5, image.alpha = 2, alpha = 5)

#--------------------------- Older Stuff --------------------------------------











  
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
