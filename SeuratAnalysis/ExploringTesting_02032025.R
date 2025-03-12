library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(EBImage)
library(rstudioapi)


setwd(dirname(getActiveDocumentContext()$path))

F07837_loc = "../../data/SpacerangerCountOutput/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07837_22WJCYLT3/outs/"

data_F07837 = Load10X_Spatial(F07837_loc, bin.size = c(8))



#------------------------First exploration -----------------------------
Assays(data_F07837)


vln.plot <- VlnPlot(data_F07837, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(data_F07837, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot


SpatialFeaturePlot(
  data_F07837, 
  features = "nCount_Spatial.008um", 
  min.cutoff = "q10",  
  max.cutoff = "q90"   
) + theme(legend.position = "right")


SpatialFeaturePlot(
  data_F07837, 
  features = "nFeature_Spatial.008um", 
  min.cutoff = "q10",  
  max.cutoff = "q90"
) + theme(legend.position = "right")

#---------------------------Testing subsetting------------------------------

selected_barcodes <- read.csv("../../data/ExportedPractice/f07837.csv")
selected_barcodes <- subset(selected_barcodes, f07837 == "bottomLeft")
head(selected_barcodes)
barcodes_characters <- as.character(selected_barcodes$Barcode)

#subsetting

BottomLeftData <- subset(data_F07837, cells = barcodes_characters)
dim(selected_barcodes)
dim(BottomLeftData)
dim(data_F07837)

SpatialFeaturePlot(BottomLeftData, features = "nCount_Spatial.008um",
                   min.cutoff = "q10",  
                   max.cutoff = "q90",
                   image.alpha = 0.2,
                   pt.size.factor = 5,
                   crop = TRUE)



Images(BottomLeftData)
dim(BottomLeftData)


image_array <- BottomLeftData@images$slice1.008um@image

image_array[,,3]


#------------------------Normalization--------------------------------------
# normalize both 8um and 16um bins
DefaultAssay(BottomLeftData) <- "Spatial.008um"
BottomLeftData <- NormalizeData(BottomLeftData)


#-----------------------------Exploring Imaging------------------------------

image_array <- BottomLeftData@images$slice1.008um@image
display(image_array)

# Extract the H&E image
hne_image <- BottomLeftData@images$slice1.008um@image
hne_image <- hne_image / max(hne_image)

# Convert and display
library(EBImage)
display(hne_image)

dim(hne_image)

str(BottomLeftData@images$slice1.008um)


# Extract RGB channels
r_channel <- hne_image[,,1]
g_channel <- hne_image[,,2]
b_channel <- hne_image[,,3]

# Combine into an EBImage object
hne_rgb <- array(c(r_channel, g_channel, b_channel), dim = c(dim(hne_image)[1:2], 3))

# Normalize
hne_rgb <- hne_rgb / max(hne_rgb)

# Convert to Image object
hne_rgb <- Image(hne_rgb, colormode = "Color")

# Display
display(hne_rgb)


#--------------------------Visualizing gene expression----------------------

# Androgen Receptor
SpatialFeaturePlot(BottomLeftData, features = "Myc",
                   image.alpha = 0.3,
                   pt.size.factor = 5,
                   label = TRUE) + ggtitle("Androgen receptor 8µm")



#-------------------------Finding Spatially Variable Features------------------

# normalize both 8um and 16um bins
DefaultAssay(BottomLeftData) <- "Spatial.008um"
BottomLeftData <- NormalizeData(BottomLeftData)

BottomLeftData <- ScaleData(BottomLeftData,
                            features = VariableFeatures(BottomLeftData),
                            verbose = TRUE)

FindSpatiallyVariableFeatures(
  object = BottomLeftData,  # Your Seurat object
  assay = "Spatial.008um",   # Specify the correct assay (check `DefaultAssay(F07837_topright)`)
  selection.method = "moransi",  # OR "markvariogram"
  features = NULL,            # Analyze all genes or specify a subset
  nfeatures = 2000,           # Number of top spatially variable genes to select
  image = "slice1.008um",     # Spatial image name (check `F07837_topright@images`)
  verbose = TRUE
)















