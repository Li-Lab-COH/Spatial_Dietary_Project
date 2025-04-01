library(Seurat)
library(Matrix)
library(data.table)

data_location <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Original_data/SCP859/expression/5e9c0404771a5b0f0fbe2beb/"
metadata_file <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Original_data/SCP859/metadata/mmProstate10x_scPortal_metadata.txt"
save_file <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Karthaus_reference.rds"
# Load expression matrix
counts <- readMM(file.path(data_location, "matrix.mtx.gz"))
genes <- fread(file.path(data_location, "matrix.genes.tsv.gz"), header = FALSE)
barcodes <- fread(file.path(data_location, "batch_barcode.tsv.gz"), header = FALSE)

rownames(counts) <- make.unique(genes$V2)
colnames(counts) <- barcodes$V1

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

# Load metadata
metadata <- fread(metadata_file)
rownames(metadata) <- metadata$NAME
seurat_obj <- AddMetaData(seurat_obj, metadata)

# Add celltype as identity
Idents(seurat_obj) <- seurat_obj$celltype # replace 'celltype' with actual column name

saveRDS(seurat_obj, save_file)
#------------------------- Checking gene name matches --------------------------
# Get gene names from spatial data (Seurat object)
spatial_genes <- rownames(prostate_ST)

# Get gene names from the reference data
ref_genes <- make.unique(genes$V2)  # ref_expr = GetAssayData(reference_obj, slot = "data")

# Check overlap
shared_genes <- intersect(spatial_genes, ref_genes)
length(shared_genes)  # 18k genes

# How many spatial genes are not in the reference?
length(setdiff(spatial_genes, ref_genes))# 769


# How many reference genes are not in spatial?
length(setdiff(ref_genes, spatial_genes))

# Optional: see examples of mismatches
head(setdiff(ref_genes, spatial_genes))

a <- setdiff(ref_genes, spatial_genes)

#-------------------- Further processing the reference ------------------------

# Load the Karthaus reference
karthaus_ref <- readRDS("C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Karthaus_reference.rds")

# Rename the identity: assign FullTypeMerged to a new metadata column 'celltype'
karthaus_ref$celltype <- karthaus_ref$FullTypeMerged
Idents(karthaus_ref) <- karthaus_ref$celltype

# Define the doublet labels to remove
doublet_labels <- c("PredDoublet_Str_Epi", "PredDoublet_Epi_Imm", "PredDoublet_Str_Imm")

# Remove all cells with any of these doublet labels
karthaus_ref <- subset(karthaus_ref, idents = doublet_labels, invert = TRUE)

# Drop any unused factor levels in the 'celltype' metadata
# karthaus_ref$celltype <- droplevels(karthaus_ref$celltype)
table(karthaus_ref$celltype)
# Log-normalize the data (this stores the normalized data in the "data" slot)
# karthaus_ref <- NormalizeData(karthaus_ref, normalization.method = "LogNormalize", scale.factor = 1e4)

karthaus_ref@misc$note <- "This object uses log2(TP10K + 1) values in the counts slot. No raw UMI counts present."

karthaus_ref@assays[["RNA_logTP10K"]] <- karthaus_ref@assays[["RNA"]]
karthaus_ref@assays[["RNA"]] <- NULL
DefaultAssay(karthaus_ref) <- "RNA_logTP10K"




# Save the cleaned and normalized object
save_file_updated <- "C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Karthaus_reference_cleaned.rds"
saveRDS(karthaus_ref, save_file_updated)
