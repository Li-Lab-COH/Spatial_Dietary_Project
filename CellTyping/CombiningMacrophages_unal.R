MyC_CaP_ref <- readRDS("~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_unal_reference.rds")

table(MyC_CaP_ref$celltype)

# Removing 0 annotation

MyC_CaP_ref$celltype <- droplevels(factor(MyC_CaP_ref$celltype))

MyC_CaP_ref$celltype <- factor(MyC_CaP_ref$celltype)
MyC_CaP_ref$celltype <- droplevels(MyC_CaP_ref$celltype)

Idents(MyC_CaP_ref) <- "celltype"


# Merging

# Define the cell types to merge
macrophage_types <- c(
  "ISG-high Macrophages",
  "Macrophages",
  "Proliferating Macrophages",
  "Resident Macrophages",
  "M2 Macrophages"
)

# Replace them in the metadata
MyC_CaP_ref$celltype <- as.character(MyC_CaP_ref$celltype)  # ensure character
MyC_CaP_ref$celltype[MyC_CaP_ref$celltype %in% macrophage_types] <- "Macrophages"

# Convert back to factor and clean up levels
MyC_CaP_ref$celltype <- factor(MyC_CaP_ref$celltype)
MyC_CaP_ref$celltype <- droplevels(MyC_CaP_ref$celltype)

# Set as Idents again
Idents(MyC_CaP_ref) <- "celltype"

DimPlot(MyC_CaP_ref, group.by = "celltype", reduction = "umap", label = TRUE, pt.size = 0.5)



#Combining nk, t, and tregs

# Define the cell types to merge
nk_tcell_types <- c("NK Cells", "T Cells", "Tregs")

# Convert to character (if needed)
MyC_CaP_ref$celltype <- as.character(MyC_CaP_ref$celltype)

# Relabel those types
MyC_CaP_ref$celltype[MyC_CaP_ref$celltype %in% nk_tcell_types] <- "NK_n_Tcells"

# Convert back to factor and drop unused levels
MyC_CaP_ref$celltype <- factor(MyC_CaP_ref$celltype)
MyC_CaP_ref$celltype <- droplevels(MyC_CaP_ref$celltype)

# Set as Idents again
Idents(MyC_CaP_ref) <- "celltype"
table(Idents(MyC_CaP_ref))

# Visualizing LumP

# Create a new column for plotting
MyC_CaP_ref$LumP_highlight <- ifelse(MyC_CaP_ref$celltype == "LumP", "LumP", "Other")

# Plot with custom colors
DimPlot(MyC_CaP_ref, group.by = "LumP_highlight", reduction = "umap", 
        cols = c("red", "gray80"), pt.size = 0.5)

# Dropping LumP

MyC_CaP_ref <- subset(MyC_CaP_ref, subset = celltype != "LumP")

MyC_CaP_ref$celltype <- droplevels(MyC_CaP_ref$celltype)
Idents(MyC_CaP_ref) <- "celltype"

table(Idents(MyC_CaP_ref))


#Visualizing

# Create a data frame from the identity table
cell_counts <- as.data.frame(table(Idents(MyC_CaP_ref)))
colnames(cell_counts) <- c("CellType", "Count")

# Filter out "MyC-CaP"
cell_counts <- subset(cell_counts, CellType != "MyC-CaP")

# Plot
library(ggplot2)
ggplot(cell_counts, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Cell Count per Cell Type (Excluding MyC-CaP)",
       x = "Cell Type",
       y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



table(Idents(MyC_CaP_ref))

# Standardizing the names
# Replace spaces with underscores in the celltype metadata
MyC_CaP_ref$celltype <- gsub(" ", "_", MyC_CaP_ref$celltype)

# Set the Idents to the updated celltype column
Idents(MyC_CaP_ref) <- MyC_CaP_ref$celltype

# Confirm the change
unique(Idents(MyC_CaP_ref))


saveRDS(MyC_CaP_ref, "~/Roselab/Spatial/dietary_project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_condensed.rds")
