library(Seurat)
library(dplyr)
library(tibble)  # for rownames_to_column()


input_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/data/Rogelio/Analysis_031625/BANKSY_Normalized_QC_Filtered"

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)




tumor_clusters_list <- list(
  "F07833_5_RT" = c(6, 9),
  "F07834_28_LFRT" = c(3, 16),
  "F07833_8_CRRT" = c(1, 13),
  "F07833_9_KDRT" = c(2, 10),
  "F07833_28_LFRT" = c(3),
  "F07834_5_RT" = c(3, 6, 9),
  "F07834_8_CRRT" = c(3, 4),
  "F07834_9_KDRT" = c(3, 4, 0),
  "F07834_28_LFRT" = c(3, 16, 7),
  "F07835_7_CRRT" = c(2, 5, 7, 9),
  "F07836_11_RT" = c(8, 9),
  "F07837_31_RT" = c(0, 2, 6, 8, 9),
  "F07837_37_CR_RT" = c(8, 9, 11),
  "F07838_32_RT" = c(2, 4, 10, 11),
  "F07838_36_CR_RT" = c(4, 5),
  "F07838_38_CR_RT" = c(2, 7)
)

### 1. Set up object and tumor cell identification ###
# Set your Seurat object and assay
obj <- F07834_28_LFRT  # update with your object name if needed
DefaultAssay(obj) <- "BANKSY"

# Identify cells from clusters 3, 16, and 7 (for tumor classification)
target_clusters <- c("3", "16", "7")
cells_to_check <- rownames(obj@meta.data[obj@meta.data$BANKSY_snn_res.0.5 %in% target_clusters, ])

# Define tumor gene list and fetch normalized data
tumor_gene_list <- c("Myc", "Ar", "Ccna2", "Ccnd1", "Cdk1", "Odc1", "Slc7a5", "Ldha", "Cad")
expr_data <- FetchData(obj, vars = tumor_gene_list, slot = "data")
# In case some genes come with a prefix (e.g., "spatial008um_"), remove it:
colnames(expr_data) <- gsub("^spatial008um_", "", colnames(expr_data))

# Subset to cells from the target clusters
expr_data_subset <- expr_data[cells_to_check, ]

# Apply your tumor criteria (here simplified to any expression > 0 as an example)
tumor_cells <- expr_data_subset %>%
  filter(
    (Myc > 0 | Ar > 0) &    # at least one of these markers
      rowSums(.[, c("Ccna2", "Ccnd1", "Cdk1", "Odc1", "Slc7a5", "Ldha", "Cad")] > 0) >= 1
  ) %>%
  rownames()

# Create a new metadata column "cell_types" and assign "Tumor" to these cells;
# all others get "Other" by default.
obj@meta.data$cell_types <- "Other"
obj@meta.data[tumor_cells, "cell_types"] <- "Tumor"

### 2. Expand gene list to include immune markers ###
# Extend the gene list to include markers for immune cell subtypes.
immune_gene_list <- c(
  "Trac", "Trbc1", "Trbc2", "Cd3e", "Cd4", "Cd8a", "Foxp3",      # T cell markers (fix Trbc1)
  "Cd19", "Cd22", "Ighm",                                        # B cell markers
  "Itgax", "Batf3", "H2-Ab1",                                    # Dendritic cell markers (fix Batf3, H2-Ab1)
  "Ncr1",                                                        # NK cell marker
  "Cxcr3", "Il12b", "Nos2", "Tnf", "Ccl5",                       # N1 neutrophil markers
  "Tgfb1", "Vegfa", "Arg1", "Cd274", "Mmp9",                     # pro-tumor markers for neutrophils (fix Tgfb1)
  "Retnlg", "Cxcr2",                                             # N2 markers
  "Cd163", "Mrc1", "Il10", "Ccl22"                               # M2 markers
)


# Fetch the extended expression data.
# (If these genes are present with a prefix, the same gsub() cleaning applies.)
expr_data_ext <- FetchData(obj, vars = immune_gene_list, slot = "data")
colnames(expr_data_ext) <- gsub("^spatial008um_", "", colnames(expr_data_ext))

### 3. Classify immune cells among non-tumor cells ###
# Work only with cells not already labeled as "Tumor"
non_tumor_ids <- rownames(obj@meta.data)[obj@meta.data$cell_types == "Other"]
immune_expr <- expr_data_ext[non_tumor_ids, ] %>% rownames_to_column("cell")

# Use rowwise() and case_when() to assign immune subtypes.
# The order in case_when determines priority.
immune_expr <- immune_expr %>%
  rowwise() %>%
  mutate(
    immune_type = case_when(
      # T cell criteria: any one of T cell markers >0, then subtype based on further markers.
      (Foxp3 > 0) ~ "Treg",
      (Cd4 > 0) ~ "Cd4",
      (Cd8a > 0) ~ "Cd8",
      
      # B cells: require at least one marker above 0
      (Cd19 > 0) ~ "B cell",
      
      # Dendritic cells: either Itgax or Batf3 > 0
      (Itgax > 0) ~ "Dendritic cell",
      
      # NK cells: Ncr1 > 0.5 (lowered from 1.5)
      (Ncr1 > 0.5) ~ "NK cell",
      
      # Neutrophils:
      # N1 (anti-tumor, pro-inflammatory): all five markers above new thresholds
      (Cxcr3 > 1 & Il12b > 1 & Nos2 > 1.5 & Tnf > 1 & Ccl5 > 1) ~ "N1",
      
      # Pro-tumorigenic neutrophils: at least 2 of the set [Tgfb1, Vegfa, Arg1, Cd274, Mmp9] are hit
      (( (Tgfb1 > 2) + (Vegfa > 1.75) + (Arg1 > 1.75) + (Cd274 > 1.75) + (Mmp9 > 1.5) ) >= 2 & 
         (Retnlg > 1 | Cxcr2 > 1.75)) ~ "N2",
      
      # M2 Neutrophils: at least 2 of [Cd163, Mrc1, Il10, Ccl22] required
      (( (Tgfb1 > 2) + (Vegfa > 1.75) + (Arg1 > 1.75) + (Cd274 > 1.75) + (Mmp9 > 1.5) ) >= 2 & 
         ((Cd163 > 1.75) + (Mrc1 > 2) + (Il10 > 1.75) + (Ccl22 > 1)) >= 2) ~ "M2",
      
      # If at least 2 pro-tumor markers are hit but no further specification, call it "Suppressor"
      (( (Tgfb1 > 2) + (Vegfa > 1.75) + (Arg1 > 1.75) + (Cd274 > 1.75) + (Mmp9 > 1.5) ) >= 3) ~ "Suppressor",
      
      # Default if no immune criteria match
      TRUE ~ "Other"
    )
  ) %>%
  ungroup()


# Now, update the metadata for these non-tumor cells.
# We use the cell IDs (in column "cell") to assign the immune type.
immune_assignments <- immune_expr %>% select(cell, immune_type)
obj@meta.data$cell_types[match(immune_assignments$cell, rownames(obj@meta.data))] <- immune_assignments$immune_type

### 4. Update active identities and create a summary table ###
# Set the Seurat object's active identity to the new cell_types
Idents(obj) <- obj@meta.data$cell_types

# Create a data frame that counts cells in each group
cell_counts <- as.data.frame(table(obj@meta.data$cell_types))
colnames(cell_counts) <- c("Cell_Type", "Count")
print(cell_counts)

# Define the expected list of all possible cell types
expected_cell_types <- c(
  "B cell", "Cd4", "Cd8", "Dendritic cell", "NK cell", 
  "N1", "N2", "M2", "Suppressor", "Treg", "Tumor", "Other"
)

# Get actual counts from the metadata
cell_counts <- as.data.frame(table(obj@meta.data$cell_types))
colnames(cell_counts) <- c("Cell_Type", "Count")

# Create a complete dataframe ensuring all expected cell types are present
all_cell_counts <- data.frame(Cell_Type = expected_cell_types) 

# Merge with actual counts, filling missing ones with zero
all_cell_counts <- merge(all_cell_counts, cell_counts, by = "Cell_Type", all.x = TRUE)

# Replace NAs (missing counts) with zeros
all_cell_counts$Count[is.na(all_cell_counts$Count)] <- 0

# Print the updated cell count table
print(all_cell_counts)








