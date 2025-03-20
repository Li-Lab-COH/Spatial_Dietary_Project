library(Seurat)
library(dplyr)
library(tibble)
library(pryr)  # for mem_used()
library(future)
library(FSA)  # For Dunn’s post-hoc test
library(tidyverse)
library(openxlsx) 

# Prevent unwanted parallelization
plan("sequential")

# Define the expected immune cell types for the output table
expected_cell_types <- c("B cell", "Cd4", "Cd8", "Dendritic cell", "Macrophage", "M1", "M2", "NK cell", 
                         "N1", "N2", "Suppressor", "Treg", "Luminal", "Prostate_Cell", "Active_Cell", "Other")

# Specify the input directory containing your .rds files
input_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/data/Rogelio/Analysis_031625/BANKSY_Normalized_QC_Filtered"

# Initialize a master counts data frame
master_counts <- data.frame(Cell_Type = expected_cell_types)

# List of sample names
sample_names <- c(
  "F07833_5_RT", "F07834_28_LFRT", "F07833_8_CRRT", "F07833_9_KDRT", 
  "F07834_5_RT", "F07834_8_CRRT", "F07834_9_KDRT", "F07835_7_CRRT", 
  "F07836_11_RT", "F07837_31_RT", "F07837_37_CR_RT", "F07838_32_RT", 
  "F07838_36_CR_RT", "F07838_38_CR_RT",
  "F07833_28_LFRT", "F07835_27_LFRT", "F07835_2_KDRT", "F07835_3_RT",
  "F07836_10_CRRT", "F07836_13_KDRT",
  "F07836_29_LFRT", "F07837_33_LGI_RT", "F07837_35_LGI_RT",
  "F07837_37_CR_RT", "F07838_34_LGI_RT", "F07838_36_CR_RT", "F07838_38_CR_RT"
)
sample_names
# sample_names <- sample_names[1:2]

# Loop over each sample
for (sample_name in sample_names) {
  
  file_path <- file.path(input_dir, paste0(sample_name, ".rds"))
  
  if (!file.exists(file_path)) {
    cat("File not found for sample", sample_name, "\n")
    next
  }
  
  cat("Processing sample:", sample_name, "\n")
  seurat_obj <- readRDS(file_path)
  
  ### 1. Immune & Tumor Cell Classification ###
  DefaultAssay(seurat_obj) <- "BANKSY"
  
  immune_gene_list <- c(
    # T cell markers
    "Trac", "Trbc1", "Trbc2", "Cd3e", "Cd4", "Cd8a", "Foxp3",      
    
    # B cell markers
    "Cd19", "Cd22", "Ighm",                                         
    
    # Dendritic cell markers
    "Itgax", "Batf3", "H2-Ab1",                                     
    
    # NK cell marker
    "Ncr1",                                                         
    
    # N1 neutrophil markers (strengthened specificity)
    "Cxcr3", "Il12b", "Nos2", "Tnf", "Ccl5", "S100a8", "S100a9", "Cd177",
    
    # N2 neutrophil markers (improved exclusivity)
    "Tgfb1", "Vegfa", "Arg1", "Cd274", "Mmp9", "Retnlg", "Cxcr2", "Camp", "Il1rn",
    
    # M2 macrophage markers (new specificity markers added)
    "Cd163", "Mrc1", "Il10", "Ccl22", "C1qa", "C1qb", "C1qc", "Mgl2",
    
    # Macrophage markers
    "Adgre1", "Cd68", "Itgam", "Fcgr1", "Lyz2", "Csf1r", "Marco",   
    
    # M1 macrophage marker (for better differentiation)
    "Il1b",
    
    # Suppressor/MDSC markers
    "Arg1", "Nos2",
    
    # Luminal markers
    "Cldn10", "Ppp1r1b",
    
    # Prostate cell markers
    "Ar", "Myc",
    
    # Active cell markers (proliferative cells)
    "Ccna2", "Rcc1", "Cdk1", "Ldha", "Cad", "Cdkn2a", "Mki67"
  )
  
  expr_data <- FetchData(seurat_obj, vars = immune_gene_list, slot = "data")
  colnames(expr_data) <- gsub("^spatial008um_", "", colnames(expr_data))
  
  # Convert row names to a column for easier manipulation
  immune_expr <- expr_data %>%
    rownames_to_column("cell")
  
  immune_expr <- expr_data %>%
    rownames_to_column("cell") %>%
    rowwise() %>%
    mutate(
      immune_type = case_when(
        # **M1 Macrophages**
        (Adgre1 > 0 | Cd68 > 0 | Itgam > 0) & (Csf1r > 0 | Marco > 0) & (Tnf > 0 | Il1b > 0) ~ "M1",
        
        # **General Macrophages (Including M2)**
        (Adgre1 > 0 | Cd68 > 0 | Itgam > 0) & (Csf1r > 0 | Marco > 0) ~ "Macrophage",
        
        # Regulatory T Cells (Tregs)
        (Foxp3 > 0) ~ "Treg",
        
        # CD4 and CD8 T Cells
        (Cd4 > 0) ~ "Cd4",
        (Cd8a > 0) ~ "Cd8",
        
        # B Cells
        (Cd19 > 0) ~ "B cell",
        
        # Dendritic Cells
        (Itgax > 0) ~ "Dendritic cell",
        
        # Natural Killer Cells
        (Ncr1 > 0.5) ~ "NK cell",
        
        # **N1 Neutrophils** - More inclusive thresholds
        (S100a8 > 0.1 | S100a9 > 0.1) & (Cd177 > 0.1 | Nos2 > 0.25) & 
          (Cxcr3 > 0.1 | Il12b > 0.1 | Tnf > 0.1 | Ccl5 > 0.1) ~ "N1",
        
        # **N2 Neutrophils** - Looser thresholds, still requiring at least one granulocyte marker
        ((Tgfb1 > 0.5) + (Vegfa > 0.5) + (Arg1 > 0.5) + (Cd274 > 0.5) + (Mmp9 > 0.5)) >= 1 &
          (Retnlg > 0.1 | Cxcr2 > 0.25) & (Camp > 0.1 | Il1rn > 0.1) ~ "N2",
        
        # **M2 Macrophages** - More flexible, allowing for weaker expression
        ((Tgfb1 > 1) + (Vegfa > 1) + (Arg1 > 1) + (Cd274 > 1) + (Mmp9 > 1)) >= 2 &
          ((Cd163 > 1) + (Mrc1 > 1.5) + (Il10 > 1) + (Ccl22 > 1) + (C1qa > 0.5 | C1qb > 0.5 | C1qc > 0.5) + (Mgl2 > 0.5)) >= 2 ~ "M2",
        
        # **Suppressor Cells (MDSCs, Tregs, TAMs)** - More inclusive
        (Arg1 > 1 & Nos2 > 0.5) | ((Tgfb1 > 1) + (Vegfa > 1) + (Cd274 > 1)) >= 2 ~ "Suppressor",
        
        # **Luminal Cells**
        (Cldn10 > 0 & Ppp1r1b > 0) ~ "Luminal",
        
        # **Prostate Cells**
        (Ar > 0 & Myc > 0) ~ "Prostate_Cell",
        
        # **Active Cells (Highly proliferative cells)**
        (Ccna2 > 0 | Rcc1 > 0 | Cdk1 > 0 | Ldha > 0 | Cad > 0 | Cdkn2a > 0 | Mki67 > 0) ~ "Active_Cell",
        
        TRUE ~ "Other"
      )
    ) %>%
    ungroup()
  
  immune_assignments <- immune_expr %>% select(cell, immune_type)
  seurat_obj@meta.data$cell_types <- "Other"
  seurat_obj@meta.data$cell_types[match(immune_assignments$cell, rownames(seurat_obj@meta.data))] <- immune_assignments$immune_type
  
  Idents(seurat_obj) <- seurat_obj@meta.data$cell_types
  
  cell_counts <- as.data.frame(table(seurat_obj@meta.data$cell_types))
  colnames(cell_counts) <- c("Cell_Type", "Count")
  
  all_cell_counts <- data.frame(Cell_Type = expected_cell_types)
  all_cell_counts <- merge(all_cell_counts, cell_counts, by = "Cell_Type", all.x = TRUE)
  all_cell_counts$Count[is.na(all_cell_counts$Count)] <- 0
  
  sample_counts <- all_cell_counts %>% select(Cell_Type, Count)
  colnames(sample_counts)[2] <- sample_name
  master_counts <- merge(master_counts, sample_counts, by = "Cell_Type", all = TRUE)
  
  ### 3. Clean Up Memory ###
  cat("Memory before clearing object for", sample_name, ":", pryr::mem_used(), "\n")
  
  # Explicitly remove all large objects
  rm(seurat_obj, expr_data, immune_expr, immune_assignments, cell_counts, all_cell_counts)
  
  # Aggressive garbage collection
  gc(reset = TRUE)
  invisible(gc())
  
  cat("Memory after clearing object for", sample_name, ":", pryr::mem_used(), "\n\n")
}

print(master_counts)
# Identify base column names without .x or .y
base_names <- gsub("\\.x$|\\.y$", "", colnames(master_counts))  

# Find duplicate base names
duplicate_cols <- names(which(table(base_names) > 1))

# Keep only the first occurrence of each duplicate column
master_counts_cleaned <- master_counts[, !duplicated(base_names)]
# If needed, save the cleaned dataframe
# write.csv(master_counts_cleaned, "cleaned_master_counts.csv", row.names = FALSE)

# Define the file path
output_path_master_counts <- "/Users/janzules/Roselab/Spatial/dietary_project/data/master_counts.csv"

# Write the dataframe to a CSV file
write.csv(master_counts_cleaned, file = output_path_master_counts, row.names = FALSE)

master_counts <- master_counts_cleaned

#---------------------------- Creating immune fraction df ---------------------

# List of immune cell types to keep
immune_cell_types <- c("B cell", "Cd4", "Cd8", "Dendritic cell", 
                       "M1", "M2", "Macrophage", "N1", "N2", 
                       "NK cell", "Suppressor", "Treg")

# Filter master_counts to keep only immune cell types
immune_counts <- master_counts[master_counts$Cell_Type %in% immune_cell_types, ]

# Calculate total immune cells for each column
total_immune_cells <- colSums(immune_counts[, -1], na.rm = TRUE)

# Convert counts into fractions of total immune cells
immune_fractions <- immune_counts

for (col in colnames(immune_fractions)[-1]) {
  immune_fractions[[col]] <- immune_counts[[col]] / total_immune_cells[col]
}

# Print the transformed dataframe
print(immune_fractions)



# Define the file path
output_path_immune_fractions <- "/Users/janzules/Roselab/Spatial/Dietary_project/Results/immune_cell_fractions/immune_fractions.csv"

# Write the dataframe to a CSV file
write.csv(immune_fractions, file = output_path_immune_fractions, row.names = FALSE)
#-------------------------- Updated stats ----------------------------------
# Define sample groups
control_samples <- c("F07833_5_RT", "F07834_5_RT", "F07836_11_RT", "F07837_31_RT", "F07838_32_RT", "F07835_3_RT", "F07837_37_CR_RT")
LFRT_samples <- c("F07834_28_LFRT", "F07833_28_LFRT", "F07835_27_LFRT", "F07836_29_LFRT")
LGI_RT_samples <- c("F07837_33_LGI_RT", "F07837_35_LGI_RT", "F07838_34_LGI_RT")
CRRT_samples <- c("F07833_8_CRRT", "F07834_8_CRRT", "F07835_7_CRRT", "F07836_10_CRRT", "F07837_37_CR_RT", "F07838_36_CR_RT", "F07838_38_CR_RT")
KDRT_samples <- c("F07833_9_KDRT", "F07834_9_KDRT", "F07835_2_KDRT", "F07836_13_KDRT")



# List of immune cell types
immune_cell_types <- c("B cell", "Cd4", "Cd8", "Dendritic cell", 
                       "M1", "M2", "Macrophage", "N1", "N2", 
                       "NK cell", "Suppressor", "Treg")

# Create a condition mapping for each sample
sample_conditions <- data.frame(
  Sample = c(control_samples, LFRT_samples, LGI_RT_samples, CRRT_samples, KDRT_samples),
  Condition = c(rep("Control", length(control_samples)),
                rep("LFRT", length(LFRT_samples)),
                rep("LGI_RT", length(LGI_RT_samples)),
                rep("CRRT", length(CRRT_samples)),
                rep("KDRT", length(KDRT_samples)))
)

# Convert master_counts to long format
immune_fractions_long <- immune_fractions %>%
  pivot_longer(cols = -Cell_Type, names_to = "Sample", values_to = "Fraction") %>%
  inner_join(sample_conditions, by = "Sample")  # Add condition labels

# Initialize results storage
kruskal_results <- list()
dunn_results <- list()

# Perform Kruskal-Wallis test for each immune cell type
for (cell in immune_cell_types) {
  subset_data <- immune_fractions_long %>% filter(Cell_Type == cell)
  
  # Kruskal-Wallis test
  kw_test <- kruskal.test(Fraction ~ Condition, data = subset_data)
  
  # Store results
  kruskal_results[[cell]] <- kw_test
  
  # If significant, perform Dunn's test for pairwise comparisons
  if (kw_test$p.value < 0.05) {
    dunn_test <- dunnTest(Fraction ~ Condition, data = subset_data, method = "bonferroni")
    dunn_results[[cell]] <- dunn_test
  }
}

# Print Kruskal-Wallis results
print(kruskal_results)

# Print Dunn’s post-hoc results (only if significant)
print(dunn_results)

# Visualize immune fractions across conditions
ggplot(immune_fractions_long, aes(x = Condition, y = Fraction, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Immune Cell Fractions Across Conditions",
       x = "Condition",
       y = "Fraction of Total Immune Cells")

subset_nk <- immune_fractions_long %>% filter(Cell_Type == "NK cell")
dunn_nk <- dunnTest(Fraction ~ Condition, data = subset_nk, method = "bonferroni")

# Print results
print(dunn_nk)

# Define file paths
output_path_dunn_results <- "/Users/janzules/Roselab/Spatial/Dietary_project/Results/immune_cell_fractions/dunn_results.csv"
output_path_dunn_nk <- "/Users/janzules/Roselab/Spatial/Dietary_project/Results/immune_cell_fractions/dunn_nk.csv"

# Convert Dunn's test results to a data frame
dunn_results_df <- do.call(rbind, lapply(dunn_results, function(x) x$res))
dunn_nk_df <- as.data.frame(dunn_nk$res)

# Write to CSV
write.csv(dunn_results_df, file = output_path_dunn_results, row.names = FALSE)
write.csv(dunn_nk_df, file = output_path_dunn_nk, row.names = FALSE)


# MORE
subset_treg <- immune_fractions_long %>% filter(Cell_Type == "Treg")
dunnTest(Fraction ~ Condition, data = subset_treg, method = "bonferroni")

subset_m2 <- immune_fractions_long %>% filter(Cell_Type == "M2")
dunnTest(Fraction ~ Condition, data = subset_m2, method = "bonferroni")



subset_cd4 <- immune_fractions_long %>% filter(Cell_Type == "Cd4")
dunnTest(Fraction ~ Condition, data = subset_cd4, method = "bonferroni")



immune_of_interest <- c("Bcell", "Cd4", "Cd8", "Dendritic cell",
                        "M1", "M2", "NK cell", "Suppressor", "Treg",
                        "Macrophage")




# Define the list of immune cells of interest
immune_of_interest <- c("B cell", "Cd4", "Cd8", "Dendritic cell",
                        "M1", "M2", "NK cell", "Suppressor", "Treg",
                        "Macrophage")

# Define output directory
output_dir <- "/Users/janzules/Roselab/Spatial/Dietary_project/Results/immune_cell_fractions/"

# Ensure directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each immune cell type and perform Dunn's test
for (immune_cell in immune_of_interest) {
  
  # Subset the data for the specific immune cell type
  subset_data <- immune_fractions_long %>% filter(Cell_Type == immune_cell)
  
  # Perform Dunn's test
  dunn_result <- dunnTest(Fraction ~ Condition, data = subset_data, method = "bonferroni")
  
  # Convert results to a data frame
  dunn_df <- as.data.frame(dunn_result$res)
  
  # Define output file path
  output_file <- file.path(output_dir, paste0("dunn_", immune_cell, ".xlsx"))
  
  # Save the result as an Excel file
  write.xlsx(dunn_df, file = output_file, rowNames = FALSE)
  
  # Print confirmation
  cat("Saved:", output_file, "\n")
}





















# ----------------- Stats now----------------#

# Vectors for your groups
control_samples <- c("F07833_5_RT", "F07834_5_RT", "F07836_11_RT", "F07837_31_RT", "F07838_32_RT")
control_samples2 <- c("F07833_5_RT", "F07834_5_RT", "F07836_11_RT", "F07837_31_RT", "F07838_32_RT", "F07834_8_CRRT", "F07835_7_CRRT")
cr_samples <- c("F07837_37_CR_RT", "F07838_36_CR_RT", "F07838_38_CR_RT")
cr_samples <- c("F07833_9_KDRT", "F07833_9_KDRT")
cr_samples <- c("F07834_8_CRRT", "F07835_7_CRRT", "F07837_37_CR_RT", "F07838_36_CR_RT", "F07838_38_CR_RT")
# Get all column names excluding "Cell_Type"
all_samples <- colnames(master_counts)[colnames(master_counts) != "Cell_Type"]

# Identify non-control (CR) samples
cr_samples <- setdiff(all_samples, control_samples2)

library(dplyr)
library(tidyr)

# Assuming your dataframe is called master_counts
# 1) Subset columns: Cell_Type plus the 10 relevant samples
cols_of_interest <- c("Cell_Type", control_samples, cr_samples)

df_subset <- master_counts[, cols_of_interest]

# 2) Pivot to long format
df_long <- df_subset %>%
  pivot_longer(
    cols = -Cell_Type,
    names_to = "Sample",
    values_to = "Count"
  )

# 3) Add a "Group" column (Control vs CR)
df_long <- df_long %>%
  mutate(
    Group = case_when(
      Sample %in% control_samples ~ "Control",
      Sample %in% cr_samples     ~ "CR",
      TRUE                       ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))  # keep only the 10 relevant samples



# Extract tumor counts for each sample
tumor_counts <- df_long %>%
  filter(Cell_Type == "Tumor") %>%
  select(Sample, Tumor_Count = Count, Group)

# Merge tumor counts back with df_long
df_long <- df_long %>%
  left_join(tumor_counts, by = c("Sample", "Group")) %>%
  filter(Cell_Type != "Tumor")  # We won't do ratio for Tumor itself

# Compute ratio & log-transform
df_long <- df_long %>%
  mutate(
    Ratio = (Count + 0.0) / (Tumor_Count + 1),  # +1 to avoid divide by zero
    Log_Ratio = log2(Ratio + 1)
  )

# run test
wilcox_results <- df_long %>%
  group_by(Cell_Type) %>%
  summarise(
    p_value = wilcox.test(Log_Ratio ~ Group)$p.value
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")  # multiple-testing correction
  ) %>%
  arrange(p_value)

wilcox_results


# Define the immune cell types you want to sum
immune_types <- c("B cell", "Cd4", "Cd8", "Dendritic cell", "NK cell")

# Summarize the sum of these immune cells per sample
df_total_immune <- df_long %>%
  filter(Cell_Type %in% immune_types) %>%
  group_by(Sample, Group, Tumor_Count) %>%
  summarise(Sum_Immune = sum(Count)) %>%
  ungroup() %>%
  # Create ratio & log-ratio
  mutate(
    Immune_Ratio = Sum_Immune / (Tumor_Count + 1),
    Log_Immune_Ratio = log2(Immune_Ratio + 1),
    Cell_Type = "Total_Immune"  # label it
  )

# We can do a Wilcoxon test just on this total_immune
wilcox_total_immune <- wilcox.test(Immune_Ratio ~ Group, data = df_total_immune)
wilcox_total_immune



library(ggplot2)

# Boxplot for each cell type
ggplot(df_long, aes(x = Group, y = Ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  theme_minimal() +
  ggtitle("Immune Cell Enrichment (Log2 Ratio to Tumor)")

# Boxplot just for the total_immune measure
ggplot(df_total_immune, aes(x = Group, y = Log_Immune_Ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_minimal() +
  ggtitle("Total Immune Cell Enrichment (Log2 Ratio to Tumor)")



#----------------------kruskal wallis----------------------------
# Define sample groups
control_samples <- c("F07833_5_RT", "F07834_5_RT", "F07836_11_RT", "F07837_31_RT", "F07838_32_RT")

# All non-control samples get a group label based on sample names
df_long <- df_long %>%
  mutate(
    Detailed_Group = case_when(
      Sample %in% control_samples ~ "Control",
      grepl("CRRT", Sample) ~ "CRRT",
      grepl("KDRT", Sample) ~ "KDRT",
      grepl("LFRT", Sample) ~ "LFRT",
      grepl("CR_RT", Sample) ~ "CR_RT",
      TRUE ~ "Other"
    )
  )


# Run Kruskal-Wallis test for each immune cell type
kruskal_results <- df_long %>%
  group_by(Cell_Type) %>%
  summarise(
    p_value = kruskal.test(Ratio ~ Detailed_Group)$p.value
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

# Display results
print(kruskal_results)

kruskal_total_immune <- kruskal.test(Immune_Ratio ~ Detailed_Group, data = df_total_immune)
print(kruskal_total_immune)


#-----------------------------processing cell numbers--------------------------
# Define immune-related cell types
immune_cell_types <- c("B cell", "Cd4", "Cd8", "Dendritic cell", "M1", "M2", 
                       "Macrophage", "N1", "N2", "NK cell", "Suppressor", "Treg")

# Define "Other" cells (immune + other specific non-immune cells)
other_cell_types <- c(immune_cell_types, "Active_Cell", "Prostate_Cell", "Luminal")

# Compute the Total_Immune row
total_immune <- colSums(master_counts[master_counts$Cell_Type %in% immune_cell_types, -1], na.rm = TRUE)

# Compute the Total_Other_Cells row (immune cells + Active_Cell + Prostate_Cell)
total_other_cells <- colSums(master_counts[master_counts$Cell_Type %in% other_cell_types, -1], na.rm = TRUE)

# Create new rows for total immune and total other cells
total_immune_row <- c("Total_Immune", total_immune)
total_other_row <- c("Total_Other_Cells", total_other_cells)

# Add new rows to the master_counts data frame
master_counts <- rbind(master_counts, total_immune_row, total_other_row)

# Ensure numeric columns remain numeric
master_counts[, -1] <- lapply(master_counts[, -1], as.numeric)

# Print the updated table
print(master_counts)


# cat("Processing complete. Immune and tumor cell counts saved.\n")
# print(master_counts)

