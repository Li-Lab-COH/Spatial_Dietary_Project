library(Seurat)
library(dplyr)
library(readr)
library(stringr)

# Example metadata: spatial folders and their associated barcode CSVs
metadata <- data.frame(
  sample_id = c(
    "F07833",
    "F07834",
    "F07835",
    "F07836",
    "F07837",
    "F07838"
  ),
  spatial_dir = c(
    "/Users/janzules/Roselab/Spatial/dietary_project/data/SpacerangerOut/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07833_22WJCYLT3/outs",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/SpacerangerOut/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07834_22WJCYLT3/outs",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/SpacerangerOut/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07835_22WJCYLT3/outs",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/SpacerangerOut/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07836_22WJCYLT3/outs",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/SpacerangerOut/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07837_22WJCYLT3/outs",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/SpacerangerOut/BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_F07838_22WJCYLT3/outs"
  ),
  barcode_csv = c(
    "/Users/janzules/Roselab/Spatial/dietary_project/data/ExportedSectionLabels/F07833.csv",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/ExportedSectionLabels/F07834.csv",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/ExportedSectionLabels/F07835.csv",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/ExportedSectionLabels/F07836.csv",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/ExportedSectionLabels/F07837.csv",
    "/Users/janzules/Roselab/Spatial/dietary_project/data/ExportedSectionLabels/F07838.csv"
  ),
  stringsAsFactors = FALSE
)

# Output folder for the Seurat objects
output_dir <- "/Users/janzules/Roselab/Spatial/dietary_project/data/Raw_objects"
dir.create(output_dir, showWarnings = FALSE)

# Initialize an empty QC metrics data frame.
# Each row will be labelled with "sampleID_tissueLabel" and contain various metrics.
qc_metrics_all <- data.frame(
  label               = character(),
  total_bins          = integer(),
  fraction_under_100  = numeric(),
  fraction_zero       = numeric(),
  mean_umi            = numeric(),
  median_umi          = numeric(),
  sd_umi              = numeric(),
  mean_percent_mt     = numeric(),
  median_percent_mt   = numeric(),
  bins_gt_100         = integer(),
  stringsAsFactors    = FALSE
)

# Loop through each sample
for (i in seq_len(nrow(metadata))) {
  
  sample_id <- metadata$sample_id[i]
  spatial_path <- metadata$spatial_dir[i]
  barcode_path <- metadata$barcode_csv[i]
  
  print(sample_id)
  print(spatial_path)
  print(barcode_path)
  cat("Processing sample:", sample_id, "\n")
  
  # Load the full spatial object
  print("Load10X_Spatial()")
  seurat_obj <- Load10X_Spatial(
    data.dir = spatial_path,
    bin.size = 8,
    slice = sample_id,
  )
  # Add percent.mt column to metadata (assumes mouse gene names)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # Load barcode table (assumed to have at least two columns: 'barcode' and 'tissue')
  cat("reading: ", barcode_path)
  barcode_df <- read_csv(barcode_path)
  
  # Loop over unique tissues defined in the barcode file
  for (tissue_label in unique(barcode_df$Treatment)) {
    
    cat("  Subsetting tissue:", tissue_label, "\n")
    
    tissue_barcodes <- barcode_df %>%
      filter(Treatment == tissue_label) %>%
      pull(Barcode)
    
    if (length(tissue_barcodes) == 0) {
      warning(paste0("No barcodes matched for ", subset_label, ". Skipping."))
      next
    }
    
    # Subset the Seurat object based on the barcodes for the current tissue
    subset_obj <- subset(seurat_obj, cells = tissue_barcodes)
    
    # Save the subset object with a descriptive name
    out_filename <- file.path(
      output_dir,
      paste0(sample_id, "_", tissue_label, ".rds")
    )
    saveRDS(subset_obj, file = out_filename)
    
    # Extract metadata for QC metrics. We use the nCount_Spatial column as umi_values.
    meta_data <- subset_obj@meta.data
    umi_values <- meta_data$nCount_Spatial
    
    # Compute QC metrics
    total_bins <- nrow(meta_data)
    fraction_under_100 <- mean(umi_values < 100)
    fraction_zero <- mean(umi_values == 0)
    mean_umi <- mean(umi_values)
    median_umi <- median(umi_values)
    sd_umi <- sd(umi_values)
    bins_gt_100 <- sum(umi_values > 100)
    
    if ("percent.mt" %in% colnames(meta_data)) {
      mean_percent_mt <- mean(meta_data$percent.mt)
      median_percent_mt <- median(meta_data$percent.mt)
    } else {
      mean_percent_mt <- NA
      median_percent_mt <- NA
    }
    
    # Create a label for this subset
    subset_label <- paste0(sample_id, "_", tissue_label)
    
    # Append these metrics as a new row to the QC metrics dataframe
    qc_metrics_all <- rbind(
      qc_metrics_all,
      data.frame(
        label = subset_label,
        total_bins = total_bins,
        fraction_under_100 = fraction_under_100,
        fraction_zero = fraction_zero,
        mean_umi = mean_umi,
        median_umi = median_umi,
        sd_umi = sd_umi,
        mean_percent_mt = mean_percent_mt,
        median_percent_mt = median_percent_mt,
        bins_gt_100 = bins_gt_100,
        stringsAsFactors = FALSE
      )
    )
    print(qc_metrics_all)
  }
}

# Print out the combined QC metrics
print(qc_metrics_all)

# Optionally, save the QC metrics to a CSV file
write.csv(qc_metrics_all, file = file.path(output_dir, "qc_metrics_all.csv"), row.names = FALSE)
