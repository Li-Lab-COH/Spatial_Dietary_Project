# runSpatialAnnotation.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)         # for ImmGenData()
library(SummarizedExperiment)
library(data.table)
library(Matrix)
library(openxlsx)

# Define the function
runSpatialAnnotation <- function(seurat_obj, output_folder, 
                                 k_params = c(10, 20, 30), 
                                 resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2)) {
  # Create output folder if it doesn't exist
  ensure_dir <- function(dir_path) {
    dir_path <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
  ensure_dir(output_folder)
  
  # Loop over clustering parameters
  for (k in k_params) {
    for (res in resolutions) {
      
      # Create folder for this parameter combination
      comb_folder <- file.path(output_folder, paste0("neighbors_", k, "_resolution_", res))
      ensure_dir(comb_folder)
      
      # Work on a fresh copy of the object
      st_obj <- seurat_obj
      
      # Run neighbor search and clustering using your PCA reduction (assumed precomputed)
      st_obj <- FindNeighbors(
        st_obj,
        reduction = "pca.prostate.full",   # Ensure your object has this reduction
        dims = 1:15,
        k.param = k,
        prune.SNN = 1/15,
        verbose = FALSE
      )
      st_obj <- FindClusters(
        st_obj,
        resolution = res,
        verbose = FALSE
      )
      
      # Run UMAP using the PCA reduction
      st_obj <- RunUMAP(
        st_obj,
        dims = 1:15,
        reduction = "pca.prostate.full",
        verbose = FALSE
      )
      
      #----------------------- UMAP Plot -----------------------
      p_umap <- DimPlot(
        st_obj,
        reduction = "umap",
        label = TRUE
      ) + ggtitle(paste("UMAP: k =", k, "resolution =", res))
      
      ggsave(
        filename = file.path(comb_folder, paste0("UMAP_clusters_k", k, "_res", res, ".png")),
        plot = p_umap,
        width = 8,
        height = 6
      )
      
      #----------------------- Violin Plot -----------------------
      Idents(st_obj) <- "seurat_clusters"
      p_vln <- VlnPlot(
        st_obj,
        features = "nCount_Spatial.008um",
        group.by = "seurat_clusters",
        pt.size = 0
      ) + ggtitle(paste("Violin plot (nCount) k =", k, "resolution =", res))
      
      ggsave(
        filename = file.path(comb_folder, paste0("VlnPlot_nCount_Spatial_k", k, "_res", res, ".png")),
        plot = p_vln,
        width = 8,
        height = 6
      )
      
      #----------------------- Marker Identification -----------------------
      markers <- FindAllMarkers(
        st_obj,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25
      )
      top20_markers_all <- markers %>%
        group_by(cluster) %>%
        top_n(n = 20, wt = avg_log2FC)
      
      write.csv(
        top20_markers_all,
        file = file.path(comb_folder, paste0("top20_markers_k", k, "_res", res, ".csv")),
        row.names = FALSE
      )
      
      #----------------------- SingleR Annotation -----------------------
      # Extract normalized expression matrix and cluster identities
      expr_mat <- GetAssayData(st_obj, slot = "data")
      clusters <- Idents(st_obj)
      
      ## Helper function to process SingleR scores into top-n predictions
      process_top_predictions <- function(preds, n = 7) {
        scores <- preds$scores
        rownames(scores) <- rownames(preds)
        top_matches <- lapply(rownames(scores), function(cluster_id) {
          x <- scores[cluster_id, ]
          sorted <- sort(x, decreasing = TRUE)
          data.frame(
            Top1 = names(sorted)[1], Score1 = sorted[1],
            Top2 = names(sorted)[2], Score2 = sorted[2],
            Top3 = names(sorted)[3], Score3 = sorted[3],
            Top4 = names(sorted)[4], Score4 = sorted[4],
            Top5 = names(sorted)[5], Score5 = sorted[5],
            Top6 = names(sorted)[6], Score6 = sorted[6],
            Top7 = names(sorted)[7], Score7 = sorted[7],
            row.names = cluster_id,
            stringsAsFactors = FALSE
          )
        })
        top_matches_df <- do.call(rbind, top_matches)
        top_matches_df$Cluster <- rownames(top_matches_df)
        # Reorder columns to have Cluster first
        top_matches_df <- top_matches_df[, c("Cluster",
                                             "Top1", "Score1", "Top2", "Score2", "Top3", "Score3",
                                             "Top4", "Score4", "Top5", "Score5", "Top6", "Score6",
                                             "Top7", "Score7")]
        return(top_matches_df)
      }
      
      # --- (A) ImmGen Reference ---
      ref_immgen <- ImmGenData()
      preds_immgen <- SingleR(
        test = expr_mat,
        ref = ref_immgen,
        labels = ref_immgen$label.fine,
        clusters = clusters
      )
      top_matches_df_immgen <- process_top_predictions(preds_immgen, n = 7)
      write.csv(
        top_matches_df_immgen,
        file = file.path(comb_folder, paste0("singleR_top7_ImmGen_k", k, "_res", res, ".csv")),
        row.names = FALSE
      )
      
      # --- (B) Karhtaus Reference ---
      ref_kar <- readRDS("C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Karthaus_2020/Karthaus_reference_cleaned.rds")
      # For Karhtaus, the log-normalized data is stored in the counts slot
      ref_expr <- GetAssayData(ref_kar, slot = "counts")
      ref_labels <- ref_kar$celltype
      shared_genes <- intersect(rownames(expr_mat), rownames(ref_expr))
      expr_mat_sub <- expr_mat[shared_genes, , drop = FALSE]
      ref_expr_sub <- ref_expr[shared_genes, , drop = FALSE]
      
      preds_kar <- SingleR(
        test = expr_mat_sub,
        ref = ref_expr_sub,
        labels = ref_labels,
        clusters = clusters
      )
      top_matches_df_kar <- process_top_predictions(preds_kar, n = 7)
      write.csv(
        top_matches_df_kar,
        file = file.path(comb_folder, paste0("singleR_top7_Karthaus_k", k, "_res", res, ".csv")),
        row.names = FALSE
      )
      
      # --- (C) Unal Reference ---
      ref_unal <- readRDS("C://Users/jonan/Documents/1Work/RoseLab/Spatial/Dietary_Project/data/cell_typing_reference/Unal_2024_Myc_CaP/Fully_annotated_condensed.rds")
      # For Unal, the normalized data is stored in the 'data' slot
      ref_expr <- GetAssayData(ref_unal, slot = "data")
      ref_labels <- ref_unal$celltype
      shared_genes <- intersect(rownames(expr_mat), rownames(ref_expr))
      expr_mat_sub <- expr_mat[shared_genes, , drop = FALSE]
      ref_expr_sub <- ref_expr[shared_genes, , drop = FALSE]
      
      preds_unal <- SingleR(
        test = expr_mat_sub,
        ref = ref_expr_sub,
        labels = ref_labels,
        clusters = clusters
      )
      top_matches_df_unal <- process_top_predictions(preds_unal, n = 7)
      write.csv(
        top_matches_df_unal,
        file = file.path(comb_folder, paste0("singleR_top7_Unal_k", k, "_res", res, ".csv")),
        row.names = FALSE
      )
      
      #----------------------- Combine CSV files into an Excel Workbook -----------------------
      csv_files <- list.files(comb_folder, pattern = "\\.csv$", full.names = TRUE)
      wb <- createWorkbook()
      for (file_path in csv_files) {
        sheet_name <- tools::file_path_sans_ext(basename(file_path))
        sheet_name <- substr(sheet_name, 1, 31)  # Limit sheet name to 31 characters
        data <- read.csv(file_path)
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, data)
      }
      excel_output_path <- file.path(comb_folder, "combined_output.xlsx")
      saveWorkbook(wb, file = excel_output_path, overwrite = TRUE)
      
      message(paste("Completed parameter combination: k =", k, "resolution =", res))
    }
  }
}
