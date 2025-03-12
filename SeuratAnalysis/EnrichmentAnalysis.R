library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot) 
library(ggplot2)
library(stringr)
library(dplyr)



set.seed(1337)
################################################################
# Loops the clusters of a given file and performs an enrichment analysis
# for each cluster. This analysis looks for GO and KEGG pathways that are
# over represented in each cluster. Previous scripts should have generated
# 


#-----------------------------Analysis Parameters------------------------------

# Top how many pathways to display
top_number = 15

# Input location and file name
# Output location
# TODO: Need to make sure that plots have the mouse label on the title
# Define the directory to save figures

# save_dir_GO <- "../../figures/ms_31RT/GO/"
# save_dir_KEGG <- "../../figures/ms_31RT/KEGG/"

save_dir_GO <- "../../figures/ms_5_RT/GO/"
save_dir_KEGG <- "../../figures/ms_5_RT/KEGG/"

# Create directories if they do not exist
dir.create(save_dir_GO, recursive = TRUE, showWarnings = FALSE)
dir.create(save_dir_KEGG, recursive = TRUE, showWarnings = FALSE)



#---------------Preparing data-----------

# marker_genes <- read.csv("../../data/ClusterGenes/ms_31RT.csv")
marker_genes <- read.csv("../../data/ClusterGenes/ms_5_RT.csv")
all_clusters <- sort(unique(as.integer(markers$cluster)))[1]


dim(marker_genes)
dim(marker_genes[marker_genes$p_val_adj < 0.05, ])

summary(marker_genes$avg_log2FC)

#-----------------------------Loop 2 with KEGG--------------------------------


# TODO: Make a check for both folders, include a failure of the entire scripts
# if it doesn't exist
#
# Ensure the directory exists
# if (!dir.exists(save_dir)) {
#   dir.create(save_dir, recursive = TRUE)  
# }

# Extract all unique clusters
all_clusters <- unique(marker_genes$cluster)


# Loop through each cluster
for (clust in all_clusters) {
  
  print("")
  print("")
  print(paste("Processing Cluster:", clust))
  
  # Subset marker genes for the current cluster
  cluster_subset <- subset(marker_genes, cluster == clust)
  
  # Filter for significant genes
  sig_genes <- cluster_subset[cluster_subset$p_val_adj < 0.05, ]
  
  # Extract gene symbols and clean them
  gene_symbols <- sig_genes$gene
  cleaned_gene_symbols <- gsub("\\.m0$", "", gene_symbols)
  
  # Convert Gene Symbols to Entrez IDs for KEGG analysis
  entrez_ids <- mapIds(org.Mm.eg.db, 
                       keys = cleaned_gene_symbols, 
                       column = "ENTREZID", 
                       keytype = "SYMBOL", 
                       multiVals = "first")
  
  # Remove NAs
  entrez_ids <- na.omit(entrez_ids)
  
  # Ensure there are genes to analyze
  if (length(cleaned_gene_symbols) > 0) {
    
    print(paste("Performing GO analysis for cluster: ", clust))
    # Perform GO enrichment
    ego <- enrichGO(
      gene = cleaned_gene_symbols,  
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",           
      ont = "BP", 
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05
    )
    
    # Perform KEGG enrichment (Only if there are valid Entrez IDs)
    if (length(entrez_ids) > 0) {
      print(paste("Performing KEGG analysis for cluster: ", clust))
      kegg <- enrichKEGG(
        gene = as.character(entrez_ids),  # KEGG requires Entrez IDs
        organism = "mmu",  # mmu = Mus musculus (mouse)
        pvalueCutoff = 0.05
      )
    } else{ message(paste("No Entrez ID - Cluster: ", clust))   }
    
    # Ensure GO results are not empty before plotting
    if (nrow(ego@result) > 0) {
      
      # Keep only the top top_number results
      top_ego <- ego
      top_ego@result <- top_ego@result %>%
        arrange(p.adjust, desc(FoldEnrichment), desc(Count)) %>%
        head(top_number)
      
      # Wrap the text in the Description column BEFORE plotting. Removes annoying warning
      top_ego@result$Description <- str_wrap(top_ego@result$Description, width = 30)
      
      # Generate GO dotplot
      dot_plot_go <- dotplot(top_ego, showCategory = top_number) +
          ggtitle(paste("GO Enrichment - Cluster", clust)) +
          theme(
            axis.text.y = element_text(size = 5)  # Reduce text size
          )
      
      # # Generate GO dotplot
      # dot_plot_go <- dotplot(top_ego, showCategory = top_number) +
      #   ggtitle(paste("GO Enrichment - Cluster", clust)) +
      #   theme(
      #     axis.text.y = element_text(size = 5)  # Reduce text size
      #   ) +
      #   scale_y_discrete(labels = function(x) str_wrap(x, width = 30))  # Wrap text at 30 characters
      
      
      
      # Save GO plot
      ggsave(filename = paste0(save_dir_GO, "GO_Cluster_", clust, ".png"), 
             plot = dot_plot_go, width = 6, height = 4)
    } else{message(paste("No ego Results - Cluster: ", clust))}
    
    # Ensure KEGG results are not empty before plotting
    if (exists("kegg") && dim(kegg)[1] > 0) {
      
      # Remove " - Mus musculus (house mouse)" from pathway names
      kegg@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg@result$Description)
      
      # Keep only the top top_number results
      top_kegg <- kegg
      top_kegg@result <- top_kegg@result %>%
        arrange(p.adjust) %>%
        head(top_number)
      
      # # Generate KEGG dotplot
      # dot_plot_kegg <- dotplot(top_kegg, showCategory = top_number) +
      #   ggtitle(paste("KEGG Pathway Enrichment - Cluster", clust)) +
      #   theme(
      #     axis.text.y = element_text(size = 5)  # Reduce text size
      #   ) +
      #   scale_y_discrete(labels = function(x) str_wrap(x, width = 30))  # Wrap text at 30 characters
      
      # Modify the labels directly in `top_kegg` (if possible)
      top_kegg@result$Description <- str_wrap(top_kegg@result$Description, width = 30)
      
      dot_plot_kegg <- dotplot(top_kegg, showCategory = top_number) +
          ggtitle(paste("KEGG Pathway Enrichment - Cluster", clust)) +
          theme(
            axis.text.y = element_text(size = 5)  # Reduce text size
          )
      
      # Save KEGG plot
      ggsave(filename = paste0(save_dir_KEGG, "KEGG_Cluster_", clust, ".png"), 
             plot = dot_plot_kegg, width = 6, height = 4)
    } else{message(paste("No KEGG Results - Cluster: ", clust))}
    
  } else {
    message(paste("No significant genes for Cluster", clust, "- Skipping..."))
  }
}


