library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot) 
library(ggplot2)
library(stringr)

#---------------Preparing data-----------

marker_genes <- read.csv("../../data/ClusterGenes/ms_31RT.csv")
all_clusters <- sort(unique(as.integer(markers$cluster)))[1]


dim(marker_genes)
dim(marker_genes[marker_genes$p_val_adj < 0.05, ])

summary(marker_genes$avg_log2FC)


#----------------Testing loop-------------
test_group <- subset(marker_genes, cluster == all_clusters )
dim(test_group)
dim(test_group[test_group$p_val_adj < 0.05, ])

summary(test_group$p_val_adj)
gene_symbols <- test_group$gene
cleaned_gene_symbols <- gsub("\\.m0$", "", gene_symbols)

ego <- enrichGO(
  gene = cleaned_gene_symbols,  # No need to manually convert
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",           # Tell enrichGO that input genes are symbols
  ont = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)


head(ego)
top_ego1 <- ego
top_ego1@result <- top_ego1@result %>%
  arrange(p.adjust, desc(FoldEnrichment), desc(Count)) %>%
  head(15)

# Now dotplot() will work because top_ego1 is still an enrichResult object
dotplot(top_ego1, showCategory = 15)
dotplot(ego, showCategory = 15)

barplot(ego, showCategory = 15)
#--------------------Actual Loop---------------------

for (clust in all_clusters) {
  cluster_subset <- subset(markers, cluster == clust & p_val_adj < 0.05 & avg_log2FC > 0.25)
  gene_symbols   <- cluster_subset$gene
  
  entrez_ids <- mapIds(
    org.Mm.eg.db,
    keys = gene_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  entrez_ids <- na.omit(entrez_ids)
  
  # Overrepresentation analysis
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05
  )
  
  # Print or store the result
  print(paste("Cluster", clust))
  print(head(ego, 10))  # top 10 terms
}


#---------------------Actual Loop 2-------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

# Define the directory to save figures
save_dir <- "../../figures/ms_31RT/"

# # Ensure the directory exists
# if (!dir.exists(save_dir)) {
#   dir.create(save_dir, recursive = TRUE)  # Create nested directories if needed
# }

# Extract all unique clusters
all_clusters <- unique(marker_genes$cluster)

# Initialize list to store dotplots
dotPlot_GO_cluster <- list()

# Loop through each cluster
for (clust in all_clusters) {
  
  # Subset marker genes for the current cluster
  cluster_subset <- subset(marker_genes, cluster == clust)
  
  # Filter for significant genes
  sig_genes <- cluster_subset[cluster_subset$p_val_adj < 0.05, ]
  
  # Extract gene symbols and clean them
  gene_symbols <- sig_genes$gene
  cleaned_gene_symbols <- gsub("\\.m0$", "", gene_symbols)
  
  # Check if there are genes to analyze
  if (length(cleaned_gene_symbols) > 0) {
    
    # Perform GO enrichment analysis
    ego <- enrichGO(
      gene = cleaned_gene_symbols,  
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",           
      ont = "BP", 
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05
    )
    
    # Ensure GO results are not empty before plotting
    if (nrow(ego@result) > 0) {
      
      # Keep only the top 15 results
      top_ego <- ego
      top_ego@result <- top_ego@result %>%
        arrange(p.adjust, desc(FoldEnrichment), desc(Count)) %>%
        head(15)
      
      # Generate dotplot
      
      
      dot_plot <- dotplot(top_ego, showCategory = 15) +
        ggtitle(paste("GO Enrichment - Cluster", clust)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +  # Wrap text at 30 characters
        theme(
          axis.text.y = element_text(size = 5)  # Reduce text size
        )
      
      
      # Store plot in list
      dotPlot_GO_cluster[[paste0("dotPlot_GO_cluster", clust)]] <- dot_plot
      
      # Define full save path
      file_path <- paste0(save_dir, "GO_Cluster_", clust, ".png")
      
      # Save each plot to the specified directory
      ggsave(filename = file_path, plot = dot_plot, width = 6, height = 4)
      
    } else {
      message(paste("No enriched GO terms found for Cluster", clust, "- Skipping..."))
    }
    
  } else {
    message(paste("No significant genes for Cluster", clust, "- Skipping..."))
  }
}

# View stored plots in R
dotPlot_GO_cluster


