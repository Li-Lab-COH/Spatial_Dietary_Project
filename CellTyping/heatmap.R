# Load necessary libraries
library(pheatmap) # for simple heatmap
# or library(ComplexHeatmap) if you want a more flexible one

# Load the data
cell_counts <- read.csv("~/Roselab/Spatial/dietary_project/data/CellTyping_objects/Cell_counts_before_condensing.csv", row.names = 1)

# Convert to matrix
count_mat <- as.matrix(cell_counts)

# Row-wise z-score normalization
zscore_mat <- t(scale(t(count_mat)))  # scale() works column-wise, so we transpose twice

# Plot heatmap
pheatmap(zscore_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
