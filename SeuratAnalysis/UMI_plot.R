library(Seurat)
library(dplyr)
library(ggplot2)

files <- list.files("~/Roselab/Spatial/dietary_project/data/Raw_objects",
                    pattern="\\.rds$", full.names=TRUE)

df_list <- lapply(files, function(f) {
  obj <- readRDS(f)
  if ("nCount_Spatial" %in% colnames(obj@meta.data)) {
    umi <- obj@meta.data$nCount_Spatial
  } else if ("nCount_RNA" %in% colnames(obj@meta.data)) {
    umi <- obj@meta.data$nCount_RNA
  } else {
    umi <- colSums(GetAssayData(obj, slot="counts"))
  }
  data.frame(sample = tools::file_path_sans_ext(basename(f)), umi = umi)
})

umi_df <- bind_rows(df_list)


umi_df$sample <- sub("^F\\d{5}_", "", umi_df$sample)# example filter

umi_df <- umi_df %>% filter(umi >= 50)

desired_order <- c(
  "2_KDRT", "3_RT",   "5_RT",   "7_CRRT", "8_CRRT", "9_KDRT",
  "10_CRRT","11_RT",  "13_KDRT","27_LFRT","28_LFRT","29_LFRT",
  "31_RT",  "32_RT",  "33_LGIRT","34_LGIRT","35_LGIRT",
  "36_CRRT","37_CRRT","38_CRRT"
)

# apply to your data frame
umi_df$sample <- factor(umi_df$sample, levels = desired_order)

# then plot
ggplot(umi_df, aes(x = sample, y = umi)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    draw_quantiles = c(0.25, 0.5, 0.75),
    fill = "lightblue",
    color = "black"
  ) +
  theme_bw() +
  ylim(0, 1000) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Total UMI counts")

# ggplot(umi_df, aes(x = sample, y = umi)) +
#   geom_violin(trim = FALSE, scale = "width") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   ylim(0, 1000)
