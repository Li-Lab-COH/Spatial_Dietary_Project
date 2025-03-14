library(dplyr)
library(Seurat)
library(ggplot2)

ST_data <- F07836_13_KDRT
rm(F07836_13_KDRT)




metadata <- ST_data@meta.data

nCount_summary <- metadata %>%
  group_by(cluster = BANKSY_snn_res.0.5) %>%
  summarise(mean = mean(nCount_Spatial.008um), 
            min = min(nCount_Spatial.008um),
            max = max(nCount_Spatial.008um),
            Q1 = quantile(nCount_Spatial.008um, 0.25, na.rm = TRUE),
            Q3 = quantile(nCount_Spatial.008um, 0.75, na.rm = TRUE)
            ) %>%
  arrange(min, mean)



VlnPlot(ST_data, features = "nCount_Spatial.008um", 
        group.by = "BANKSY_snn_res.0.5", pt.size = 0, raster=FALSE)+
  geom_hline(yintercept = quantile(metadata$nCount_Spatial.008um, 0.05),
             color = "#DC143C", linewidth = 1.5)

quantile(metadata$nCount_Spatial.008um, 0.05)
summary(ST_data$nCount_Spatial.008um)
View(nCount_summary)

VlnPlot(ST_data, features = "nCount_Spatial.008um", 
        group.by = "BANKSY_snn_res.0.5", pt.size = 0, raster=FALSE)+
  geom_hline(yintercept = 70, 
             color = "#DC143C", linewidth = 1.5)

ls()
