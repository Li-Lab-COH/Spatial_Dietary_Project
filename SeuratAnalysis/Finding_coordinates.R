# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(Seurat)  # for GetTissueCoordinates()

# Set paths
rds_dir <- "C:/Users/jonan/Documents/1Work/RoseLab/Spatial/dietary_droject/data/Raw_objects"
qc_file <- file.path(rds_dir, "qc_metrics_all_edited.csv")
qc_metrics <- read_csv(qc_file)

# Function to get bounding box from Seurat Visium object
get_bbox <- function(obj) {
  coords <- SeuratObject::GetTissueCoordinates(obj, scale = NULL)
  tibble(
    xmin = min(coords$x, na.rm = TRUE),
    xmax = max(coords$x, na.rm = TRUE),
    ymin = min(coords$y, na.rm = TRUE),
    ymax = max(coords$y, na.rm = TRUE)
  )
}

# List all .rds files
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# Loop through RDS files and extract information
coords_list <- lapply(rds_files, function(path) {
  fn      <- basename(path) %>% str_remove("\\.rds$")
  fnumber <- str_extract(fn, "^F\\d{5}")
  sample  <- str_extract(fn, "(?<=F\\d{5}_).+$")
  
  obj  <- readRDS(path)
  bbox <- get_bbox(obj)
  
  tibble(
    label       = fn,
    Fnumber     = fnumber,
    sample_name = sample,
    xmin        = bbox$xmin,
    xmax        = bbox$xmax,
    ymin        = bbox$ymin,
    ymax        = bbox$ymax
  )
})

# Combine and join with QC table
coords_df <- bind_rows(coords_list)
final_df <- qc_metrics %>%
  left_join(coords_df, by = "label")

# Save result
write_csv(final_df, file.path(rds_dir, "sample_coords.csv"))



obj_test <- readRDS(rds_files[3])
a <- SeuratObject::GetTissueCoordinates(obj_test, scale = NULL)
head(a)
get 

fn      <- basename(rds_files[3]) %>% str_remove("\\.rds$")
fnumber <- str_extract(fn, "^F\\d{5}")
sample  <- str_extract(fn, "(?<=F\\d{5}_).+$")
