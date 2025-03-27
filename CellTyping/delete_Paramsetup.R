# How to grab the job array value
args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])  # This is your SLURM_ARRAY_TASK_ID
array_id <- 0

# Read CSV, skip the first column (param names), keep column headers
param_df <- read.csv("/home/janzules/Spatial/dietary_project/data/cell_typing_reference/rctd_parameter_grid.csv",
               row.names = 1, check.names = FALSE)

# param_df <- read.csv("~/Roselab/Spatial/dietary_project/data/cell_typing_reference/rctd_parameter_grid.csv",
#                      row.names = 1, check.names = FALSE)


# Assign parameter values
test_name <- colnames(param_df)[array_id + 1]
gene_cutoff <- as.numeric(param_df["gene_cutoff", array_id +1])
fc_cutoff <- as.numeric(param_df["fc_cutoff", array_id +1])
confidence_threshold <- as.numeric(param_df["confidence_threshold", array_id +1])
doublet_threshold <- as.numeric(param_df["doublet_threshold", array_id +1])

# Example: Print or use the parameters
cat("Running test", test_name, "\n")
cat("gene_cutoff:", gene_cutoff, "\n")
cat("fc_cutoff:", fc_cutoff, "\n")
cat("confidence_threshold:", confidence_threshold, "\n")
cat("doublet_threshold:", doublet_threshold, "\n")



# Create a named character vector of parameters
param_lines <- c(
  paste0("test_name: ", colnames(param_df)[array_id + 1]),
  paste0("gene_cutoff: ", gene_cutoff),
  paste0("fc_cutoff: ", fc_cutoff),
  paste0("confidence_threshold: ", confidence_threshold),
  paste0("doublet_threshold: ", doublet_threshold)
)

# Save to a file
writeLines(param_lines, paste0("~/Roselab/Spatial/dietary_project/data/" ,"rctd_params", ".txt"))





# Continue with RCTD setup using these values...

files <- readLines("~/Roselab/Spatial/dietary_project/code/addresses/Myc_Cap_cells.txt")

for (file in files){
  # print(file)
  sample_name <- gsub("\\.rds$", "", basename(file))
  print(sample_name)
}
