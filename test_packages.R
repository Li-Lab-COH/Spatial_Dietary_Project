# Function to check if packages are available in the HPC cluster
check_packages_in_hpc <- function(package_list) {
  unavailable_packages <- c()  # Store missing packages
  
  for (package in package_list) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message(paste("Package", package, "is NOT available in the HPC cluster."))
      unavailable_packages <- c(unavailable_packages, package)
    } else {
      message(paste("Package", package, "is available in the HPC cluster."))
    }
  }
  
  # Return missing packages
  if (length(unavailable_packages) > 0) {
    message("The following packages are NOT available in the HPC cluster:")
    print(unavailable_packages)
  } else {
    message("All packages are available in the HPC cluster.")
  }
}

# Example package list (Replace with your actual list)
package_list <- c("Seurat", "ggplot2", "dplyr", "BiocManager", "remotes", 
                  "SeuratData", "SeuratWrappers", "DESeq2", "edgeR", "limma")

# Run the function to check package availability
check_packages_in_hpc(package_list)

