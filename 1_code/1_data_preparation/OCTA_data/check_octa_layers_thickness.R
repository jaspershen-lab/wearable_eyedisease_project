library(tidyverse)
library(r4projects)
setwd(get_project_wd())

# Function to check all possible layer names
check_layers <- function() {
  cat("Starting layer check...\n")
  
  # Get all subjects
  all_subjects <- dir("2_data/OCTA/OCTA_thickness/")
  cat("Found", length(all_subjects), "subjects\n")
  
  all_layers <- c()
  
  # Process each subject
  for(subject in all_subjects) {
    cat("Checking subject:", subject, "\n")
    
    # Get all files in subject folder and subfolders
    files <- dir(
      file.path("2_data/OCTA/OCTA_thickness/", subject), 
      recursive = TRUE, 
      full.names = TRUE
    )
    
    # Filter for thickness files
    thickness_files <- files[grep("Thickness.*ETDRS", files)]
    
    if(length(thickness_files) > 0) {
      cat("  Found", length(thickness_files), "thickness files\n")
      
      for(file in thickness_files) {
        # Extract layer name
        layer <- str_extract(basename(file), "(?<=Thickness_)[^_]+(?=_ETDRS)")
        if(!is.na(layer)) {
          all_layers <- c(all_layers, layer)
        }
      }
    }
  }
  
  # Get unique layer names
  unique_layers <- unique(all_layers)
  
  cat("\n=== RESULTS ===\n")
  cat("Found", length(unique_layers), "unique layer names:\n")
  for(layer in sort(unique_layers)) {
    cat("- ", layer, "\n")
  }
  
  # Count frequency of each layer
  layer_freq <- table(all_layers)
  cat("\nLayer frequencies:\n")
  print(layer_freq)
  
  return(list(
    unique_layers = unique_layers,
    frequencies = layer_freq
  ))
}

# Run the check
results <- check_layers()