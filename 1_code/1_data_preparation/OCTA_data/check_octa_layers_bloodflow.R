library(tidyverse)
library(r4projects)
setwd(get_project_wd())

# Function to check all possible measurement types and layers in bloodflow data
check_bloodflow_patterns <- function() {
  cat("Starting bloodflow pattern check...\n")
  
  # Get all subjects
  all_subjects <- dir("2_data/OCTA/OCTA_bloodflow_2/")
  cat("Found", length(all_subjects), "subjects\n")
  
  # Initialize vectors to store patterns
  all_measures <- c()
  all_layers <- c()
  all_patterns <- c()
  
  # Process each subject
  for(subject in all_subjects) {
    cat("Checking subject:", subject, "\n")
    
    # Get all files in subject folder and subfolders
    files <- dir(
      file.path("2_data/OCTA/OCTA_bloodflow_2/", subject), 
      recursive = TRUE, 
      full.names = TRUE
    )
    
    # Filter for ETDRS files
    blood_files <- files[grep("ETDRS", files)]
    
    if(length(blood_files) > 0) {
      cat("  Found", length(blood_files), "bloodflow files\n")
      
      for(file in blood_files) {
        # Get base filename
        base_file <- basename(file)
        
        # Extract measure type (including parentheses part)
        measure <- str_extract(base_file, "[^_]+\\([^)]+\\)(?=_[^_]+_ETDRS)")
        if(!is.na(measure)) {
          all_measures <- c(all_measures, measure)
        }
        
        # Extract layer (part before ETDRS)
        layer <- str_extract(base_file, "(?<=_Angio)[^_]+(?=_ETDRS)")
        if(!is.na(layer)) {
          all_layers <- c(all_layers, layer)
        }
        
        # Store complete pattern for verification
        pattern <- paste0(measure, "_Angio", layer)
        if(!is.na(pattern)) {
          all_patterns <- c(all_patterns, pattern)
        }
      }
    }
  }
  
  # Get unique values
  unique_measures <- unique(all_measures)
  unique_layers <- unique(all_layers)
  unique_patterns <- unique(all_patterns)
  
  cat("\n=== RESULTS ===\n")
  cat("\nFound", length(unique_measures), "unique measurement types:\n")
  for(measure in sort(unique_measures)) {
    cat("- ", measure, "\n")
  }
  
  cat("\nFound", length(unique_layers), "unique layers:\n")
  for(layer in sort(unique_layers)) {
    cat("- ", layer, "\n")
  }
  
  # Count frequencies
  measure_freq <- table(all_measures)
  layer_freq <- table(all_layers)
  
  cat("\nMeasurement type frequencies:\n")
  print(measure_freq)
  
  cat("\nLayer frequencies:\n")
  print(layer_freq)
  
  return(list(
    measures = unique_measures,
    layers = unique_layers,
    patterns = unique_patterns,
    measure_frequencies = measure_freq,
    layer_frequencies = layer_freq
  ))
}

# Run the check
results <- check_bloodflow_patterns()
