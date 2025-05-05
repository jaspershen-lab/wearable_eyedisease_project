library(tidyverse)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Function to extract measurement info from filename
extract_measure_info <- function(filename) {
  # Extract layer name using all possible layers we found
  layer <- str_extract(
    filename,
    paste0(
      "(",
      paste(
        "Choroid",
        "GCL\\+IPL",
        "ILM\\+BM",
        "INL",
        "OuterRetina",
        "PED",
        "Retina",
        "RNFL",
        "RNFL\\+GCL\\+IPL",
        sep = "|"
      ),
      ")"
    )
  )
  
  # Extract measurement type (always Thickness in this case)
  measure_type <- "Thickness"
  
  list(measure_type = measure_type, layer = layer)
}

# Initialize empty list for all results
result_list <- list()

# Main processing
all_subjects <- dir("2_data/OCTA/OCTA_thickness_2/")

for(i in 1:length(all_subjects)) {
  cat("\nProcessing subject", i, "(", all_subjects[i], ")\n")
  current_id <- paste0("SH", stringr::str_extract(all_subjects[i], "[0-9]{3}"))
  
  # Get all 6x6 files for this subject
  all_files <- dir(file.path("2_data/OCTA/OCTA_thickness_2/", all_subjects[i]))
  all_files <- all_files[stringr::str_detect(all_files, "26x21")]
  
  if(length(all_files) == 0) next
  
  # Create time mapping for this subject
  dates <- str_extract(all_files, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  unique_dates <- sort(unique(dates))
  time_points <- paste0("T", seq_along(unique_dates) - 1)
  names(time_points) <- unique_dates
  
  subject_results <- list()
  
  for(temp_folder in all_files) {
    current_date <- str_extract(temp_folder, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    which_eye <- str_extract(temp_folder, "O[DS]{1}")
    current_time <- time_points[current_date]
    
    # Get all thickness files
    quant_path <- file.path("2_data/OCTA/OCTA_thickness_2/", 
                            all_subjects[i], temp_folder, "Quantization")
    all_data <- dir(quant_path)
    all_data <- all_data[str_detect(all_data, "Thickness") & str_detect(all_data, "ETDRS")]
    
    if(length(all_data) == 0) next
    
    # Process each data file
    for(temp_data in all_data) {
      tryCatch({
        # Extract measurement info
        measure_info <- extract_measure_info(temp_data)
        if(is.na(measure_info$layer)) next
        
        # Read and process data
        data <- read_csv(
          file.path(quant_path, temp_data),
          locale = locale(encoding = "UTF-8"),
          col_names = FALSE,
          show_col_types = FALSE
        )
        
        # Process data
        result <- data[-1,] %>%
          dplyr::select(1:2) %>%
          setNames(c("region", "value")) %>%
          mutate(
            region = case_when(
              grepl("0-21", region) ~ "0_21",
              grepl("0-3", region) ~ "0_3",
              grepl("0-6", region) ~ "0_6",
              grepl("0-9", region) ~ "0-9",
              grepl("0-12", region) ~ "0_12",
              grepl("0-15", region) ~ "0-15",
              TRUE ~ NA_character_
            ),
            value = as.numeric(as.character(value))
          ) %>%
          filter(!is.na(region), !is.na(value))
        
        if (nrow(result) == 0) next
        
        # Create wide format column names with time point
        for(j in 1:nrow(result)) {
          col_name <- paste(
            measure_info$measure_type,
            measure_info$layer,
            which_eye,
            result$region[j],
            current_time,
            sep = "_"
          )
          
          # Store result with subject ID
          if(length(subject_results) == 0) {
            subject_results[[1]] <- tibble(id = current_id)
          }
          subject_results[[1]][[col_name]] <- result$value[j]
        }
        
      }, error = function(e) {
        cat("Error processing file:", e$message, "\n")
      })
    }
  }
  
  # Add subject results to main list
  if(length(subject_results) > 0) {
    result_list[[i]] <- subject_results[[1]]
  }
}

# Combine all results
final_data <- bind_rows(result_list) %>%
  arrange(.data$id)

###create working directory
if (!dir.exists(
  "2_data/analysis_data"
)) {
  dir.create(
    "2_data/analysis_data",
    showWarnings = FALSE,
    recursive = TRUE
  )
}

setwd(
  "2_data/analysis_data/"
)

# Save results
write_csv(final_data, "octa_data_thickness.csv")

# Display results
cat("\nFinal data dimensions:", nrow(final_data), "x", ncol(final_data), "\n")
cat("\nFirst few rows:\n")
print(head(final_data))

colnames(final_data)
