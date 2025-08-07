# ----------------------------------------------------
# OCTA + Vision Combined Clustering Analysis for PPV Group
# Using Mfuzz Fuzzy Clustering with OCTA (blood flow + thickness) and vision parameters
# MODIFIED: Only using WF region (0_21) OCTA parameters
# MODIFIED: Added LogMAR conversion for vision parameters before visualization
# ----------------------------------------------------

# Load required libraries
library(tidyverse)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(factoextra)
library(corrplot)
library(r4projects)
library(pheatmap)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# -------------------- 0. Vision conversion functions --------------------
# LogMAR conversion functions
convert_to_logmar <- function(decimal_vision) {
  # Convert decimal vision (1.0, 0.8, etc.) to LogMAR
  # LogMAR = -log10(decimal_vision)
  # Handle edge cases: values <= 0 or > 2.0
  
  logmar_vision <- ifelse(decimal_vision > 0 & decimal_vision <= 2.0,
                          -log10(decimal_vision),
                          NA)
  
  # For very poor vision (< 0.05 decimal), cap at LogMAR 1.3
  logmar_vision <- ifelse(!is.na(logmar_vision) & logmar_vision > 1.3,
                          1.3,
                          logmar_vision)
  
  return(logmar_vision)
}

convert_vision_parameters <- function(data) {
  # Convert vision parameters to LogMAR before analysis
  cat("Converting vision parameters from decimal to LogMAR...\n")
  
  vision_cols <- c("pre_vision", "post_vision_1m", "vision_improvement_1m")
  
  # Store original values for reference
  for(col in vision_cols) {
    if(col %in% names(data)) {
      data[[paste0(col, "_decimal")]] <- data[[col]]
    }
  }
  
  # Convert baseline and post-operative vision to LogMAR
  if("pre_vision" %in% names(data)) {
    data$pre_vision_logmar <- convert_to_logmar(data$pre_vision)
    cat("Pre-vision: converted", sum(!is.na(data$pre_vision)), "values to LogMAR\n")
  }
  
  if("post_vision_1m" %in% names(data)) {
    data$post_vision_1m_logmar <- convert_to_logmar(data$post_vision_1m)
    cat("Post-vision: converted", sum(!is.na(data$post_vision_1m)), "values to LogMAR\n")
  }
  
  # Recalculate improvement in LogMAR (negative change = improvement)
  if("pre_vision_logmar" %in% names(data) & "post_vision_1m_logmar" %in% names(data)) {
    data$vision_improvement_1m_logmar <- data$pre_vision_logmar - data$post_vision_1m_logmar
    cat("Vision improvement: recalculated", sum(!is.na(data$vision_improvement_1m_logmar)), "LogMAR improvement values\n")
  }
  
  # Update vision parameters list to use LogMAR versions
  # Replace original parameters with LogMAR versions for analysis
  data$pre_vision <- data$pre_vision_logmar
  data$post_vision_1m <- data$post_vision_1m_logmar  
  data$vision_improvement_1m <- data$vision_improvement_1m_logmar
  
  cat("✓ Vision parameters converted to LogMAR for analysis\n\n")
  
  return(data)
}

# Function to add proper labels for vision plots
get_vision_plot_labels <- function(parameter_name) {
  # Return proper labels for vision parameters
  labels <- list(
    "pre_vision" = list(
      title = "Pre-operative Vision (LogMAR)",
      subtitle = "Lower LogMAR = better vision",
      y_label = "LogMAR (Lower = Better)"
    ),
    "post_vision_1m" = list(
      title = "Post-operative Vision at 1 Month (LogMAR)",
      subtitle = "Lower LogMAR = better vision",
      y_label = "LogMAR (Lower = Better)"
    ),
    "vision_improvement_1m" = list(
      title = "Vision Improvement at 1 Month (LogMAR)",
      subtitle = "Positive LogMAR change = improvement",
      y_label = "LogMAR Change (Positive = Improvement)"
    ),
    "age" = list(
      title = "Age",
      subtitle = "Patient age distribution",
      y_label = "Age (years)"
    )
  )
  
  if(parameter_name %in% names(labels)) {
    return(labels[[parameter_name]])
  } else {
    return(list(
      title = parameter_name,
      subtitle = "",
      y_label = parameter_name
    ))
  }
}

# -------------------- 1. Load data --------------------
# Load baseline information and OCTA data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster")

# -------------------- 2. Process OCTA data --------------------
# Function to process OCTA data
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

# Process blood flow and thickness data
octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

# Function to process patient data
process_patient_data <- function(patient_data, time_points = c("T0", "T2")) {
  current_eye <- patient_data$surgery_eye_1[1]
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  for(suffix in time_points) {
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Function to process all patients
process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# Process OCTA data
octa_bloodflow_processed <- process_all_patients(octa_bloodflow_features)
octa_thickness_processed <- process_all_patients(octa_thickness_features)

# -------------------- 3. Extract PPV group data --------------------
# Get PPV patient IDs
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

# Filter for PPV patients
ppv_bloodflow <- octa_bloodflow_processed %>%
  filter(ID %in% ppv_patients)

ppv_thickness <- octa_thickness_processed %>%
  filter(ID %in% ppv_patients)

# -------------------- 4. Process vision data --------------------
# Create vision dataset
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,
      surgery_eye_1 == 1 ~ os_corrected_bas,
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,
      TRUE ~ NA_real_
    ),
    post_vision_1m = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1m,
      surgery_eye_1 == 1 ~ os_corrected_1m,
      surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,
      TRUE ~ NA_real_
    ),
    vision_improvement_1m = post_vision_1m - pre_vision
  ) %>%
  dplyr::select(ID, vision_improvement_1m, pre_vision, post_vision_1m, age, gender)

# Filter vision data for PPV patients
ppv_vision <- vision_data %>%
  filter(ID %in% ppv_patients)

# *** CRITICAL: Convert vision data to LogMAR BEFORE further analysis ***
ppv_vision <- convert_vision_parameters(ppv_vision)

# -------------------- 5. Filter OCTA parameters --------------------
# Function to filter blood flow layers (enhanced for 0_21 only)
filter_bloodflow_layers <- function(data) {
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  # MODIFIED: Only include 0_21 region (removing 0_6)
  regions_of_interest <- c("0_21")
  
  # Create pattern for 0_21 region only
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target blood flow layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "target blood flow layer T0 parameters:\n")
    
    # Separate macular parameters for better reporting
    macular_params <- params_T0[grep("0_21_T0$", params_T0)]
    
    cat("- Macular region (0_21):", length(macular_params), "parameters\n")
    cat("- Total:", length(params_T0), "parameters\n\n")
    
    if(length(macular_params) > 0) {
      cat("Macular parameters:\n")
      cat(paste(gsub("_T0$", "", macular_params), collapse = "\n"), "\n\n")
    }
  }
  
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  if(length(params_T2) < length(params_T0)) {
    missing_params <- gsub("_T0$", "_T2", params_T0[!(gsub("_T0$", "_T2", params_T0) %in% params_T2)])
    warning("Missing T2 parameters: ", paste(missing_params, collapse = ", "))
  }
  
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params,
    macular_params = gsub("_T0$", "", macular_params),
    widefield_params = character(0)  # Empty since we removed 0_6
  ))
}

# Function to filter thickness layers (enhanced for 0_21 only)
filter_thickness_layers <- function(data) {
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  # MODIFIED: Only include 0_21 region (removing 0_6)
  regions_of_interest <- c("0_21")
  
  # Create pattern for 0_21 region only
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target thickness layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "target thickness layer T0 parameters:\n")
    
    # Separate macular parameters for better reporting
    macular_params <- params_T0[grep("0_21_T0$", params_T0)]
    
    cat("- Macular region (0_21):", length(macular_params), "parameters\n")
    cat("- Total:", length(params_T0), "parameters\n\n")
    
    if(length(macular_params) > 0) {
      cat("Macular thickness parameters:\n")
      cat(paste(gsub("_T0$", "", macular_params), collapse = "\n"), "\n\n")
    }
  }
  
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  if(length(params_T2) < length(params_T0)) {
    missing_params <- gsub("_T0$", "_T2", params_T0[!(gsub("_T0$", "_T2", params_T0) %in% params_T2)])
    warning("Missing T2 parameters: ", paste(missing_params, collapse = ", "))
  }
  
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params,
    macular_params = gsub("_T0$", "", macular_params),
    widefield_params = character(0)  # Empty since we removed 0_6
  ))
}

# Apply filters
ppv_bloodflow_filtered <- filter_bloodflow_layers(ppv_bloodflow)
ppv_thickness_filtered <- filter_thickness_layers(ppv_thickness)

# -------------------- 6. Calculate OCTA improvements --------------------
calculate_improvement <- function(data, params_T0, params_T2) {
  result <- data %>% dplyr::select(ID)
  
  for(i in 1:length(params_T0)) {
    t0_param <- params_T0[i]
    t2_param <- params_T2[i]
    base_param <- gsub("_T0$", "", t0_param)
    
    result[[paste0(base_param, "_improvement")]] <- data[[t2_param]] - data[[t0_param]]
  }
  
  return(result)
}

# Calculate improvements
ppv_bloodflow_improvement <- calculate_improvement(
  ppv_bloodflow_filtered$data, 
  ppv_bloodflow_filtered$params_T0, 
  ppv_bloodflow_filtered$params_T2
)

ppv_thickness_improvement <- calculate_improvement(
  ppv_thickness_filtered$data, 
  ppv_thickness_filtered$params_T0, 
  ppv_thickness_filtered$params_T2
)

# -------------------- 7. Combine OCTA and Vision data --------------------
# Merge blood flow and thickness improvements
ppv_octa_combined <- ppv_bloodflow_improvement %>%
  full_join(ppv_thickness_improvement, by = "ID", suffix = c("_bloodflow", "_thickness"))

# Merge with vision data (now in LogMAR)
comprehensive_data <- ppv_vision %>%
  dplyr::select(ID, vision_improvement_1m, pre_vision, age) %>%
  full_join(ppv_octa_combined, by = "ID")

# Print data combination summary (enhanced for regions)
cat("\n===== Comprehensive Data Summary (with LogMAR Vision) =====\n")
cat("Total patients in comprehensive dataset:", nrow(comprehensive_data), "\n")
cat("Vision parameters: vision_improvement_1m (LogMAR), pre_vision (LogMAR), age\n")
cat("Note: LogMAR conversion applied - positive improvement = vision gain\n")

if(length(ppv_bloodflow_filtered$macular_params) > 0) {
  cat("OCTA blood flow - Macular region (0_21):", length(ppv_bloodflow_filtered$macular_params), "parameters\n")
}
cat("Total OCTA blood flow parameters:", ncol(ppv_bloodflow_improvement) - 1, "\n")

if(length(ppv_thickness_filtered$macular_params) > 0) {
  cat("OCTA thickness - Macular region (0_21):", length(ppv_thickness_filtered$macular_params), "parameters\n")
}
cat("Total OCTA thickness parameters:", ncol(ppv_thickness_improvement) - 1, "\n")
cat("Total parameters:", ncol(comprehensive_data) - 1, "\n")

# -------------------- 8. Prepare data for clustering --------------------
# Analyze missing values
original_ids <- comprehensive_data$ID
clustering_data <- comprehensive_data %>% dplyr::select(-ID)

# Detailed missing value analysis
na_count_by_row <- rowSums(is.na(clustering_data))
na_count_table <- table(na_count_by_row)
na_count_by_col <- colSums(is.na(clustering_data))

cat("\n===== Missing Values Analysis =====\n")
cat("Total rows (patients):", nrow(clustering_data), "\n")
cat("Total columns (parameters):", ncol(clustering_data), "\n")

cat("\nPatient missing value distribution:\n")
for(na_count in names(na_count_table)) {
  cat("Patients with", na_count, "missing values:", na_count_table[na_count], "\n")
}

# Create summary of missing values by parameter type
na_summary <- data.frame(
  Parameter = names(clustering_data),
  Missing_Count = na_count_by_col,
  Missing_Percent = round(na_count_by_col/nrow(clustering_data)*100, 1),
  Data_Type = case_when(
    grepl("vision_improvement|pre_vision|age", names(clustering_data)) ~ "Vision",
    grepl("SVP|ICP|DCP|Choroid", names(clustering_data)) ~ "BloodFlow",
    grepl("GCL|INL|Retina", names(clustering_data)) ~ "Thickness",
    TRUE ~ "Other"
  )
)

cat("\nTop 15 parameters with most missing values:\n")
na_sorted <- na_summary %>% arrange(desc(Missing_Percent))
print(head(na_sorted, 15))

# Filter for complete cases
complete_rows <- complete.cases(clustering_data)
complete_data <- clustering_data[complete_rows, ]
complete_ids <- original_ids[complete_rows]

cat("\n===== Complete Data Statistics =====\n")
cat("Original patient count:", nrow(clustering_data), "\n")
cat("Complete data patient count:", nrow(complete_data), 
    "(", round(nrow(complete_data)/nrow(clustering_data)*100, 1), "%)\n")
cat("Removed", sum(!complete_rows), "patients with missing values\n")

# Count parameters by type
vision_params <- names(complete_data)[grepl("vision_improvement|pre_vision|age", names(complete_data))]
bloodflow_params <- names(complete_data)[grepl("SVP|ICP|DCP|Choroid", names(complete_data))]
thickness_params <- names(complete_data)[grepl("GCL|INL|Retina", names(complete_data))]

# Count parameters by type and region
bloodflow_macular <- bloodflow_params[grepl("0_21", bloodflow_params)]
bloodflow_widefield <- character(0)  # Empty since we removed 0_6
thickness_macular <- thickness_params[grepl("0_21", thickness_params)]
thickness_widefield <- character(0)  # Empty since we removed 0_6

cat("Vision parameters (LogMAR):", length(vision_params), "\n")
cat("Blood flow parameters - Total:", length(bloodflow_params), "\n")
cat("  - Macular (0_21):", length(bloodflow_macular), "\n")
cat("Thickness parameters - Total:", length(thickness_params), "\n")
cat("  - Macular (0_21):", length(thickness_macular), "\n")

# Check data adequacy
if(nrow(complete_data) < 5) {
  stop("Not enough patients with complete data. Consider imputation strategies.")
} else if(nrow(complete_data) < 10) {
  cat("⚠️ Warning: Small sample size (", nrow(complete_data), " patients)\n")
}

# Standardize data
comprehensive_data_std <- as.data.frame(scale(complete_data))

cat("\n===== Final Clustering Dataset =====\n")
cat("Total parameters used:", ncol(comprehensive_data_std), "\n")
cat("- Vision parameters (LogMAR):", length(vision_params), "\n")
cat("- Blood flow parameters:", length(bloodflow_params), "\n")
cat("  * Macular region:", length(bloodflow_macular), "\n")
cat("- Thickness parameters:", length(thickness_params), "\n")
cat("  * Macular region:", length(thickness_macular), "\n")
cat("Patients used:", nrow(comprehensive_data_std), "\n")
cat("Patient IDs:", paste(complete_ids, collapse=", "), "\n")

# -------------------- 9. Perform Mfuzz clustering --------------------
# Prepare ExpressionSet
prep_data_for_mfuzz <- function(data, row_ids) {
  data_matrix <- as.matrix(data)
  rownames(data_matrix) <- row_ids
  eset <- ExpressionSet(assayData = data_matrix)
  return(eset)
}

# Create ExpressionSet
comprehensive_eset <- prep_data_for_mfuzz(comprehensive_data_std, complete_ids)
comprehensive_eset_std <- standardise(comprehensive_eset)

# Estimate optimal fuzzy coefficient
comprehensive_m <- mestimate(comprehensive_eset_std)
cat("\nEstimated optimal fuzzy coefficient m:", comprehensive_m, "\n")

# Perform clustering with 2 clusters
set.seed(123)
comprehensive_clustering <- mfuzz(comprehensive_eset_std, c = 2, m = comprehensive_m)

# Get results
comprehensive_membership <- comprehensive_clustering$membership
comprehensive_main_clusters <- apply(comprehensive_membership, 1, which.max)

# Create results dataframe
comprehensive_clusters_result <- data.frame(
  subject_id = complete_ids,
  max_cluster = comprehensive_main_clusters,
  max_membership = apply(comprehensive_membership, 1, max)
)

# Save clustering results
write.csv(comprehensive_clusters_result, "ppv_comprehensive_cluster_results.csv", row.names = FALSE)

# -------------------- 10. Visualize clustering results (MODIFIED for LogMAR) --------------------
# Create comprehensive dataset with cluster assignments
comprehensive_data_with_clusters <- comprehensive_data %>%
  inner_join(comprehensive_clusters_result, by = c("ID" = "subject_id"))

# Function to create enhanced visualizations with LogMAR labels
create_comprehensive_visualizations <- function(data, vision_params, bloodflow_params, thickness_params, 
                                                bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield) {
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  data$max_cluster <- as.factor(data$max_cluster)
  all_params <- c(vision_params, bloodflow_params, thickness_params)
  
  # 0. Overall comprehensive heatmap for all parameters
  cat("Creating overall comprehensive heatmap (with LogMAR vision)...\n")
  
  # Calculate cluster means for ALL parameters
  all_params_means <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(all_params), mean, na.rm = TRUE), .groups = 'drop')
  
  # Prepare data for overall heatmap
  overall_plot_data <- all_params_means %>%
    pivot_longer(
      cols = all_of(all_params),
      names_to = "Parameter",
      values_to = "Mean_Value"
    ) %>%
    mutate(
      # Determine parameter category and region
      Data_Type = case_when(
        Parameter %in% vision_params ~ "Vision (LogMAR)",
        Parameter %in% bloodflow_macular ~ "Blood Flow - Macular",
        Parameter %in% thickness_macular ~ "Thickness - Macular",
        TRUE ~ "Other"
      ),
      # Clean parameter names for display
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter),
      Parameter_Clean = gsub("_0_21|_0_6", "", Parameter_Clean),
      Parameter_Clean = gsub("_", " ", Parameter_Clean),
      # Add region suffix for OCTA parameters
      Parameter_Display = case_when(
        grepl("0_21", Parameter) ~ paste0(Parameter_Clean, " (Mac)"),
        TRUE ~ Parameter_Clean
      )
    ) %>%
    # Order parameters by type and name
    arrange(Data_Type, Parameter_Display)
  
  # Create factor levels for proper ordering
  overall_plot_data$Parameter_Display <- factor(overall_plot_data$Parameter_Display, 
                                                levels = unique(overall_plot_data$Parameter_Display))
  
  # Prepare LogMAR-specific color mapping for different parameter types
  # For LogMAR vision parameters, we need different interpretation
  overall_plot_data <- overall_plot_data %>%
    mutate(
      # Create adjusted values for heatmap coloring
      Adjusted_Value = case_when(
        # For pre-vision (LogMAR): lower is better, so invert for color mapping
        Parameter == "pre_vision" ~ -Mean_Value,
        # For vision improvement (LogMAR): positive is better (keep as is)
        Parameter == "vision_improvement" ~ Mean_Value,
        # For age: neutral (center around mean)
        Parameter == "age" ~ Mean_Value - mean(Mean_Value, na.rm = TRUE),
        # For OCTA parameters: positive improvement is better (keep as is)
        TRUE ~ Mean_Value
      ),
      # Add interpretation labels
      Value_Interpretation = case_when(
        Parameter == "pre_vision" & Mean_Value < 0 ~ paste0(sprintf("%.3f", Mean_Value), "*"),
        Parameter == "pre_vision" & Mean_Value >= 0 ~ sprintf("%.3f", Mean_Value),
        Parameter == "vision_improvement" & Mean_Value > 0 ~ paste0(sprintf("%.3f", Mean_Value), "+"),
        Parameter == "vision_improvement" & Mean_Value <= 0 ~ sprintf("%.3f", Mean_Value),
        TRUE ~ sprintf("%.3f", Mean_Value)
      )
    )
  
  # Create the overall comprehensive heatmap with LogMAR-aware coloring
  p_overall_heatmap <- ggplot(overall_plot_data, aes(x = Parameter_Display, y = as.factor(max_cluster))) +
    geom_tile(aes(fill = Adjusted_Value), color = "white", size = 0.5) +
    # Add text labels with interpretation
    geom_text(aes(label = Value_Interpretation), 
              color = "black", size = 2.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#542788", 
      mid = "white", 
      high = "#f1a340", 
      midpoint = 0,
      name = "Adjusted\nValue*",
      guide = guide_colorbar(barwidth = 1, barheight = 10)
    ) +
    facet_wrap(~ Data_Type, scales = "free_x", ncol = 1, strip.position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      strip.background = element_rect(fill = "lightgray"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0, size = 10),
      legend.position = "right"
    ) +
    labs(
      title = "Comprehensive Cluster Comparison Heatmap (LogMAR-Adjusted)",
      subtitle = paste("Vision (LogMAR) + OCTA (Macular only) | n =", nrow(data)),
      x = "Parameters",
      y = "Cluster",
      caption = paste(
        "Green = Better outcomes, Red = Worse outcomes\n",
        "LogMAR Baseline: Lower values inverted for display (*)\n", 
        "LogMAR Improvement: Positive values marked (+)\n",
        "Mac = Macular (0_21)"
      )
    )
  
  # Save the overall heatmap
  ggsave("plots/overall_comprehensive_heatmap_logmar.pdf", p_overall_heatmap, 
         width = 20, height = 12, dpi = 300)
  ggsave("plots/overall_comprehensive_heatmap_logmar.png", p_overall_heatmap, 
         width = 20, height = 12, dpi = 300)
  
  cat("Overall comprehensive heatmap (LogMAR) saved to plots/overall_comprehensive_heatmap_logmar.pdf/png\n")
  
  # Additional summary statistics table
  cluster_summary_table <- overall_plot_data %>%
    group_by(max_cluster, Data_Type) %>%
    summarise(
      Parameter_Count = n(),
      Mean_Improvement = mean(Mean_Value, na.rm = TRUE),
      Positive_Changes = sum(Mean_Value > 0, na.rm = TRUE),
      Negative_Changes = sum(Mean_Value < 0, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(max_cluster, Data_Type)
  
  cat("\n===== Cluster Summary by Parameter Type (LogMAR) =====\n")
  print(cluster_summary_table)
  
  # Save summary table
  write.csv(cluster_summary_table, "plots/cluster_summary_by_type_logmar.csv", row.names = FALSE)
  
  # 1. Vision parameters by cluster (MODIFIED for LogMAR)
  for(param in vision_params) {
    if(param %in% names(data)) {
      # Get proper labels for the parameter
      plot_labels <- get_vision_plot_labels(param)
      
      p <- ggplot(data, aes(x = max_cluster, y = .data[[param]], fill = max_cluster)) +
        geom_boxplot() +
        geom_jitter(alpha = 0.5, width = 0.2) +
        scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
        theme_bw() +
        labs(
          title = plot_labels$title, 
          subtitle = plot_labels$subtitle,
          x = "Cluster", 
          y = plot_labels$y_label
        )
      
      param_clean <- gsub("_improvement|_1m", "", param)
      ggsave(paste0("plots/vision_", param_clean, "_logmar_boxplot.pdf"), p, width = 8, height = 6)
      ggsave(paste0("plots/vision_", param_clean, "_logmar_boxplot.png"), p, width = 8, height = 6, dpi = 300)
    }
  }
  
  # 1.2. Age boxplot by cluster
  if("age" %in% names(data)) {
    age_labels <- get_vision_plot_labels("age")
    p_age <- ggplot(data, aes(x = factor(max_cluster), y = age, fill = factor(max_cluster))) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(aes(color = factor(max_cluster)), width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster") +
      scale_color_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster") +
      theme_bw() +
      labs(
        title = age_labels$title,
        subtitle = age_labels$subtitle,
        x = "Cluster", 
        y = age_labels$y_label
      )
    
    ggsave("plots/age_boxplot.pdf", p_age, width = 8, height = 6)
    ggsave("plots/age_boxplot.png", p_age, width = 8, height = 6, dpi = 300)
  }
  
  # 1.3. Pre-vision boxplot by cluster (LogMAR)
  if("pre_vision" %in% names(data)) {
    pre_vision_labels <- get_vision_plot_labels("pre_vision")
    p_pre_vision <- ggplot(data, aes(x = factor(max_cluster), y = pre_vision, fill = factor(max_cluster))) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(aes(color = factor(max_cluster)), width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster") +
      scale_color_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster") +
      theme_bw() +
      labs(
        title = pre_vision_labels$title,
        subtitle = pre_vision_labels$subtitle,
        x = "Cluster", 
        y = pre_vision_labels$y_label
      )
    
    ggsave("plots/pre_vision_logmar_boxplot.pdf", p_pre_vision, width = 8, height = 6)
    ggsave("plots/pre_vision_logmar_boxplot.png", p_pre_vision, width = 8, height = 6, dpi = 300)
  }
  
  # 1.4. Vision improvement boxplot by cluster (LogMAR) 
  if("vision_improvement_1m" %in% names(data)) {
    improvement_labels <- get_vision_plot_labels("vision_improvement_1m")
    p_improvement <- ggplot(data, aes(x = factor(max_cluster), y = vision_improvement_1m, fill = factor(max_cluster))) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(aes(color = factor(max_cluster)), width = 0.2, alpha = 0.6, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
      scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster") +
      scale_color_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster") +
      theme_bw() +
      labs(
        title = improvement_labels$title,
        subtitle = improvement_labels$subtitle,
        x = "Cluster", 
        y = improvement_labels$y_label,
        caption = "Red dashed line = no change"
      )
    
    ggsave("plots/vision_improvement_logmar_boxplot.pdf", p_improvement, width = 8, height = 6)
    ggsave("plots/vision_improvement_logmar_boxplot.png", p_improvement, width = 8, height = 6, dpi = 300)
  }
  
  # 1.5. OCTA specific boxplots for PA, VD, and Thickness indicators
  cat("Creating OCTA-specific boxplots for PA, VD, and Thickness indicators...\n")
  
  # Function to create OCTA indicator boxplots
  create_octa_boxplots <- function(data, all_params) {
    # Define OCTA indicators of interest
    pa_params <- all_params[grepl("PA.*_improvement", all_params)]
    vd_params <- all_params[grepl("VD.*_improvement", all_params)]
    thickness_params <- all_params[grepl("(GCL|INL|Retina).*_improvement", all_params)]
    
    # Function to create boxplot for a group of parameters
    create_grouped_boxplot <- function(params, indicator_name, color_palette = "Set2") {
      if(length(params) == 0) {
        cat("No", indicator_name, "parameters found\n")
        return(NULL)
      }
      
      cat("Creating", indicator_name, "boxplot with", length(params), "parameters\n")
      
      # Prepare data for plotting
      plot_data <- data %>%
        dplyr::select(max_cluster, all_of(params)) %>%
        pivot_longer(cols = all_of(params), names_to = "Parameter", values_to = "Value") %>%
        mutate(
          Parameter_Clean = gsub("_improvement", "", Parameter),
          Parameter_Clean = gsub("_0_21", "", Parameter_Clean),
          Parameter_Clean = gsub("_", " ", Parameter_Clean)
        )
      
      # Create boxplot
      p <- ggplot(plot_data, aes(x = Parameter_Clean, y = Value, fill = max_cluster)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
        geom_jitter(aes(color = max_cluster), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
                    alpha = 0.6, size = 1) +
        scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster")+
        scale_color_manual(values = c("#1f77b4", "#ff7f0e"), name = "Cluster")+
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "top"
        ) +
        labs(
          title = paste(indicator_name, "Improvements by Cluster"),
          subtitle = paste("Comparison across", length(params), "parameters"),
          x = paste(indicator_name, "Parameters"),
          y = "Improvement Value"
        )
      
      return(p)
    }
    
    # Create PA boxplot
    if(length(pa_params) > 0) {
      p_pa <- create_grouped_boxplot(pa_params, "Perfusion Area (PA)", "Set1")
      if(!is.null(p_pa)) {
        ggsave("plots/octa_PA_improvements_boxplot.pdf", p_pa, width = 12, height = 8)
        ggsave("plots/octa_PA_improvements_boxplot.png", p_pa, width = 12, height = 8, dpi = 300)
        cat("PA boxplot saved to plots/octa_PA_improvements_boxplot.pdf/png\n")
      }
    }
    
    # Create VD boxplot
    if(length(vd_params) > 0) {
      p_vd <- create_grouped_boxplot(vd_params, "Vessel Density (VD)", "Set1")
      if(!is.null(p_vd)) {
        ggsave("plots/octa_VD_improvements_boxplot.pdf", p_vd, width = 12, height = 8)
        ggsave("plots/octa_VD_improvements_boxplot.png", p_vd, width = 12, height = 8, dpi = 300)
        cat("VD boxplot saved to plots/octa_VD_improvements_boxplot.pdf/png\n")
      }
    }
    
    # Create Thickness boxplot
    if(length(thickness_params) > 0) {
      p_thickness <- create_grouped_boxplot(thickness_params, "Retinal Thickness", "Set1")
      if(!is.null(p_thickness)) {
        ggsave("plots/octa_Thickness_improvements_boxplot.pdf", p_thickness, width = 12, height = 8)
        ggsave("plots/octa_Thickness_improvements_boxplot.png", p_thickness, width = 12, height = 8, dpi = 300)
        cat("Thickness boxplot saved to plots/octa_Thickness_improvements_boxplot.pdf/png\n")
      }
    }
    
    # Create combined summary boxplot for all three indicators
    combined_params <- c(pa_params, vd_params, thickness_params)
    if(length(combined_params) > 0) {
      combined_plot_data <- data %>%
        dplyr::select(max_cluster, all_of(combined_params)) %>%
        pivot_longer(cols = all_of(combined_params), names_to = "Parameter", values_to = "Value") %>%
        mutate(
          Indicator_Type = case_when(
            grepl("PA", Parameter) ~ "Perfusion Area (PA)",
            grepl("VD", Parameter) ~ "Vessel Density (VD)",
            grepl("GCL|INL|Retina", Parameter) ~ "Retinal Thickness",
            TRUE ~ "Other"
          ),
          Parameter_Clean = gsub("_improvement", "", Parameter),
          Parameter_Clean = gsub("_0_21", "", Parameter_Clean),
          Parameter_Clean = gsub("_", " ", Parameter_Clean)
        )
      
      p_combined <- ggplot(combined_plot_data, aes(x = Indicator_Type, y = Value, fill = max_cluster)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
        geom_jitter(aes(color = max_cluster), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), 
                    alpha = 0.5, size = 1.5) +
        scale_fill_manual(values = c("#2ca02c", "#d62728"), name = "Cluster")+
        scale_color_manual(values = c("#2ca02c", "#d62728"), name = "Cluster")+
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 15, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "top",
          strip.text = element_text(size = 12, face = "bold")
        ) +
        labs(
          title = "OCTA Indicators Improvement Comparison by Cluster",
          subtitle = paste("PA, VD, and Thickness improvements | n =", nrow(data)),
          x = "OCTA Indicator Type",
          y = "Improvement Value"
        ) +
        facet_wrap(~ Indicator_Type, scales = "free_x", ncol = 3)
      
      ggsave("plots/octa_combined_indicators_boxplot.pdf", p_combined, width = 16, height = 10)
      ggsave("plots/octa_combined_indicators_boxplot.png", p_combined, width = 16, height = 10, dpi = 300)
      cat("Combined OCTA indicators boxplot saved to plots/octa_combined_indicators_boxplot.pdf/png\n")
    }
    
    # Summary statistics for OCTA indicators by cluster
    if(length(combined_params) > 0) {
      octa_summary <- combined_plot_data %>%
        group_by(max_cluster, Indicator_Type) %>%
        summarise(
          Parameter_Count = n(),
          Mean_Improvement = mean(Value, na.rm = TRUE),
          Median_Improvement = median(Value, na.rm = TRUE),
          SD_Improvement = sd(Value, na.rm = TRUE),
          Positive_Improvements = sum(Value > 0, na.rm = TRUE),
          Negative_Improvements = sum(Value < 0, na.rm = TRUE),
          .groups = 'drop'
        )
      
      cat("\n===== OCTA Indicators Summary by Cluster =====\n")
      print(octa_summary)
      
      # Save OCTA summary
      write.csv(octa_summary, "plots/octa_indicators_cluster_summary.csv", row.names = FALSE)
      cat("OCTA indicators summary saved to plots/octa_indicators_cluster_summary.csv\n")
    }
  }
  
  # Execute OCTA boxplot creation
  create_octa_boxplots(data, all_params)
  
  # 2. Enhanced heatmap with regional breakdown (LogMAR-aware coloring)
  cluster_means <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(all_params), mean, na.rm = TRUE))
  
  plot_data <- cluster_means %>%
    pivot_longer(
      cols = all_of(all_params),
      names_to = "Parameter",
      values_to = "Mean_Value"
    ) %>%
    mutate(
      Data_Type = case_when(
        Parameter %in% vision_params ~ "Vision (LogMAR)",
        Parameter %in% bloodflow_macular ~ "Blood Flow - Macular",
        Parameter %in% thickness_macular ~ "Thickness - Macular",
        TRUE ~ "Other"
      ),
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter),
      # LogMAR-specific adjustments for heatmap display
      Adjusted_Value = case_when(
        # For pre-vision (LogMAR): lower is better, so invert for color mapping
        Parameter == "pre_vision" ~ -Mean_Value,
        # For vision improvement (LogMAR): positive is better (keep as is)  
        Parameter == "vision_improvement_1m" ~ Mean_Value,
        # For age: center around overall mean
        Parameter == "age" ~ Mean_Value - mean(Mean_Value, na.rm = TRUE),
        # For OCTA parameters: positive improvement is better (keep as is)
        TRUE ~ Mean_Value
      ),
      # Create display labels with LogMAR interpretation
      Display_Value = case_when(
        Parameter == "pre_vision" ~ paste0(sprintf("%.3f", Mean_Value), "*"),
        Parameter == "vision_improvement_1m" & Mean_Value > 0 ~ paste0(sprintf("%.3f", Mean_Value), "+"),
        TRUE ~ sprintf("%.3f", Mean_Value)
      )
    )
  
  # Create comprehensive heatmap with LogMAR-aware coloring
  p_heatmap <- ggplot(plot_data, aes(x = Parameter_Clean, y = max_cluster)) +
    geom_tile(aes(fill = Adjusted_Value), color = "white", size = 0.5) +
    # Add text with LogMAR interpretation
    geom_text(aes(label = Display_Value), color = "black", size = 2.8, fontface = "bold") +
    scale_fill_gradient2(
      low = "#542788", 
      mid = "white", 
      high = "#f1a340", 
      midpoint = 0,
      name = "Adjusted\nValue"
    ) +
    facet_wrap(~ Data_Type, scales = "free_x", ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      strip.text = element_text(size = 11, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9)
    ) +
    labs(
      title = "Comprehensive Parameter Overview by Cluster\n(LogMAR-Adjusted Vision + OCTA Macular)",
      subtitle = "Color-coded for clinical interpretation: Green = Better outcomes",
      x = "Parameters",
      y = "Cluster", 
      fill = "Adjusted\nValue",
      caption = paste(
        "LogMAR Baseline (*): Lower values inverted for display\n",
        "LogMAR Improvement (+): Positive values indicate improvement\n",
        "Green = Better clinical outcomes, Red = Worse outcomes"
      )
    )
  
  ggsave("plots/comprehensive_regional_heatmap_logmar.pdf", p_heatmap, width = 24, height = 12)
  ggsave("plots/comprehensive_regional_heatmap_logmar.png", p_heatmap, width = 24, height = 12, dpi = 300)
  
  # 3. PCA visualization
  if(length(all_params) > 2) {
    pca_data <- data %>%
      dplyr::select(all_of(all_params))
    
    pca_result <- prcomp(pca_data, scale. = TRUE)
    
    # Create PCA plot data frame - THIS IS THE CORRECT SYNTAX
    pca_plot_data <- data.frame(
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2],
      Cluster = data$max_cluster,
      Membership = data$max_membership
    )
    
    # Create PCA plot
    p_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Cluster, alpha = Membership)) +
      geom_point(size = 3) +
      stat_ellipse(aes(group = Cluster), level = 0.95) +
      scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
      theme_bw() +
      labs(
        title = "PCA of Comprehensive Parameters\n(Vision LogMAR + OCTA Macular)",
        subtitle = "LogMAR vision parameters included",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    ggsave("plots/comprehensive_regional_pca_logmar.pdf", p_pca, width = 10, height = 8)
    ggsave("plots/comprehensive_regional_pca_logmar.png", p_pca, width = 10, height = 8, dpi = 300)
    
    # 4. Enhanced variable contribution plot with regional colors
    loadings <- pca_result$rotation[, 1:2]
    loadings_df <- data.frame(
      Variable = rownames(loadings),
      PC1 = loadings[, 1],
      PC2 = loadings[, 2],
      Data_Type = case_when(
        rownames(loadings) %in% vision_params ~ "Vision (LogMAR)",
        rownames(loadings) %in% bloodflow_macular ~ "BF - Macular",
        rownames(loadings) %in% thickness_macular ~ "TH - Macular",
        TRUE ~ "Other"
      )
    )
    
    p_loadings <- ggplot(loadings_df, aes(x = PC1, y = PC2, color = Data_Type)) +
      geom_point(size = 3) +
      geom_text(aes(label = gsub("_improvement|_1m", "", Variable)), 
                vjust = -0.5, hjust = 0.5, size = 2.5) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      scale_color_brewer(palette = "Set2") +
      theme_bw() +
      labs(
        title = "Parameter Contributions to Principal Components\n(LogMAR Vision + Regional OCTA Analysis)",
        subtitle = "LogMAR vision parameters included in analysis",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
        color = "Parameter Type"
      )
    
    ggsave("plots/comprehensive_regional_loadings_logmar.pdf", p_loadings, width = 14, height = 12)
    ggsave("plots/comprehensive_regional_loadings_logmar.png", p_loadings, width = 14, height = 12, dpi = 300)
  }
}

# Create visualizations (with LogMAR labels)
create_comprehensive_visualizations(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                                    bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield)

# -------------------- 11. Statistical analysis --------------------
# Function to analyze differences between clusters (enhanced for LogMAR vision)
analyze_comprehensive_differences <- function(data, vision_params, bloodflow_params, thickness_params,
                                              bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield) {
  all_params <- c(vision_params, bloodflow_params, thickness_params)
  
  results <- data.frame(
    Parameter = character(),
    Data_Type = character(),
    Region = character(),
    Cluster1_Mean = numeric(),
    Cluster2_Mean = numeric(),
    Mean_Difference = numeric(),
    P_Value = numeric(),
    Significant = character(),
    Parameter_Type = character(),  # Add parameter type for LogMAR interpretation
    stringsAsFactors = FALSE
  )
  
  for(param in all_params) {
    if(param %in% names(data) && !all(is.na(data[[param]]))) {
      # Determine data type and region
      data_type <- case_when(
        param %in% vision_params ~ "Vision (LogMAR)",
        param %in% bloodflow_params ~ "Blood Flow",
        param %in% thickness_params ~ "Thickness",
        TRUE ~ "Other"
      )
      
      region <- case_when(
        param %in% vision_params ~ "N/A",
        param %in% bloodflow_macular ~ "Macular",
        param %in% thickness_macular ~ "Macular",
        TRUE ~ "Unknown"
      )
      
      # Determine parameter type for interpretation
      param_type <- case_when(
        param == "pre_vision" ~ "Baseline (Lower=Better)",
        param == "vision_improvement_1m" ~ "Improvement (Positive=Better)",
        param %in% vision_params ~ "Vision (LogMAR)",
        TRUE ~ "OCTA"
      )
      
      param_data <- data[, c("max_cluster", param)]
      param_data <- param_data[!is.na(param_data[[param]]), ]
      
      if(nrow(param_data) > 0) {
        means <- tapply(param_data[[param]], param_data$max_cluster, mean, na.rm = TRUE)
        
        test_result <- try(t.test(reformulate("max_cluster", param), data = param_data), silent = TRUE)
        
        if(class(test_result) != "try-error") {
          results <- rbind(results, data.frame(
            Parameter = gsub("_improvement|_1m", "", param),
            Data_Type = data_type,
            Region = region,
            Cluster1_Mean = ifelse("1" %in% names(means), means["1"], NA),
            Cluster2_Mean = ifelse("2" %in% names(means), means["2"], NA),
            Mean_Difference = ifelse(length(means) >= 2, means["2"] - means["1"], NA),
            P_Value = test_result$p.value,
            Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
            Parameter_Type = param_type,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Add adjusted p-values
  if(nrow(results) > 0) {
    results$P_Adjusted <- p.adjust(results$P_Value, method = "fdr")
    results$Significant_Adjusted <- ifelse(results$P_Adjusted < 0.05, "Yes", "No")
    results <- results %>% arrange(Data_Type, Region, P_Value)
  }
  
  return(results)
}

# Perform statistical analysis (updated function call for LogMAR)
comprehensive_stats <- analyze_comprehensive_differences(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                                                         bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield)

# Save results
write.csv(comprehensive_stats, "ppv_comprehensive_cluster_statistics_logmar.csv", row.names = FALSE)

# Print significant results by data type and region (with LogMAR interpretation)
cat("\n===== Significant Parameters by Data Type and Region (LogMAR Vision) =====\n")

for(data_type in c("Vision (LogMAR)", "Blood Flow", "Thickness")) {
  cat("\n", data_type, "parameters:\n")
  
  if(data_type == "Vision (LogMAR)") {
    sig_params <- comprehensive_stats %>% 
      filter(Data_Type == data_type & Significant == "Yes") %>%
      arrange(P_Value)
    
    if(nrow(sig_params) > 0) {
      cat("LogMAR Interpretation Guide:\n")
      cat("- Pre-vision (LogMAR): Lower values = better vision\n")
      cat("- Vision improvement (LogMAR): Positive values = improvement\n")
      cat("- Mean Difference: Cluster2 - Cluster1\n\n")
      print(sig_params %>% dplyr::select(Parameter, Parameter_Type, Mean_Difference, P_Value, Significant))
    } else {
      cat("No significant parameters found\n")
    }
  } else {
    # For OCTA parameters, show by region
    for(region in c("Macular")) {  # Only Macular since we removed 0_6
      sig_params <- comprehensive_stats %>% 
        filter(Data_Type == data_type & Region == region & Significant == "Yes") %>%
        arrange(P_Value)
      
      if(nrow(sig_params) > 0) {
        cat("  ", region, "region:\n")
        print(sig_params %>% dplyr::select(Parameter, Mean_Difference, P_Value, Significant))
      } else {
        cat("  ", region, "region: No significant parameters found\n")
      }
    }
  }
}

# -------------------- 12. Interpret clusters (with LogMAR considerations) --------------------
# Comprehensive cluster interpretation function (modified for LogMAR)
interpret_comprehensive_clusters <- function(stats_results) {
  if(nrow(stats_results) == 0) {
    cat("No statistical results available for interpretation\n")
    return(NULL)
  }
  
  # Analyze by data type and region
  vision_stats <- stats_results %>% filter(Data_Type == "Vision (LogMAR)")
  bloodflow_macular_stats <- stats_results %>% filter(Data_Type == "Blood Flow" & Region == "Macular")
  thickness_macular_stats <- stats_results %>% filter(Data_Type == "Thickness" & Region == "Macular")
  
  # Count significant improvements by cluster and region
  # For LogMAR vision improvement: positive difference = cluster 2 better
  # For LogMAR pre-vision: negative difference = cluster 2 better (lower LogMAR = better vision)
  vision_cluster2_better <- 0
  vision_cluster1_better <- 0
  
  if(nrow(vision_stats) > 0) {
    for(i in 1:nrow(vision_stats)) {
      param <- vision_stats$Parameter[i]
      mean_diff <- vision_stats$Mean_Difference[i]
      is_significant <- vision_stats$Significant[i] == "Yes"
      
      if(is_significant) {
        if(param == "pre_vision") {
          # For pre-vision (LogMAR), negative difference means cluster 2 has better vision
          if(mean_diff < 0) vision_cluster2_better <- vision_cluster2_better + 1
          else vision_cluster1_better <- vision_cluster1_better + 1
        } else if(param == "vision_improvement") {
          # For vision improvement (LogMAR), positive difference means cluster 2 has better improvement
          if(mean_diff > 0) vision_cluster2_better <- vision_cluster2_better + 1
          else vision_cluster1_better <- vision_cluster1_better + 1
        }
      }
    }
  }
  
  bf_mac_cluster2_better <- sum(bloodflow_macular_stats$Mean_Difference > 0 & bloodflow_macular_stats$Significant == "Yes", na.rm = TRUE)
  bf_mac_cluster1_better <- sum(bloodflow_macular_stats$Mean_Difference < 0 & bloodflow_macular_stats$Significant == "Yes", na.rm = TRUE)
  
  th_mac_cluster2_better <- sum(thickness_macular_stats$Mean_Difference > 0 & thickness_macular_stats$Significant == "Yes", na.rm = TRUE)
  th_mac_cluster1_better <- sum(thickness_macular_stats$Mean_Difference < 0 & thickness_macular_stats$Significant == "Yes", na.rm = TRUE)
  
  # Total significant improvements
  total_cluster2_better <- vision_cluster2_better + bf_mac_cluster2_better + th_mac_cluster2_better
  total_cluster1_better <- vision_cluster1_better + bf_mac_cluster1_better + th_mac_cluster1_better
  
  # Determine better cluster
  if(total_cluster2_better > total_cluster1_better) {
    better_cluster <- 2
    worse_cluster <- 1
  } else if(total_cluster1_better > total_cluster2_better) {
    better_cluster <- 1
    worse_cluster <- 2
  } else {
    # If tied, prioritize vision improvement
    if(vision_cluster2_better > vision_cluster1_better) {
      better_cluster <- 2
      worse_cluster <- 1
    } else if(vision_cluster1_better > vision_cluster2_better) {
      better_cluster <- 1
      worse_cluster <- 2
    } else {
      # If still tied, use total number of significant parameters
      total_sig_cluster2 <- sum(stats_results$Mean_Difference > 0 & stats_results$Significant == "Yes", na.rm = TRUE)
      total_sig_cluster1 <- sum(stats_results$Mean_Difference < 0 & stats_results$Significant == "Yes", na.rm = TRUE)
      
      if(total_sig_cluster2 >= total_sig_cluster1) {
        better_cluster <- 2
        worse_cluster <- 1
      } else {
        better_cluster <- 1
        worse_cluster <- 2
      }
    }
  }
  
  # Print interpretation results (with LogMAR considerations)
  cat("\n===== Comprehensive Cluster Interpretation (LogMAR Vision) =====\n")
  cat("Overall better cluster:", better_cluster, "\n")
  cat("Overall worse cluster:", worse_cluster, "\n\n")
  
  cat("Detailed breakdown (considering LogMAR interpretation):\n")
  cat("Vision (LogMAR) - Cluster 2 advantages:", vision_cluster2_better, ", Cluster 1 advantages:", vision_cluster1_better, "\n")
  cat("  Note: For LogMAR, lower baseline = better, positive improvement = better\n")
  cat("Blood Flow Macular - Cluster 2 advantages:", bf_mac_cluster2_better, ", Cluster 1 advantages:", bf_mac_cluster1_better, "\n")
  cat("Thickness Macular - Cluster 2 advantages:", th_mac_cluster2_better, ", Cluster 1 advantages:", th_mac_cluster1_better, "\n")
  cat("Total - Cluster 2 advantages:", total_cluster2_better, ", Cluster 1 advantages:", total_cluster1_better, "\n\n")
  
  # Additional LogMAR-specific interpretation
  if(nrow(vision_stats) > 0) {
    cat("LogMAR Vision Parameter Details:\n")
    for(i in 1:nrow(vision_stats)) {
      param <- vision_stats$Parameter[i]
      mean_diff <- vision_stats$Mean_Difference[i]
      p_val <- vision_stats$P_Value[i]
      is_sig <- vision_stats$Significant[i] == "Yes"
      
      if(is_sig) {
        cat("- ", param, ": Mean difference =", round(mean_diff, 3))
        if(param == "pre_vision") {
          if(mean_diff < 0) {
            cat(" (Cluster 2 has better baseline vision)")
          } else {
            cat(" (Cluster 1 has better baseline vision)")
          }
        } else if(param == "vision_improvement") {
          if(mean_diff > 0) {
            cat(" (Cluster 2 shows greater improvement)")
          } else {
            cat(" (Cluster 1 shows greater improvement)")
          }
        }
        cat(" | p =", round(p_val, 3), "\n")
      }
    }
    cat("\n")
  }
  
  # Additional summary
  cat("Summary:\n")
  if(better_cluster == 2) {
    cat("Cluster 2 patients show significantly better outcomes across multiple parameters\n")
  } else {
    cat("Cluster 1 patients show significantly better outcomes across multiple parameters\n")
  }
  
  cat("This indicates that PPV patients can be stratified into distinct outcome groups\n")
  cat("based on their comprehensive vision (LogMAR) and OCTA response patterns.\n")
  
  # Return results
  return(list(
    better_cluster = better_cluster,
    worse_cluster = worse_cluster,
    vision_cluster2_better = vision_cluster2_better,
    vision_cluster1_better = vision_cluster1_better,
    bf_mac_cluster2_better = bf_mac_cluster2_better,
    bf_mac_cluster1_better = bf_mac_cluster1_better,
    th_mac_cluster2_better = th_mac_cluster2_better,
    th_mac_cluster1_better = th_mac_cluster1_better,
    total_cluster2_better = total_cluster2_better,
    total_cluster1_better = total_cluster1_better
  ))
}

# Execute cluster interpretation (with LogMAR considerations)
cluster_interpretation <- interpret_comprehensive_clusters(comprehensive_stats)

# Add outcome labels to cluster results
if(!is.null(cluster_interpretation)) {
  comprehensive_clusters_result$outcome_quality <- ifelse(
    comprehensive_clusters_result$max_cluster == cluster_interpretation$better_cluster,
    "Better",
    "Worse"
  )
  
  # Save final results with outcome labels
  write.csv(comprehensive_clusters_result, "ppv_comprehensive_cluster_results_with_outcomes_logmar.csv", row.names = FALSE)
  
  cat("\nCluster outcome labels have been added and saved to:\n")
  cat("'ppv_comprehensive_cluster_results_with_outcomes_logmar.csv'\n")
} else {
  cat("\nWarning: Could not determine cluster interpretation due to insufficient statistical results.\n")
}

# -------------------- 13. Parameter importance analysis (LogMAR version) --------------------
# Analyze which parameters contribute most to cluster separation (with LogMAR consideration)
analyze_parameter_importance <- function(stats_results) {
  if(nrow(stats_results) == 0) {
    return(data.frame())
  }
  
  importance_scores <- stats_results %>%
    mutate(
      Effect_Size = abs(Mean_Difference),
      Importance_Score = Effect_Size * (-log10(P_Value + 1e-10)),
      # Add LogMAR interpretation
      Clinical_Interpretation = case_when(
        Parameter == "pre_vision" & Mean_Difference < 0 ~ "Cluster 2 has better baseline vision",
        Parameter == "pre_vision" & Mean_Difference > 0 ~ "Cluster 1 has better baseline vision",
        Parameter == "vision_improvement" & Mean_Difference > 0 ~ "Cluster 2 shows greater improvement",
        Parameter == "vision_improvement" & Mean_Difference < 0 ~ "Cluster 1 shows greater improvement",
        Data_Type == "Vision (LogMAR)" ~ "LogMAR parameter",
        Mean_Difference > 0 ~ "Cluster 2 higher values",
        Mean_Difference < 0 ~ "Cluster 1 higher values",
        TRUE ~ "No clear pattern"
      )
    ) %>%
    arrange(desc(Importance_Score))
  
  # Get top parameters by data type
  top_vision <- importance_scores %>%
    filter(Data_Type == "Vision (LogMAR)") %>%
    head(3)
  
  top_bloodflow <- importance_scores %>%
    filter(Data_Type == "Blood Flow") %>%
    head(5)
  
  top_thickness <- importance_scores %>%
    filter(Data_Type == "Thickness") %>%
    head(5)
  
  cat("\n===== Parameter Importance Analysis (LogMAR Vision) =====\n")
  
  if(nrow(top_vision) > 0) {
    cat("\nTop Vision parameters (LogMAR):\n")
    print(top_vision %>% dplyr::select(Parameter, Parameter_Type, Mean_Difference, P_Value, Importance_Score, Clinical_Interpretation, Significant))
  }
  
  if(nrow(top_bloodflow) > 0) {
    cat("\nTop Blood Flow parameters:\n")
    print(top_bloodflow %>% dplyr::select(Parameter, Data_Type, Mean_Difference, P_Value, Importance_Score, Significant))
  }
  
  if(nrow(top_thickness) > 0) {
    cat("\nTop Thickness parameters:\n")
    print(top_thickness %>% dplyr::select(Parameter, Data_Type, Mean_Difference, P_Value, Importance_Score, Significant))
  }
  
  # Save importance analysis
  write.csv(importance_scores, "ppv_comprehensive_parameter_importance_logmar.csv", row.names = FALSE)
  
  return(importance_scores)
}

# Perform parameter importance analysis (LogMAR version)
parameter_importance <- analyze_parameter_importance(comprehensive_stats)

# -------------------- 14. Create comprehensive summary (LogMAR version) --------------------
# Final summary statistics (with LogMAR interpretation)
final_summary <- comprehensive_data_with_clusters %>%
  group_by(max_cluster) %>%
  summarise(
    Count = n(),
    Mean_Membership = mean(max_membership, na.rm = TRUE),
    # Vision parameters (LogMAR)
    Mean_Vision_Improvement_LogMAR = mean(vision_improvement_1m, na.rm = TRUE),
    Mean_Pre_Vision_LogMAR = mean(pre_vision, na.rm = TRUE),
    Mean_Age = mean(age, na.rm = TRUE),
    # Summary of OCTA improvements (count of positive improvements)
    Positive_BF_Improvements = rowSums(across(all_of(bloodflow_params), ~ .x > 0), na.rm = TRUE),
    Positive_TH_Improvements = rowSums(across(all_of(thickness_params), ~ .x > 0), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Outcome_Quality = ifelse(max_cluster == cluster_interpretation$better_cluster, "Better", "Worse"),
    # Add LogMAR interpretation columns
    Vision_Quality_Baseline = case_when(
      Mean_Pre_Vision_LogMAR < 0.3 ~ "Good (LogMAR < 0.3)",
      Mean_Pre_Vision_LogMAR < 0.7 ~ "Moderate (LogMAR 0.3-0.7)",
      TRUE ~ "Poor (LogMAR > 0.7)"
    ),
    Vision_Improvement_Category = case_when(
      Mean_Vision_Improvement_LogMAR > 0.1 ~ "Significant Improvement (>0.1 LogMAR)",
      Mean_Vision_Improvement_LogMAR > 0.0 ~ "Mild Improvement (0-0.1 LogMAR)",
      TRUE ~ "No/Negative Change"
    )
  )

cat("\n===== Final Comprehensive Summary (LogMAR Vision) =====\n")
print(final_summary)

write.csv(final_summary, "ppv_comprehensive_cluster_summary_logmar.csv", row.names = FALSE)

# -------------------- 15. Generate comprehensive report (LogMAR version) --------------------
generate_comprehensive_report <- function() {
  vision_sig <- sum(comprehensive_stats$Data_Type == "Vision (LogMAR)" & comprehensive_stats$Significant == "Yes")
  bf_mac_sig <- sum(comprehensive_stats$Data_Type == "Blood Flow" & comprehensive_stats$Region == "Macular" & comprehensive_stats$Significant == "Yes")
  th_mac_sig <- sum(comprehensive_stats$Data_Type == "Thickness" & comprehensive_stats$Region == "Macular" & comprehensive_stats$Significant == "Yes")
  total_sig <- sum(comprehensive_stats$Significant == "Yes")
  
  report <- paste0(
    "========================================\n",
    "PPV Group Comprehensive Clustering Analysis Report\n",
    "Vision (LogMAR) + OCTA (Macular Region 0_21 Only)\n",
    "========================================\n\n",
    
    "⭐ KEY MODIFICATION: Vision parameters converted to LogMAR before analysis\n",
    "   - Pre-vision: LogMAR scale (lower = better vision)\n",
    "   - Vision improvement: LogMAR change (positive = improvement)\n",
    "   - Age: unchanged (years)\n\n",
    
    "1. Data Overview:\n",
    "   - Total patients analyzed: ", nrow(complete_data), "\n",
    "   - Vision parameters (LogMAR): ", length(vision_params), " (improvement, baseline vision, age)\n",
    "   - OCTA blood flow - Macular: ", length(bloodflow_macular), " parameters\n",
    "   - OCTA thickness - Macular: ", length(thickness_macular), " parameters\n",
    "   - Total parameters: ", ncol(complete_data), "\n",
    "   - Focus: Central retinal region (0_21) analysis with LogMAR vision\n\n",
    
    "2. LogMAR Conversion Details:\n",
    "   - Formula: LogMAR = -log10(decimal_vision)\n",
    "   - Interpretation: Lower LogMAR = better vision\n",
    "   - Improvement: Positive LogMAR change = vision gain\n",
    "   - Capped at LogMAR 1.3 for very poor vision\n\n",
    
    "3. Clustering Results:\n",
    "   - Number of clusters: 2\n",
    "   - Fuzzy coefficient (m): ", round(comprehensive_m, 3), "\n",
    "   - Better outcome group (Cluster ", cluster_interpretation$better_cluster, "): ", 
    sum(comprehensive_clusters_result$max_cluster == cluster_interpretation$better_cluster), " patients\n",
    "   - Worse outcome group (Cluster ", cluster_interpretation$worse_cluster, "): ", 
    sum(comprehensive_clusters_result$max_cluster == cluster_interpretation$worse_cluster), " patients\n\n",
    
    "4. Statistical Significance:\n",
    "   - Significant vision parameters (LogMAR): ", vision_sig, "\n",
    "   - Significant BF macular parameters: ", bf_mac_sig, "\n",
    "   - Significant TH macular parameters: ", th_mac_sig, "\n",
    "   - Total significant parameters: ", total_sig, "\n\n",
    
    "5. Key Findings (LogMAR Interpretation):\n",
    "   - Vision improvements: Better cluster advantages in ", 
    cluster_interpretation$vision_cluster2_better + cluster_interpretation$vision_cluster1_better, " parameters\n",
    "   - Blood flow macular: Better cluster advantages in ", 
    ifelse(cluster_interpretation$better_cluster == 2, cluster_interpretation$bf_mac_cluster2_better, cluster_interpretation$bf_mac_cluster1_better), " parameters\n",
    "   - Thickness macular: Better cluster advantages in ", 
    ifelse(cluster_interpretation$better_cluster == 2, cluster_interpretation$th_mac_cluster2_better, cluster_interpretation$th_mac_cluster1_better), " parameters\n",
    "   - LogMAR analysis enables standardized vision comparison\n\n",
    
    "6. Clinical Implications:\n",
    "   - LogMAR provides standardized vision assessment\n",
    "   - Identifies patients with comprehensive good outcomes\n",
    "   - Focuses on macular region most relevant for visual function\n",
    "   - Provides targeted assessment of central retinal recovery\n",
    "   - Can guide macular-specific follow-up strategies\n\n",
    
    "7. LogMAR Analysis Advantages:\n",
    "   ✓ Standardized vision measurement scale\n",
    "   ✓ Linear relationship with visual function\n",
    "   ✓ Appropriate for statistical analysis\n",
    "   ✓ Clinically relevant improvement thresholds\n",
    "   ✓ Compatible with international standards\n\n",
    
    "8. Output Files (LogMAR Version):\n",
    "   - ppv_comprehensive_cluster_results_with_outcomes_logmar.csv: Main clustering results\n",
    "   - ppv_comprehensive_cluster_statistics_logmar.csv: Statistical analysis\n",
    "   - ppv_comprehensive_parameter_importance_logmar.csv: Parameter importance ranking\n",
    "   - ppv_comprehensive_cluster_summary_logmar.csv: Cluster summary statistics\n",
    "   - plots/*_logmar.*: LogMAR-specific visualizations\n\n",
    
    "9. Visualization Updates:\n",
    "   - All vision plots now use LogMAR scale\n",
    "   - Proper axis labels for LogMAR interpretation\n",
    "   - Clinical thresholds indicated where appropriate\n",
    "   - Improvement direction clearly marked\n\n",
    
    "10. Statistical Interpretation Notes:\n",
    "    - Pre-vision (LogMAR): Negative cluster difference = Cluster 2 better baseline\n",
    "    - Vision improvement (LogMAR): Positive cluster difference = Cluster 2 better improvement\n",
    "    - Clinical significance: 0.1 LogMAR ≈ 1 line on eye chart\n",
    "    - Meaningful improvement: ≥ 0.3 LogMAR (3 lines)\n\n"
  )
  
  # Save report
  writeLines(report, "PPV_Comprehensive_Clustering_Report_LogMAR.txt")
  cat(report)
}

# Generate comprehensive report (LogMAR version)
generate_comprehensive_report()

# -------------------- 16. Quality control and validation (LogMAR version) --------------------
# Perform quality control checks (with LogMAR considerations)
quality_control <- function() {
  cat("===== Quality Control Summary (LogMAR Vision) =====\n")
  
  # Check cluster stability
  membership_stability <- comprehensive_clusters_result %>%
    summarise(
      Mean_Membership = mean(max_membership),
      Min_Membership = min(max_membership),
      Low_Membership_Count = sum(max_membership < 0.6)
    )
  
  cat("Cluster stability:\n")
  cat("- Mean membership confidence: ", round(membership_stability$Mean_Membership, 3), "\n")
  cat("- Minimum membership confidence: ", round(membership_stability$Min_Membership, 3), "\n")
  cat("- Patients with low confidence (<0.6): ", membership_stability$Low_Membership_Count, "\n\n")
  
  # Check parameter distribution
  param_summary <- comprehensive_stats %>%
    group_by(Data_Type) %>%
    summarise(
      Total_Params = n(),
      Significant_Params = sum(Significant == "Yes"),
      Sig_Percentage = round(sum(Significant == "Yes")/n()*100, 1)
    )
  
  cat("Parameter significance by type (LogMAR version):\n")
  print(param_summary)
  
  # LogMAR-specific quality checks
  if("vision_improvement_1m" %in% names(comprehensive_data_with_clusters)) {
    vision_range <- range(comprehensive_data_with_clusters$vision_improvement_1m, na.rm = TRUE)
    cat("\nLogMAR Vision Improvement Range: ", round(vision_range[1], 3), " to ", round(vision_range[2], 3), "\n")
    
    meaningful_improvement <- sum(comprehensive_data_with_clusters$vision_improvement_1m > 0.1, na.rm = TRUE)
    cat("Patients with meaningful improvement (>0.1 LogMAR): ", meaningful_improvement, "\n")
  }
  
  if("pre_vision" %in% names(comprehensive_data_with_clusters)) {
    baseline_range <- range(comprehensive_data_with_clusters$pre_vision, na.rm = TRUE)
    cat("LogMAR Baseline Vision Range: ", round(baseline_range[1], 3), " to ", round(baseline_range[2], 3), "\n")
    
    good_baseline <- sum(comprehensive_data_with_clusters$pre_vision < 0.3, na.rm = TRUE)
    cat("Patients with good baseline vision (<0.3 LogMAR): ", good_baseline, "\n")
  }
  
  # Save quality control results
  write.csv(membership_stability, "quality_control_membership_logmar.csv", row.names = FALSE)
  write.csv(param_summary, "quality_control_parameters_logmar.csv", row.names = FALSE)
}

# Perform quality control (LogMAR version)
quality_control()

# -------------------- 17. LogMAR-specific summary statistics --------------------
# Create LogMAR-specific interpretation summary
create_logmar_summary <- function() {
  cat("\n===== LogMAR Conversion Summary =====\n")
  
  if("vision_improvement_1m" %in% names(comprehensive_data_with_clusters)) {
    # Vision improvement analysis
    improvement_summary <- comprehensive_data_with_clusters %>%
      group_by(max_cluster) %>%
      summarise(
        n = n(),
        Mean_Improvement_LogMAR = mean(vision_improvement_1m, na.rm = TRUE),
        SD_Improvement_LogMAR = sd(vision_improvement_1m, na.rm = TRUE),
        Patients_Improved = sum(vision_improvement_1m > 0, na.rm = TRUE),
        Patients_Meaningful_Improvement = sum(vision_improvement_1m > 0.1, na.rm = TRUE),
        Patients_Substantial_Improvement = sum(vision_improvement_1m > 0.3, na.rm = TRUE),
        .groups = 'drop'
      )
    
    cat("Vision Improvement by Cluster (LogMAR):\n")
    print(improvement_summary)
    
    # Statistical test for vision improvement
    if(length(unique(comprehensive_data_with_clusters$max_cluster)) >= 2) {
      vision_test <- t.test(vision_improvement_1m ~ max_cluster, data = comprehensive_data_with_clusters)
      cat("\nVision Improvement Comparison (t-test):\n")
      cat("p-value:", round(vision_test$p.value, 4), "\n")
      cat("Mean difference:", round(diff(vision_test$estimate), 3), "LogMAR\n")
      cat("Clinical interpretation: ", 
          ifelse(abs(diff(vision_test$estimate)) > 0.1, "Clinically meaningful difference", "Small difference"), "\n")
    }
  }
  
  if("pre_vision" %in% names(comprehensive_data_with_clusters)) {
    # Baseline vision analysis
    baseline_summary <- comprehensive_data_with_clusters %>%
      group_by(max_cluster) %>%
      summarise(
        Mean_Baseline_LogMAR = mean(pre_vision, na.rm = TRUE),
        SD_Baseline_LogMAR = sd(pre_vision, na.rm = TRUE),
        Patients_Good_Vision = sum(pre_vision < 0.3, na.rm = TRUE),
        Patients_Moderate_Vision = sum(pre_vision >= 0.3 & pre_vision < 0.7, na.rm = TRUE),
        Patients_Poor_Vision = sum(pre_vision >= 0.7, na.rm = TRUE),
        .groups = 'drop'
      )
    
    cat("\nBaseline Vision by Cluster (LogMAR):\n")
    print(baseline_summary)
  }
  
  # Save LogMAR summaries
  if(exists("improvement_summary")) {
    write.csv(improvement_summary, "logmar_improvement_summary_by_cluster.csv", row.names = FALSE)
  }
  if(exists("baseline_summary")) {
    write.csv(baseline_summary, "logmar_baseline_summary_by_cluster.csv", row.names = FALSE)
  }
  
  cat("\n✓ LogMAR analysis completed with proper statistical interpretation\n")
}

# Create LogMAR summary
create_logmar_summary()

# -------------------- Final Summary --------------------
cat("\n========================================\n")
cat("✅ OCTA + Vision (LogMAR) Clustering Analysis Complete!\n")  
cat("========================================\n\n")

cat("🔄 Key Modifications Made:\n")
cat("✓ Vision parameters converted to LogMAR before clustering\n")
cat("✓ Statistical tests interpreted correctly for LogMAR scale\n")
cat("✓ Visualizations updated with proper LogMAR labels\n")
cat("✓ Clinical interpretation adjusted for LogMAR thresholds\n")
cat("✓ Quality control includes LogMAR-specific checks\n\n")

cat("📊 LogMAR Interpretation Guide:\n")
cat("- Pre-vision LogMAR: Lower values = better vision\n")
cat("- Vision improvement LogMAR: Positive values = improvement\n")
cat("- Clinical significance: 0.1 LogMAR ≈ 1 line on eye chart\n")
cat("- Meaningful improvement: ≥ 0.3 LogMAR (3 lines)\n\n")

cat("📁 Key Output Files:\n")
cat("- PPV_Comprehensive_Clustering_Report_LogMAR.txt: Complete analysis report\n")
cat("- ppv_comprehensive_cluster_results_with_outcomes_logmar.csv: Main results\n")
cat("- plots/*_logmar.*: LogMAR-corrected visualizations\n")
cat("- logmar_*_summary_by_cluster.csv: LogMAR-specific summaries\n\n")

cat("🎯 Next Steps:\n")
cat("1. Review LogMAR-corrected visualizations\n")
cat("2. Validate clinical interpretation with domain experts\n")
cat("3. Consider vision improvement thresholds for clinical decisions\n")
cat("4. Use standardized LogMAR scale for future comparisons\n\n")

cat("========================================\n")