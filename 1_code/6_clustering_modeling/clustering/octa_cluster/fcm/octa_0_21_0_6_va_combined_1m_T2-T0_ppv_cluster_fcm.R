# ----------------------------------------------------
# OCTA + Vision Combined Clustering Analysis for PPV Group
# Using FCM (Fuzzy C-Means) Clustering with OCTA (blood flow + thickness) and vision parameters
# Comprehensive analysis integrating anatomical and functional improvements
# ----------------------------------------------------

# Load required libraries
library(tidyverse)
library(e1071)         # For FCM clustering
library(cluster)       # For additional clustering metrics
library(ggplot2)
library(factoextra)
library(corrplot)
library(r4projects)
library(pheatmap)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# -------------------- 1. Load data --------------------
# Load baseline information and OCTA data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/fcm_comprehensive", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/fcm_comprehensive")

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

# -------------------- 5. Filter OCTA parameters --------------------
# Function to filter blood flow layers (enhanced for both macular and wide-field)
filter_bloodflow_layers <- function(data) {
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  # Include both 0_21 (macular) and 0_6 (wide-field) regions
  regions_of_interest <- c("0_21", "0_6")
  
  # Create pattern for both regions
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target blood flow layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "target blood flow layer T0 parameters:\n")
    
    # Separate macular and wide-field parameters for better reporting
    macular_params <- params_T0[grep("0_21_T0$", params_T0)]
    widefield_params <- params_T0[grep("0_6_T0$", params_T0)]
    
    cat("- Macular region (0_21):", length(macular_params), "parameters\n")
    cat("- Wide-field region (0_6):", length(widefield_params), "parameters\n")
    cat("- Total:", length(params_T0), "parameters\n\n")
    
    if(length(macular_params) > 0) {
      cat("Macular parameters:\n")
      cat(paste(gsub("_T0$", "", macular_params), collapse = "\n"), "\n\n")
    }
    
    if(length(widefield_params) > 0) {
      cat("Wide-field parameters:\n")
      cat(paste(gsub("_T0$", "", widefield_params), collapse = "\n"), "\n\n")
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
    widefield_params = gsub("_T0$", "", widefield_params)
  ))
}

# Function to filter thickness layers (enhanced for both macular and wide-field)
filter_thickness_layers <- function(data) {
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  # Include both 0_21 (macular) and 0_6 (wide-field) regions
  regions_of_interest <- c("0_21", "0_6")
  
  # Create pattern for both regions
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target thickness layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "target thickness layer T0 parameters:\n")
    
    # Separate macular and wide-field parameters for better reporting
    macular_params <- params_T0[grep("0_21_T0$", params_T0)]
    widefield_params <- params_T0[grep("0_6_T0$", params_T0)]
    
    cat("- Macular region (0_21):", length(macular_params), "parameters\n")
    cat("- Wide-field region (0_6):", length(widefield_params), "parameters\n")
    cat("- Total:", length(params_T0), "parameters\n\n")
    
    if(length(macular_params) > 0) {
      cat("Macular thickness parameters:\n")
      cat(paste(gsub("_T0$", "", macular_params), collapse = "\n"), "\n\n")
    }
    
    if(length(widefield_params) > 0) {
      cat("Wide-field thickness parameters:\n")
      cat(paste(gsub("_T0$", "", widefield_params), collapse = "\n"), "\n\n")
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
    widefield_params = gsub("_T0$", "", widefield_params)
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

# Merge with vision data
comprehensive_data <- ppv_vision %>%
  dplyr::select(ID, vision_improvement_1m, pre_vision, age) %>%
  full_join(ppv_octa_combined, by = "ID")

# Print data combination summary (enhanced for regions)
cat("\n===== Comprehensive Data Summary =====\n")
cat("Total patients in comprehensive dataset:", nrow(comprehensive_data), "\n")
cat("Vision parameters: vision_improvement_1m, pre_vision, age\n")

if(length(ppv_bloodflow_filtered$macular_params) > 0) {
  cat("OCTA blood flow - Macular region (0_21):", length(ppv_bloodflow_filtered$macular_params), "parameters\n")
}
if(length(ppv_bloodflow_filtered$widefield_params) > 0) {
  cat("OCTA blood flow - Wide-field region (0_6):", length(ppv_bloodflow_filtered$widefield_params), "parameters\n")
}
cat("Total OCTA blood flow parameters:", ncol(ppv_bloodflow_improvement) - 1, "\n")

if(length(ppv_thickness_filtered$macular_params) > 0) {
  cat("OCTA thickness - Macular region (0_21):", length(ppv_thickness_filtered$macular_params), "parameters\n")
}
if(length(ppv_thickness_filtered$widefield_params) > 0) {
  cat("OCTA thickness - Wide-field region (0_6):", length(ppv_thickness_filtered$widefield_params), "parameters\n")
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

# Count parameters by type and region
vision_params <- names(complete_data)[grepl("vision_improvement|pre_vision|age", names(complete_data))]
bloodflow_params <- names(complete_data)[grepl("SVP|ICP|DCP|Choroid", names(complete_data))]
thickness_params <- names(complete_data)[grepl("GCL|INL|Retina", names(complete_data))]

# Further categorize by region
bloodflow_macular <- bloodflow_params[grepl("0_21", bloodflow_params)]
bloodflow_widefield <- bloodflow_params[grepl("0_6", bloodflow_params)]
thickness_macular <- thickness_params[grepl("0_21", thickness_params)]
thickness_widefield <- thickness_params[grepl("0_6", thickness_params)]

cat("Vision parameters:", length(vision_params), "\n")
cat("Blood flow parameters - Total:", length(bloodflow_params), "\n")
cat("  - Macular (0_21):", length(bloodflow_macular), "\n")
cat("  - Wide-field (0_6):", length(bloodflow_widefield), "\n")
cat("Thickness parameters - Total:", length(thickness_params), "\n")
cat("  - Macular (0_21):", length(thickness_macular), "\n")
cat("  - Wide-field (0_6):", length(thickness_widefield), "\n")

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
cat("- Vision parameters:", length(vision_params), "\n")
cat("- Blood flow parameters:", length(bloodflow_params), "\n")
cat("  * Macular region:", length(bloodflow_macular), "\n")
cat("  * Wide-field region:", length(bloodflow_widefield), "\n")
cat("- Thickness parameters:", length(thickness_params), "\n")
cat("  * Macular region:", length(thickness_macular), "\n")
cat("  * Wide-field region:", length(thickness_widefield), "\n")
cat("Patients used:", nrow(comprehensive_data_std), "\n")
cat("Patient IDs:", paste(complete_ids, collapse=", "), "\n")

# -------------------- 9. Perform FCM clustering --------------------
# FCM clustering function with optimal parameter selection
perform_fcm_clustering <- function(data, k_range = 2:4, m_range = seq(1.5, 3.0, 0.5)) {
  
  cat("\n===== FCM Clustering Parameter Optimization =====\n")
  
  # Convert to matrix for FCM
  data_matrix <- as.matrix(data)
  
  # Function to calculate clustering validity indices
  calculate_validity_indices <- function(data_matrix, fcm_result, k) {
    # Partition coefficient (PC)
    membership_matrix <- fcm_result$membership
    pc <- sum(membership_matrix^2) / nrow(membership_matrix)
    
    # Partition entropy (PE)
    pe <- -sum(membership_matrix * log(membership_matrix, base = 2)) / nrow(membership_matrix)
    
    # Dunn index using hard clustering
    hard_clusters <- apply(membership_matrix, 1, which.max)
    
    # Calculate silhouette if possible
    silhouette_score <- NA
    if(k > 1 && k < nrow(data_matrix) && length(unique(hard_clusters)) > 1) {
      tryCatch({
        sil <- silhouette(hard_clusters, dist(data_matrix))
        silhouette_score <- mean(sil[, 3])
      }, error = function(e) {
        silhouette_score <- NA
      })
    }
    
    return(list(
      pc = pc,
      pe = pe,
      silhouette = silhouette_score,
      hard_clusters = hard_clusters
    ))
  }
  
  # Test different parameters
  results <- data.frame()
  best_result <- NULL
  best_score <- -Inf
  
  for(k in k_range) {
    for(m in m_range) {
      cat(sprintf("Testing k=%d, m=%.1f\n", k, m))
      
      tryCatch({
        # Perform FCM clustering
        set.seed(123) # For reproducibility
        fcm_result <- cmeans(data_matrix, centers = k, m = m, iter.max = 100)
        
        # Calculate validity indices
        validity <- calculate_validity_indices(data_matrix, fcm_result, k)
        
        # Combined score (higher PC, lower PE, higher silhouette)
        combined_score <- validity$pc - validity$pe + ifelse(is.na(validity$silhouette), 0, validity$silhouette)
        
        # Store results
        results <- rbind(results, data.frame(
          k = k,
          m = m,
          pc = validity$pc,
          pe = validity$pe,
          silhouette = validity$silhouette,
          combined_score = combined_score,
          convergence = fcm_result$iter < 100
        ))
        
        # Update best result
        if(combined_score > best_score && fcm_result$iter < 100) {
          best_score <- combined_score
          best_result <- list(
            fcm = fcm_result,
            k = k,
            m = m,
            validity = validity,
            score = combined_score
          )
        }
        
      }, error = function(e) {
        cat(sprintf("Error with k=%d, m=%.1f: %s\n", k, m, e$message))
      })
    }
  }
  
  # Print optimization results
  cat("\n===== Parameter Optimization Results =====\n")
  results_sorted <- results %>% arrange(desc(combined_score))
  print(head(results_sorted, 10))
  
  if(!is.null(best_result)) {
    cat(sprintf("\nOptimal parameters: k=%d, m=%.1f\n", best_result$k, best_result$m))
    cat(sprintf("Best combined score: %.4f\n", best_result$score))
    cat(sprintf("PC: %.4f, PE: %.4f, Silhouette: %.4f\n", 
                best_result$validity$pc, best_result$validity$pe, 
                ifelse(is.na(best_result$validity$silhouette), 0, best_result$validity$silhouette)))
  }
  
  return(list(
    best_result = best_result,
    all_results = results_sorted,
    optimization_successful = !is.null(best_result)
  ))
}

# Perform FCM clustering with optimization
fcm_optimization <- perform_fcm_clustering(comprehensive_data_std)

if(!fcm_optimization$optimization_successful) {
  # Fallback to default parameters
  cat("Optimization failed, using default parameters: k=2, m=2.0\n")
  set.seed(123)
  comprehensive_fcm <- cmeans(as.matrix(comprehensive_data_std), centers = 2, m = 2.0)
  optimal_k <- 2
  optimal_m <- 2.0
} else {
  comprehensive_fcm <- fcm_optimization$best_result$fcm
  optimal_k <- fcm_optimization$best_result$k
  optimal_m <- fcm_optimization$best_result$m
}

# Get clustering results
comprehensive_membership <- comprehensive_fcm$membership
comprehensive_main_clusters <- apply(comprehensive_membership, 1, which.max)
comprehensive_max_membership <- apply(comprehensive_membership, 1, max)

# Create results dataframe
comprehensive_clusters_result <- data.frame(
  subject_id = complete_ids,
  max_cluster = comprehensive_main_clusters,
  max_membership = comprehensive_max_membership
)

# Add individual cluster memberships
for(i in 1:optimal_k) {
  comprehensive_clusters_result[[paste0("membership_c", i)]] <- comprehensive_membership[, i]
}

# Save clustering results
write.csv(comprehensive_clusters_result, "ppv_fcm_comprehensive_cluster_results.csv", row.names = FALSE)
write.csv(fcm_optimization$all_results, "ppv_fcm_optimization_results.csv", row.names = FALSE)

cat("\n===== FCM Clustering Results =====\n")
cat("Optimal number of clusters:", optimal_k, "\n")
cat("Optimal fuzziness parameter m:", optimal_m, "\n")
cat("Converged in", comprehensive_fcm$iter, "iterations\n")
cat("Final objective function value:", comprehensive_fcm$withinerror, "\n")

# Cluster size distribution
cluster_sizes <- table(comprehensive_main_clusters)
cat("Cluster sizes:\n")
print(cluster_sizes)

cat("Mean membership confidence by cluster:\n")
for(i in 1:optimal_k) {
  cluster_members <- comprehensive_main_clusters == i
  mean_conf <- mean(comprehensive_max_membership[cluster_members])
  cat(sprintf("Cluster %d: %.3f\n", i, mean_conf))
}

# -------------------- 10. Visualize clustering results --------------------
# Create comprehensive dataset with cluster assignments
comprehensive_data_with_clusters <- comprehensive_data %>%
  inner_join(comprehensive_clusters_result, by = c("ID" = "subject_id"))

# Function to create enhanced visualizations with regional analysis
create_comprehensive_visualizations <- function(data, vision_params, bloodflow_params, thickness_params, 
                                                bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield) {
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  data$max_cluster <- as.factor(data$max_cluster)
  all_params <- c(vision_params, bloodflow_params, thickness_params)
  
  # 0. Overall comprehensive heatmap for all parameters
  cat("Creating overall comprehensive heatmap...\n")
  
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
        Parameter %in% vision_params ~ "Vision",
        Parameter %in% bloodflow_macular ~ "Blood Flow - Macular",
        Parameter %in% bloodflow_widefield ~ "Blood Flow - Wide-field", 
        Parameter %in% thickness_macular ~ "Thickness - Macular",
        Parameter %in% thickness_widefield ~ "Thickness - Wide-field",
        TRUE ~ "Other"
      ),
      # Clean parameter names for display
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter),
      Parameter_Clean = gsub("_0_21|_0_6", "", Parameter_Clean),
      Parameter_Clean = gsub("_", " ", Parameter_Clean),
      # Add region suffix for OCTA parameters
      Parameter_Display = case_when(
        grepl("0_21", Parameter) ~ paste0(Parameter_Clean, " (Mac)"),
        grepl("0_6", Parameter) ~ paste0(Parameter_Clean, " (WF)"),
        TRUE ~ Parameter_Clean
      )
    ) %>%
    # Order parameters by type and name
    arrange(Data_Type, Parameter_Display)
  
  # Create factor levels for proper ordering
  overall_plot_data$Parameter_Display <- factor(overall_plot_data$Parameter_Display, 
                                                levels = unique(overall_plot_data$Parameter_Display))
  
  # Create the overall comprehensive heatmap
  p_overall_heatmap <- ggplot(overall_plot_data, aes(x = Parameter_Display, y = as.factor(max_cluster), fill = Mean_Value)) +
    geom_tile(color = "white", size = 0.5) +
    # Add text labels with values
    geom_text(aes(label = sprintf("%.2f", Mean_Value)), 
              color = "black", size = 2.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0,
      name = "Mean\nValue",
      guide = guide_colorbar(barwidth = 1, barheight = 10)
    ) +
    facet_wrap(~ Data_Type, scales = "free_x", ncol = 1, strip.position = "top") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 12, face = "bold", color = "darkblue"),
      strip.background = element_rect(fill = "lightgray"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    ) +
    labs(
      title = "FCM Comprehensive Cluster Comparison Heatmap",
      subtitle = paste("All Parameters: Vision + OCTA (Macular & Wide-field) | k =", optimal_k, ", m =", optimal_m, "| n =", nrow(data)),
      x = "Parameters",
      y = "Cluster",
      caption = "Mac = Macular (0_21), WF = Wide-field (0_6)"
    )
  
  # Save the overall heatmap
  ggsave("plots/fcm_overall_comprehensive_heatmap.pdf", p_overall_heatmap, 
         width = 20, height = 12, dpi = 300)
  ggsave("plots/fcm_overall_comprehensive_heatmap.png", p_overall_heatmap, 
         width = 20, height = 12, dpi = 300)
  
  cat("Overall comprehensive heatmap saved to plots/fcm_overall_comprehensive_heatmap.pdf/png\n")
  
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
  
  cat("\n===== FCM Cluster Summary by Parameter Type =====\n")
  print(cluster_summary_table)
  
  # Save summary table
  write.csv(cluster_summary_table, "plots/fcm_cluster_summary_by_type.csv", row.names = FALSE)
  
  # 1. Vision parameters by cluster
  for(param in vision_params) {
    if(param %in% names(data)) {
      param_clean <- gsub("_improvement|_1m", "", param)
      
      p <- ggplot(data, aes(x = max_cluster, y = .data[[param]], fill = max_cluster)) +
        geom_boxplot() +
        geom_jitter(alpha = 0.5, width = 0.2) +
        scale_fill_brewer(palette = "Set1") +
        theme_bw() +
        labs(title = paste("FCM Vision:", param_clean, "by Cluster"), 
             x = "Cluster", y = param_clean)
      
      ggsave(paste0("plots/fcm_vision_", param_clean, "_boxplot.pdf"), p, width = 8, height = 6)
    }
  }
  
  # 2. Enhanced heatmap with regional breakdown
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
        Parameter %in% vision_params ~ "Vision",
        Parameter %in% bloodflow_macular ~ "Blood Flow - Macular",
        Parameter %in% bloodflow_widefield ~ "Blood Flow - Wide-field",
        Parameter %in% thickness_macular ~ "Thickness - Macular", 
        Parameter %in% thickness_widefield ~ "Thickness - Wide-field",
        TRUE ~ "Other"
      ),
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter)
    )
  
  # Create comprehensive heatmap with regional breakdown
  p_heatmap <- ggplot(plot_data, aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    facet_wrap(~ Data_Type, scales = "free_x", ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    labs(
      title = "FCM Comprehensive Parameter Overview by Cluster\n(Including Macular and Wide-field Regions)",
      x = "Parameters",
      y = "Cluster",
      fill = "Mean\nValue"
    )
  
  ggsave("plots/fcm_comprehensive_regional_heatmap.pdf", p_heatmap, width = 24, height = 12)
  ggsave("plots/fcm_comprehensive_regional_heatmap.png", p_heatmap, width = 24, height = 12, dpi = 300)
  
  # 3. Regional comparison heatmaps
  # Macular vs Wide-field Blood Flow
  if(length(bloodflow_macular) > 0 && length(bloodflow_widefield) > 0) {
    bf_regional_data <- plot_data %>%
      filter(grepl("Blood Flow", Data_Type))
    
    p_bf_regional <- ggplot(bf_regional_data, aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      facet_wrap(~ Data_Type, scales = "free_x", ncol = 1) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "FCM Blood Flow Parameters: Macular vs Wide-field Comparison",
        x = "Blood Flow Parameters",
        y = "Cluster",
        fill = "Mean\nImprovement"
      )
    
    ggsave("plots/fcm_bloodflow_regional_comparison.pdf", p_bf_regional, width = 16, height = 8)
  }
  
  # Macular vs Wide-field Thickness
  if(length(thickness_macular) > 0 && length(thickness_widefield) > 0) {
    th_regional_data <- plot_data %>%
      filter(grepl("Thickness", Data_Type))
    
    p_th_regional <- ggplot(th_regional_data, aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      facet_wrap(~ Data_Type, scales = "free_x", ncol = 1) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "FCM Thickness Parameters: Macular vs Wide-field Comparison",
        x = "Thickness Parameters",
        y = "Cluster",
        fill = "Mean\nImprovement"
      )
    
    ggsave("plots/fcm_thickness_regional_comparison.pdf", p_th_regional, width = 16, height = 8)
  }
  
  # 4. PCA visualization
  if(length(all_params) > 2) {
    pca_data <- data %>%
      dplyr::select(all_of(all_params))
    
    pca_result <- prcomp(pca_data, scale. = TRUE)
    
    # Create PCA plot data frame
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
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      labs(
        title = "FCM PCA of Comprehensive Parameters\n(Vision + OCTA Macular + Wide-field)",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    ggsave("plots/fcm_comprehensive_regional_pca.pdf", p_pca, width = 10, height = 8)
    
    # 5. Enhanced variable contribution plot with regional colors
    loadings <- pca_result$rotation[, 1:2]
    loadings_df <- data.frame(
      Variable = rownames(loadings),
      PC1 = loadings[, 1],
      PC2 = loadings[, 2],
      Data_Type = case_when(
        rownames(loadings) %in% vision_params ~ "Vision",
        rownames(loadings) %in% bloodflow_macular ~ "BF - Macular",
        rownames(loadings) %in% bloodflow_widefield ~ "BF - Wide-field",
        rownames(loadings) %in% thickness_macular ~ "TH - Macular", 
        rownames(loadings) %in% thickness_widefield ~ "TH - Wide-field",
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
        title = "FCM Parameter Contributions to Principal Components\n(Regional Analysis)",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
        color = "Parameter Type"
      )
    
    ggsave("plots/fcm_comprehensive_regional_loadings.pdf", p_loadings, width = 14, height = 12)
  }
  
  # 6. FCM-specific membership plots
  # Membership matrix heatmap
  membership_data <- data %>%
    dplyr::select(starts_with("membership_c"))
  
  if(ncol(membership_data) > 1) {
    membership_long <- membership_data %>%
      mutate(Patient = 1:nrow(membership_data)) %>%
      pivot_longer(cols = starts_with("membership_c"), names_to = "Cluster", values_to = "Membership") %>%
      mutate(Cluster = gsub("membership_c", "Cluster ", Cluster))
    
    p_membership <- ggplot(membership_long, aes(x = Patient, y = Cluster, fill = Membership)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkblue", name = "Membership") +
      theme_minimal() +
      labs(
        title = "FCM Membership Matrix Heatmap",
        subtitle = paste("Fuzziness parameter m =", optimal_m),
        x = "Patient",
        y = "Cluster"
      )
    
    ggsave("plots/fcm_membership_heatmap.pdf", p_membership, width = 12, height = 6)
  }
  
  # Membership confidence distribution
  p_conf_dist <- ggplot(data, aes(x = max_membership)) +
    geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mean(data$max_membership), color = "red", linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(
      title = "FCM Cluster Assignment Confidence Distribution",
      subtitle = paste("Mean confidence:", round(mean(data$max_membership), 3)),
      x = "Maximum Membership Value",
      y = "Number of Patients"
    )
  
  ggsave("plots/fcm_confidence_distribution.pdf", p_conf_dist, width = 10, height = 6)
}

# Create visualizations
create_comprehensive_visualizations(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                                    bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield)

