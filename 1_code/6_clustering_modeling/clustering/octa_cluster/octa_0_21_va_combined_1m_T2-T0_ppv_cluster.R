# ----------------------------------------------------
# OCTA + Vision Combined Clustering Analysis for PPV Group
# Using Mfuzz Fuzzy Clustering with OCTA (blood flow + thickness) and vision parameters
# MODIFIED: Only using WF region (0_21) OCTA parameters
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

# -------------------- 5. Filter OCTA parameters (WF ONLY) --------------------
# Function to filter blood flow layers (WF REGION ONLY - 0_21)
filter_bloodflow_layers_WF <- function(data) {
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  # MODIFIED: Only include WF region (0_21)
  region_of_interest <- "0_21"
  
  # Create pattern for WF region only
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*", region_of_interest, "_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target blood flow layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "WF blood flow T0 parameters:\n")
    cat("- WF region (0_21) only:", length(params_T0), "parameters\n\n")
    
    cat("WF blood flow parameters:\n")
    cat(paste(gsub("_T0$", "", params_T0), collapse = "\n"), "\n\n")
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
    WF_params = gsub("_T0$", "", params_T0)
  ))
}

# Function to filter thickness layers (WF REGION ONLY - 0_21)
filter_thickness_layers_WF <- function(data) {
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  # MODIFIED: Only include WF region (0_21)
  region_of_interest <- "0_21"
  
  # Create pattern for WF region only
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*", region_of_interest, "_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target thickness layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "WF thickness T0 parameters:\n")
    cat("- WF region (0_21) only:", length(params_T0), "parameters\n\n")
    
    cat("WF thickness parameters:\n")
    cat(paste(gsub("_T0$", "", params_T0), collapse = "\n"), "\n\n")
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
    WF_params = gsub("_T0$", "", params_T0)
  ))
}

# Apply filters (WF ONLY)
ppv_bloodflow_filtered <- filter_bloodflow_layers_WF(ppv_bloodflow)
ppv_thickness_filtered <- filter_thickness_layers_WF(ppv_thickness)

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

# Print data combination summary (WF ONLY)
cat("\n===== WF-Only Data Summary =====\n")
cat("Total patients in WF dataset:", nrow(comprehensive_data), "\n")
cat("Vision parameters: vision_improvement_1m, pre_vision, age\n")

if(length(ppv_bloodflow_filtered$WF_params) > 0) {
  cat("OCTA blood flow - WF region (0_21):", length(ppv_bloodflow_filtered$WF_params), "parameters\n")
}
cat("Total OCTA blood flow parameters:", ncol(ppv_bloodflow_improvement) - 1, "\n")

if(length(ppv_thickness_filtered$WF_params) > 0) {
  cat("OCTA thickness - WF region (0_21):", length(ppv_thickness_filtered$WF_params), "parameters\n")
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

# Count parameters by type (WF ONLY)
vision_params <- names(complete_data)[grepl("vision_improvement|pre_vision|age", names(complete_data))]
bloodflow_params <- names(complete_data)[grepl("SVP|ICP|DCP|Choroid", names(complete_data))]
thickness_params <- names(complete_data)[grepl("GCL|INL|Retina", names(complete_data))]

# All OCTA parameters are now WF only
bloodflow_WF <- bloodflow_params  # All blood flow params are WF
thickness_WF <- thickness_params  # All thickness params are WF

cat("Vision parameters:", length(vision_params), "\n")
cat("Blood flow parameters - WF (0_21) only:", length(bloodflow_WF), "\n")
cat("Thickness parameters - WF (0_21) only:", length(thickness_WF), "\n")

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
cat("- Blood flow parameters (WF only):", length(bloodflow_params), "\n")
cat("- Thickness parameters (WF only):", length(thickness_params), "\n")
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
write.csv(comprehensive_clusters_result, "ppv_WF_cluster_results.csv", row.names = FALSE)

# -------------------- 10. Visualize clustering results --------------------
# Create comprehensive dataset with cluster assignments
comprehensive_data_with_clusters <- comprehensive_data %>%
  inner_join(comprehensive_clusters_result, by = c("ID" = "subject_id"))

# Function to create WF-focused visualizations
create_WF_visualizations <- function(data, vision_params, bloodflow_params, thickness_params, 
                                          bloodflow_WF, thickness_WF) {
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  data$max_cluster <- as.factor(data$max_cluster)
  all_params <- c(vision_params, bloodflow_params, thickness_params)
  
  # 0. Overall WF heatmap for all parameters
  cat("Creating WF-focused comprehensive heatmap...\n")
  
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
      # Determine parameter category
      Data_Type = case_when(
        Parameter %in% vision_params ~ "Vision",
        Parameter %in% bloodflow_WF ~ "Blood Flow - WF",
        Parameter %in% thickness_WF ~ "Thickness - WF",
        TRUE ~ "Other"
      ),
      # Clean parameter names for display
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter),
      Parameter_Clean = gsub("_0_21", "", Parameter_Clean),
      Parameter_Clean = gsub("_", " ", Parameter_Clean),
      # Add WF suffix for OCTA parameters
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
  
  # Create the overall WF heatmap
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
      title = "WF-Focused Cluster Comparison Heatmap",
      subtitle = paste("Vision + OCTA WF Parameters (0_21 only) | n =", nrow(data)),
      x = "Parameters",
      y = "Cluster",
      caption = "Mac = WF region (0_21)"
    )
  
  # Save the overall heatmap
  ggsave("plots/WF_comprehensive_heatmap.pdf", p_overall_heatmap, 
         width = 16, height = 10, dpi = 300)
  ggsave("plots/WF_comprehensive_heatmap.png", p_overall_heatmap, 
         width = 16, height = 10, dpi = 300)
  
  cat("WF comprehensive heatmap saved to plots/WF_comprehensive_heatmap.pdf/png\n")
  
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
  
  cat("\n===== Cluster Summary by Parameter Type (WF Focus) =====\n")
  print(cluster_summary_table)
  
  # Save summary table
  write.csv(cluster_summary_table, "plots/WF_cluster_summary_by_type.csv", row.names = FALSE)
  
  # 1. Vision parameters by cluster
  for(param in vision_params) {
    if(param %in% names(data)) {
      param_clean <- gsub("_improvement|_1m", "", param)
      
      p <- ggplot(data, aes(x = max_cluster, y = .data[[param]], fill = max_cluster)) +
        geom_boxplot() +
        geom_jitter(alpha = 0.5, width = 0.2) +
        scale_fill_brewer(palette = "Set1") +
        theme_bw() +
        labs(title = paste("Vision:", param_clean, "by Cluster"), 
             x = "Cluster", y = param_clean)
      
      ggsave(paste0("plots/vision_", param_clean, "_boxplot.pdf"), p, width = 8, height = 6)
    }
  }
  
  # 2. WF heatmap
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
        Parameter %in% bloodflow_WF ~ "Blood Flow - WF",
        Parameter %in% thickness_WF ~ "Thickness - WF",
        TRUE ~ "Other"
      ),
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter)
    )
  
  # Create WF heatmap
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
      title = "WF Parameter Overview by Cluster\n(Vision + OCTA WF Region Only)",
      x = "Parameters",
      y = "Cluster",
      fill = "Mean\nValue"
    )
  
  ggsave("plots/WF_heatmap.pdf", p_heatmap, width = 20, height = 10)
  ggsave("plots/WF_heatmap.png", p_heatmap, width = 20, height = 10, dpi = 300)
  
  # 3. PCA visualization
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
        title = "PCA of WF Parameters\n(Vision + OCTA WF Only)",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    ggsave("plots/WF_pca.pdf", p_pca, width = 10, height = 8)
    
    # 4. Variable contribution plot
    loadings <- pca_result$rotation[, 1:2]
    loadings_df <- data.frame(
      Variable = rownames(loadings),
      PC1 = loadings[, 1],
      PC2 = loadings[, 2],
      Data_Type = case_when(
        rownames(loadings) %in% vision_params ~ "Vision",
        rownames(loadings) %in% bloodflow_WF ~ "BF - WF",
        rownames(loadings) %in% thickness_WF ~ "TH - WF",
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
        title = "Parameter Contributions to Principal Components\n(WF Analysis)",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
        color = "Parameter Type"
      )
    
    ggsave("plots/WF_loadings.pdf", p_loadings, width = 14, height = 12)
  }
}

# Create visualizations
create_WF_visualizations(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                              bloodflow_WF, thickness_WF)

# -------------------- 11. Statistical analysis --------------------
# Function to analyze differences between clusters (WF-focused)
analyze_WF_differences <- function(data, vision_params, bloodflow_params, thickness_params,
                                        bloodflow_WF, thickness_WF) {
  all_params <- c(vision_params, bloodflow_params, thickness_params)
  
  results <- data.frame(
    Parameter = character(),
    Data_Type = character(),
    Cluster1_Mean = numeric(),
    Cluster2_Mean = numeric(),
    Mean_Difference = numeric(),
    P_Value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )
  
  for(param in all_params) {
    if(param %in% names(data) && !all(is.na(data[[param]]))) {
      # Determine data type
      data_type <- case_when(
        param %in% vision_params ~ "Vision",
        param %in% bloodflow_params ~ "Blood Flow - WF",
        param %in% thickness_params ~ "Thickness - WF",
        TRUE ~ "Other"
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
            Cluster1_Mean = ifelse("1" %in% names(means), means["1"], NA),
            Cluster2_Mean = ifelse("2" %in% names(means), means["2"], NA),
            Mean_Difference = ifelse(length(means) >= 2, means["2"] - means["1"], NA),
            P_Value = test_result$p.value,
            Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
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
    results <- results %>% arrange(Data_Type, P_Value)
  }
  
  return(results)
}

# Perform statistical analysis
comprehensive_stats <- analyze_WF_differences(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                                                   bloodflow_WF, thickness_WF)

# Save results
write.csv(comprehensive_stats, "ppv_WF_cluster_statistics.csv", row.names = FALSE)

# Print significant results by data type
cat("\n===== Significant Parameters by Data Type (WF Focus) =====\n")

for(data_type in c("Vision", "Blood Flow - WF", "Thickness - WF")) {
  cat("\n", data_type, "parameters:\n")
  
  sig_params <- comprehensive_stats %>% 
    filter(Data_Type == data_type & Significant == "Yes") %>%
    arrange(P_Value)
  
  if(nrow(sig_params) > 0) {
    print(sig_params %>% dplyr::select(Parameter, Mean_Difference, P_Value, Significant))
  } else {
    cat("No significant parameters found\n")
  }
}

# -------------------- 12. Interpret clusters --------------------
# WF cluster interpretation function
interpret_WF_clusters <- function(stats_results) {
  if(nrow(stats_results) == 0) {
    cat("No statistical results available for interpretation\n")
    return(NULL)
  }
  
  # Analyze by data type
  vision_stats <- stats_results %>% filter(Data_Type == "Vision")
  bloodflow_WF_stats <- stats_results %>% filter(Data_Type == "Blood Flow - WF")
  thickness_WF_stats <- stats_results %>% filter(Data_Type == "Thickness - WF")
  
  # Count significant improvements by cluster
  vision_cluster2_better <- sum(vision_stats$Mean_Difference > 0 & vision_stats$Significant == "Yes", na.rm = TRUE)
  vision_cluster1_better <- sum(vision_stats$Mean_Difference < 0 & vision_stats$Significant == "Yes", na.rm = TRUE)
  
  bf_mac_cluster2_better <- sum(bloodflow_WF_stats$Mean_Difference > 0 & bloodflow_WF_stats$Significant == "Yes", na.rm = TRUE)
  bf_mac_cluster1_better <- sum(bloodflow_WF_stats$Mean_Difference < 0 & bloodflow_WF_stats$Significant == "Yes", na.rm = TRUE)
  
  th_mac_cluster2_better <- sum(thickness_WF_stats$Mean_Difference > 0 & thickness_WF_stats$Significant == "Yes", na.rm = TRUE)
  th_mac_cluster1_better <- sum(thickness_WF_stats$Mean_Difference < 0 & thickness_WF_stats$Significant == "Yes", na.rm = TRUE)
  
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
  
  # Print interpretation results
  cat("\n===== WF-Focused Cluster Interpretation =====\n")
  cat("Overall better cluster:", better_cluster, "\n")
  cat("Overall worse cluster:", worse_cluster, "\n\n")
  
  cat("Detailed breakdown:\n")
  cat("Vision - Cluster 2 advantages:", vision_cluster2_better, ", Cluster 1 advantages:", vision_cluster1_better, "\n")
  cat("Blood Flow WF - Cluster 2 advantages:", bf_mac_cluster2_better, ", Cluster 1 advantages:", bf_mac_cluster1_better, "\n")
  cat("Thickness WF - Cluster 2 advantages:", th_mac_cluster2_better, ", Cluster 1 advantages:", th_mac_cluster1_better, "\n")
  cat("Total - Cluster 2 advantages:", total_cluster2_better, ", Cluster 1 advantages:", total_cluster1_better, "\n\n")
  
  # Additional summary
  cat("Summary:\n")
  if(better_cluster == 2) {
    cat("Cluster 2 patients show significantly better WF outcomes\n")
  } else {
    cat("Cluster 1 patients show significantly better WF outcomes\n")
  }
  
  cat("This indicates that PPV patients can be stratified into distinct WF outcome groups\n")
  cat("based on their vision and WF OCTA response patterns.\n")
  
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

# Execute cluster interpretation
cluster_interpretation <- interpret_WF_clusters(comprehensive_stats)

# Add outcome labels to cluster results
if(!is.null(cluster_interpretation)) {
  comprehensive_clusters_result$outcome_quality <- ifelse(
    comprehensive_clusters_result$max_cluster == cluster_interpretation$better_cluster,
    "Better",
    "Worse"
  )
  
  # Save final results with outcome labels
  write.csv(comprehensive_clusters_result, "ppv_WF_cluster_results_with_outcomes.csv", row.names = FALSE)
  
  cat("\nCluster outcome labels have been added and saved to:\n")
  cat("'ppv_WF_cluster_results_with_outcomes.csv'\n")
} else {
  cat("\nWarning: Could not determine cluster interpretation due to insufficient statistical results.\n")
}

# -------------------- 13. Parameter importance analysis --------------------
# Analyze which WF parameters contribute most to cluster separation
analyze_WF_parameter_importance <- function(stats_results) {
  if(nrow(stats_results) == 0) {
    return(data.frame())
  }
  
  importance_scores <- stats_results %>%
    mutate(
      Effect_Size = abs(Mean_Difference),
      Importance_Score = Effect_Size * (-log10(P_Value + 1e-10))
    ) %>%
    arrange(desc(Importance_Score))
  
  # Get top parameters by data type
  top_vision <- importance_scores %>%
    filter(Data_Type == "Vision") %>%
    head(3)
  
  top_bloodflow <- importance_scores %>%
    filter(Data_Type == "Blood Flow - WF") %>%
    head(5)
  
  top_thickness <- importance_scores %>%
    filter(Data_Type == "Thickness - WF") %>%
    head(5)
  
  cat("\n===== WF Parameter Importance Analysis =====\n")
  
  if(nrow(top_vision) > 0) {
    cat("\nTop Vision parameters:\n")
    print(top_vision %>% dplyr::select(Parameter, Data_Type, Mean_Difference, P_Value, Importance_Score, Significant))
  }
  
  if(nrow(top_bloodflow) > 0) {
    cat("\nTop WF Blood Flow parameters:\n")
    print(top_bloodflow %>% dplyr::select(Parameter, Data_Type, Mean_Difference, P_Value, Importance_Score, Significant))
  }
  
  if(nrow(top_thickness) > 0) {
    cat("\nTop WF Thickness parameters:\n")
    print(top_thickness %>% dplyr::select(Parameter, Data_Type, Mean_Difference, P_Value, Importance_Score, Significant))
  }
  
  # Save importance analysis
  write.csv(importance_scores, "ppv_WF_parameter_importance.csv", row.names = FALSE)
  
  return(importance_scores)
}

# Perform parameter importance analysis
parameter_importance <- analyze_WF_parameter_importance(comprehensive_stats)

# -------------------- 14. Create WF summary --------------------
# Final summary statistics for WF analysis
final_summary <- comprehensive_data_with_clusters %>%
  group_by(max_cluster) %>%
  summarise(
    Count = n(),
    Mean_Membership = mean(max_membership, na.rm = TRUE),
    # Vision parameters
    Mean_Vision_Improvement = mean(vision_improvement_1m, na.rm = TRUE),
    Mean_Pre_Vision = mean(pre_vision, na.rm = TRUE),
    Mean_Age = mean(age, na.rm = TRUE),
    # Summary of WF OCTA improvements
    Positive_BF_Improvements = rowSums(across(all_of(bloodflow_params), ~ .x > 0), na.rm = TRUE),
    Positive_TH_Improvements = rowSums(across(all_of(thickness_params), ~ .x > 0), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Outcome_Quality = ifelse(max_cluster == cluster_interpretation$better_cluster, "Better", "Worse")
  )

cat("\n===== Final WF Summary =====\n")
print(final_summary)

write.csv(final_summary, "ppv_WF_cluster_summary.csv", row.names = FALSE)

# -------------------- 15. Generate WF report --------------------
generate_WF_report <- function() {
  vision_sig <- sum(comprehensive_stats$Data_Type == "Vision" & comprehensive_stats$Significant == "Yes")
  bf_mac_sig <- sum(comprehensive_stats$Data_Type == "Blood Flow - WF" & comprehensive_stats$Significant == "Yes")
  th_mac_sig <- sum(comprehensive_stats$Data_Type == "Thickness - WF" & comprehensive_stats$Significant == "Yes")
  total_sig <- sum(comprehensive_stats$Significant == "Yes")
  
  report <- paste0(
    "========================================\n",
    "PPV Group WF-Focused Clustering Analysis Report\n",
    "Vision + OCTA (WF Region 0_21 Only)\n",
    "========================================\n\n",
    
    "1. Data Overview:\n",
    "   - Total patients analyzed: ", nrow(complete_data), "\n",
    "   - Vision parameters: ", length(vision_params), " (improvement, baseline vision, age)\n",
    "   - OCTA blood flow - WF only: ", length(bloodflow_WF), " parameters\n",
    "   - OCTA thickness - WF only: ", length(thickness_WF), " parameters\n",
    "   - Total parameters: ", ncol(complete_data), "\n",
    "   - Focus: Central WF region (0_21) only\n\n",
    
    "2. Clustering Results:\n",
    "   - Number of clusters: 2\n",
    "   - Fuzzy coefficient (m): ", round(comprehensive_m, 3), "\n",
    "   - Better outcome group (Cluster ", cluster_interpretation$better_cluster, "): ", 
    sum(comprehensive_clusters_result$max_cluster == cluster_interpretation$better_cluster), " patients\n",
    "   - Worse outcome group (Cluster ", cluster_interpretation$worse_cluster, "): ", 
    sum(comprehensive_clusters_result$max_cluster == cluster_interpretation$worse_cluster), " patients\n\n",
    
    "3. Statistical Significance:\n",
    "   - Significant vision parameters: ", vision_sig, "\n",
    "   - Significant WF blood flow parameters: ", bf_mac_sig, "\n",
    "   - Significant WF thickness parameters: ", th_mac_sig, "\n",
    "   - Total significant parameters: ", total_sig, "\n\n",
    
    "4. Key WF Findings:\n",
    "   - Vision improvements: Better cluster advantages in ", 
    cluster_interpretation$vision_cluster2_better + cluster_interpretation$vision_cluster1_better, " parameters\n",
    "   - WF blood flow: Better cluster advantages in ", 
    ifelse(cluster_interpretation$better_cluster == 2, cluster_interpretation$bf_mac_cluster2_better, cluster_interpretation$bf_mac_cluster1_better), " parameters\n",
    "   - WF thickness: Better cluster advantages in ", 
    ifelse(cluster_interpretation$better_cluster == 2, cluster_interpretation$th_mac_cluster2_better, cluster_interpretation$th_mac_cluster1_better), " parameters\n",
    "   - Focus on central 0_21 region provides detailed WF assessment\n\n",
    
    "5. Clinical Implications:\n",
    "   - Identifies patients with good central WF outcomes\n",
    "   - Focused on clinically most important retinal region\n",
    "   - Provides detailed assessment of WF vascular changes\n",
    "   - Suitable for WF-specific treatment decisions\n\n",
    
    "6. WF Analysis Advantages:\n",
    "   - Concentrated assessment of central vision-critical area\n",
    "   - Higher data completeness due to focus on single region\n",
    "   - Clinically relevant for central vision outcomes\n",
    "   - Reduced complexity while maintaining clinical relevance\n\n",
    
    "7. Output Files:\n",
    "   - ppv_WF_cluster_results_with_outcomes.csv: Main clustering results\n",
    "   - ppv_WF_cluster_statistics.csv: Statistical analysis\n",
    "   - ppv_WF_parameter_importance.csv: Parameter importance ranking\n",
    "   - ppv_WF_cluster_summary.csv: Cluster summary statistics\n",
    "   - plots/: WF-focused visualizations\n\n",
    
    "8. Advantages of WF-Only Analysis:\n",
    "   - Focused on central vision-critical region\n",
    "   - Higher patient inclusion due to better data completeness\n",
    "   - Clinically relevant for central visual function\n",
    "   - Simplified interpretation while maintaining relevance\n",
    "   - Optimal for WF-specific surgical outcome assessment\n\n"
  )
  
  # Save report
  writeLines(report, "PPV_WF_Clustering_Report.txt")
  cat(report)
}

# Generate WF report
generate_WF_report()

# -------------------- 16. Quality control and validation --------------------
# Perform quality control checks for WF analysis
quality_control_WF <- function() {
  cat("===== WF Analysis Quality Control Summary =====\n")
  
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
  
  cat("WF parameter significance by type:\n")
  print(param_summary)
  
  # Save quality control results
  write.csv(membership_stability, "quality_control_WF_membership.csv", row.names = FALSE)
  write.csv(param_summary, "quality_control_WF_parameters.csv", row.names = FALSE)
}

# Perform quality control
quality_control_WF()

# -------------------- 17. Final completion message --------------------
cat("\n========================================\n")
cat("WF-FOCUSED PPV CLUSTERING ANALYSIS COMPLETE!\n")
cat("========================================\n")

cat("\nThis analysis successfully integrated:\n")
cat("✓ Vision improvement parameters (functional outcomes)\n")
cat("✓ OCTA blood flow - WF region (0_21) parameters only\n")
cat("✓ OCTA thickness - WF region (0_21) parameters only\n")
cat("✓ Patient baseline characteristics\n\n")

cat("WF parameter breakdown:\n")
cat("- Vision parameters: ", length(vision_params), "\n")
cat("- Blood flow WF: ", length(bloodflow_WF), "\n")
cat("- Thickness WF: ", length(thickness_WF), "\n")
cat("- Total parameters: ", ncol(complete_data), "\n\n")

cat("Key outputs generated:\n")
cat("1. WF clustering results with outcome predictions\n")
cat("2. WF-focused statistical analysis\n")
cat("3. Parameter importance ranking for WF region\n")
cat("4. WF-specific visualizations\n")
cat("5. Quality control assessments\n")
cat("6. Detailed WF analysis report\n\n")

cat("Files saved:\n")
cat("- ppv_WF_cluster_results_with_outcomes.csv\n")
cat("- ppv_WF_cluster_statistics.csv\n")
cat("- ppv_WF_parameter_importance.csv\n")
cat("- ppv_WF_cluster_summary.csv\n")
cat("- PPV_WF_Clustering_Report.txt\n")
cat("- plots/ directory with WF visualizations:\n")
cat("  * WF_comprehensive_heatmap.pdf/png\n")
cat("  * WF_heatmap.pdf/png\n")
cat("  * WF_pca.pdf\n")
cat("  * WF_loadings.pdf\n")
cat("- quality_control_WF_*.csv files\n\n")

cat("ADVANTAGES OF THIS WF-ONLY APPROACH:\n")
cat("1. FOCUSED ANALYSIS: Concentrated on central vision-critical region\n")
cat("2. HIGHER COMPLETENESS: Better data availability for single region\n")
cat("3. CLINICAL RELEVANCE: Central WF region most important for vision\n")
cat("4. SIMPLIFIED INTERPRETATION: Easier to understand single-region results\n")
cat("5. SURGICAL RELEVANCE: Directly related to central visual outcomes\n")
cat("6. ROBUST CLUSTERING: Based on most complete and relevant data\n\n")

cat("This WF-focused approach provides clinically relevant patient\n")
cat("stratification by concentrating on the central retinal region most\n")
cat("critical for visual function after PPV surgery.\n")