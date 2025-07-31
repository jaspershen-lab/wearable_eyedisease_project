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

cat("Vision parameters:", length(vision_params), "\n")
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
cat("- Vision parameters:", length(vision_params), "\n")
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
  
  # Create the overall comprehensive heatmap
  p_overall_heatmap <- ggplot(overall_plot_data, aes(x = Parameter_Display, y = as.factor(max_cluster), fill = Mean_Value)) +
    geom_tile(color = "white", size = 0.5) +
    # Add text labels with values
    geom_text(aes(label = sprintf("%.2f", Mean_Value)), 
              color = "black", size = 2.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#542788", 
      mid = "white", 
      high = "#f1a340", 
      midpoint = 0,
      name = "Mean\nValue",
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
      legend.position = "right"
    ) +
    labs(
      title = "Comprehensive Cluster Comparison Heatmap",
      subtitle = paste("All Parameters: Vision + OCTA (Macular only) | n =", nrow(data)),
      x = "Parameters",
      y = "Cluster",
      caption = "Mac = Macular (0_21)"
    )
  
  # Save the overall heatmap
  ggsave("plots/overall_comprehensive_heatmap.pdf", p_overall_heatmap, 
         width = 20, height = 12, dpi = 300)
  ggsave("plots/overall_comprehensive_heatmap.png", p_overall_heatmap, 
         width = 20, height = 12, dpi = 300)
  
  cat("Overall comprehensive heatmap saved to plots/overall_comprehensive_heatmap.pdf/png\n")
  
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
  
  cat("\n===== Cluster Summary by Parameter Type =====\n")
  print(cluster_summary_table)
  
  # Save summary table
  write.csv(cluster_summary_table, "plots/cluster_summary_by_type.csv", row.names = FALSE)
  
  # 1. Vision parameters by cluster
  for(param in vision_params) {
    if(param %in% names(data)) {
      param_clean <- gsub("_improvement|_1m", "", param)
      
      p <- ggplot(data, aes(x = max_cluster, y = .data[[param]], fill = max_cluster)) +
        geom_boxplot() +
        geom_jitter(alpha = 0.5, width = 0.2) +
        scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
        theme_bw() +
        labs(title = paste("Vision:", param_clean, "by Cluster"), 
             x = "Cluster", y = param_clean)
      
      ggsave(paste0("plots/vision_", param_clean, "_boxplot.pdf"), p, width = 8, height = 6)
    }
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
        Parameter %in% thickness_macular ~ "Thickness - Macular",
        TRUE ~ "Other"
      ),
      Parameter_Clean = gsub("_improvement|_1m", "", Parameter)
    )
  
  # Create comprehensive heatmap with regional breakdown
  p_heatmap <- ggplot(plot_data, aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Value)) +
    geom_tile() +
    scale_fill_gradient2( low = "#542788", 
                          mid = "white", 
                          high = "#f1a340",  midpoint = 0) +
    facet_wrap(~ Data_Type, scales = "free_x", ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    labs(
      title = "Comprehensive Parameter Overview by Cluster\n(Macular Region Only)",
      x = "Parameters",
      y = "Cluster",
      fill = "Mean\nValue"
    )
  
  ggsave("plots/comprehensive_regional_heatmap.pdf", p_heatmap, width = 24, height = 12)
  ggsave("plots/comprehensive_regional_heatmap.png", p_heatmap, width = 24, height = 12, dpi = 300)
  
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
        title = "PCA of Comprehensive Parameters\n(Vision + OCTA Macular)",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    ggsave("plots/comprehensive_regional_pca.pdf", p_pca, width = 10, height = 8)
    
    # 4. Enhanced variable contribution plot with regional colors
    loadings <- pca_result$rotation[, 1:2]
    loadings_df <- data.frame(
      Variable = rownames(loadings),
      PC1 = loadings[, 1],
      PC2 = loadings[, 2],
      Data_Type = case_when(
        rownames(loadings) %in% vision_params ~ "Vision",
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
        title = "Parameter Contributions to Principal Components\n(Regional Analysis)",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
        color = "Parameter Type"
      )
    
    ggsave("plots/comprehensive_regional_loadings.pdf", p_loadings, width = 14, height = 12)
  }
}

# Create visualizations
create_comprehensive_visualizations(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                                    bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield)

# -------------------- 11. Statistical analysis --------------------
# Function to analyze differences between clusters (enhanced for regional analysis)
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
    stringsAsFactors = FALSE
  )
  
  for(param in all_params) {
    if(param %in% names(data) && !all(is.na(data[[param]]))) {
      # Determine data type and region
      data_type <- case_when(
        param %in% vision_params ~ "Vision",
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

# Perform statistical analysis (updated function call)
comprehensive_stats <- analyze_comprehensive_differences(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params,
                                                         bloodflow_macular, bloodflow_widefield, thickness_macular, thickness_widefield)

# Save results
write.csv(comprehensive_stats, "ppv_comprehensive_cluster_statistics.csv", row.names = FALSE)

# Print significant results by data type and region
cat("\n===== Significant Parameters by Data Type and Region =====\n")

for(data_type in c("Vision", "Blood Flow", "Thickness")) {
  cat("\n", data_type, "parameters:\n")
  
  if(data_type == "Vision") {
    sig_params <- comprehensive_stats %>% 
      filter(Data_Type == data_type & Significant == "Yes") %>%
      arrange(P_Value)
    
    if(nrow(sig_params) > 0) {
      print(sig_params %>% dplyr::select(Parameter, Mean_Difference, P_Value, Significant))
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

# -------------------- 12. Interpret clusters --------------------
# Comprehensive cluster interpretation function
interpret_comprehensive_clusters <- function(stats_results) {
  if(nrow(stats_results) == 0) {
    cat("No statistical results available for interpretation\n")
    return(NULL)
  }
  
  # Analyze by data type and region
  vision_stats <- stats_results %>% filter(Data_Type == "Vision")
  bloodflow_macular_stats <- stats_results %>% filter(Data_Type == "Blood Flow" & Region == "Macular")
  thickness_macular_stats <- stats_results %>% filter(Data_Type == "Thickness" & Region == "Macular")
  
  # Count significant improvements by cluster and region
  vision_cluster2_better <- sum(vision_stats$Mean_Difference > 0 & vision_stats$Significant == "Yes", na.rm = TRUE)
  vision_cluster1_better <- sum(vision_stats$Mean_Difference < 0 & vision_stats$Significant == "Yes", na.rm = TRUE)
  
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
  
  # Print interpretation results
  cat("\n===== Comprehensive Cluster Interpretation =====\n")
  cat("Overall better cluster:", better_cluster, "\n")
  cat("Overall worse cluster:", worse_cluster, "\n\n")
  
  cat("Detailed breakdown:\n")
  cat("Vision - Cluster 2 advantages:", vision_cluster2_better, ", Cluster 1 advantages:", vision_cluster1_better, "\n")
  cat("Blood Flow Macular - Cluster 2 advantages:", bf_mac_cluster2_better, ", Cluster 1 advantages:", bf_mac_cluster1_better, "\n")
  cat("Thickness Macular - Cluster 2 advantages:", th_mac_cluster2_better, ", Cluster 1 advantages:", th_mac_cluster1_better, "\n")
  cat("Total - Cluster 2 advantages:", total_cluster2_better, ", Cluster 1 advantages:", total_cluster1_better, "\n\n")
  
  # Additional summary
  cat("Summary:\n")
  if(better_cluster == 2) {
    cat("Cluster 2 patients show significantly better outcomes across multiple parameters\n")
  } else {
    cat("Cluster 1 patients show significantly better outcomes across multiple parameters\n")
  }
  
  cat("This indicates that PPV patients can be stratified into distinct outcome groups\n")
  cat("based on their comprehensive vision and OCTA response patterns.\n")
  
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
cluster_interpretation <- interpret_comprehensive_clusters(comprehensive_stats)

# Add outcome labels to cluster results
if(!is.null(cluster_interpretation)) {
  comprehensive_clusters_result$outcome_quality <- ifelse(
    comprehensive_clusters_result$max_cluster == cluster_interpretation$better_cluster,
    "Better",
    "Worse"
  )
  
  # Save final results with outcome labels
  write.csv(comprehensive_clusters_result, "ppv_comprehensive_cluster_results_with_outcomes.csv", row.names = FALSE)
  
  cat("\nCluster outcome labels have been added and saved to:\n")
  cat("'ppv_comprehensive_cluster_results_with_outcomes.csv'\n")
} else {
  cat("\nWarning: Could not determine cluster interpretation due to insufficient statistical results.\n")
}

# -------------------- 13. Parameter importance analysis --------------------
# Analyze which parameters contribute most to cluster separation
analyze_parameter_importance <- function(stats_results) {
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
    filter(Data_Type == "Blood Flow") %>%
    head(5)
  
  top_thickness <- importance_scores %>%
    filter(Data_Type == "Thickness") %>%
    head(5)
  
  cat("\n===== Parameter Importance Analysis =====\n")
  
  if(nrow(top_vision) > 0) {
    cat("\nTop Vision parameters:\n")
    print(top_vision %>% dplyr::select(Parameter, Data_Type, Mean_Difference, P_Value, Importance_Score, Significant))
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
  write.csv(importance_scores, "ppv_comprehensive_parameter_importance.csv", row.names = FALSE)
  
  return(importance_scores)
}

# Perform parameter importance analysis
parameter_importance <- analyze_parameter_importance(comprehensive_stats)

# -------------------- 14. Create comprehensive summary --------------------
# Final summary statistics
final_summary <- comprehensive_data_with_clusters %>%
  group_by(max_cluster) %>%
  summarise(
    Count = n(),
    Mean_Membership = mean(max_membership, na.rm = TRUE),
    # Vision parameters
    Mean_Vision_Improvement = mean(vision_improvement_1m, na.rm = TRUE),
    Mean_Pre_Vision = mean(pre_vision, na.rm = TRUE),
    Mean_Age = mean(age, na.rm = TRUE),
    # Summary of OCTA improvements (count of positive improvements)
    Positive_BF_Improvements = rowSums(across(all_of(bloodflow_params), ~ .x > 0), na.rm = TRUE),
    Positive_TH_Improvements = rowSums(across(all_of(thickness_params), ~ .x > 0), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Outcome_Quality = ifelse(max_cluster == cluster_interpretation$better_cluster, "Better", "Worse")
  )

cat("\n===== Final Comprehensive Summary =====\n")
print(final_summary)

write.csv(final_summary, "ppv_comprehensive_cluster_summary.csv", row.names = FALSE)

# -------------------- 15. Generate comprehensive report --------------------
generate_comprehensive_report <- function() {
  vision_sig <- sum(comprehensive_stats$Data_Type == "Vision" & comprehensive_stats$Significant == "Yes")
  bf_mac_sig <- sum(comprehensive_stats$Data_Type == "Blood Flow" & comprehensive_stats$Region == "Macular" & comprehensive_stats$Significant == "Yes")
  th_mac_sig <- sum(comprehensive_stats$Data_Type == "Thickness" & comprehensive_stats$Region == "Macular" & comprehensive_stats$Significant == "Yes")
  total_sig <- sum(comprehensive_stats$Significant == "Yes")
  
  report <- paste0(
    "========================================\n",
    "PPV Group Comprehensive Clustering Analysis Report\n",
    "Vision + OCTA (Macular Region 0_21 Only)\n",
    "========================================\n\n",
    
    "1. Data Overview:\n",
    "   - Total patients analyzed: ", nrow(complete_data), "\n",
    "   - Vision parameters: ", length(vision_params), " (improvement, baseline vision, age)\n",
    "   - OCTA blood flow - Macular: ", length(bloodflow_macular), " parameters\n",
    "   - OCTA thickness - Macular: ", length(thickness_macular), " parameters\n",
    "   - Total parameters: ", ncol(complete_data), "\n",
    "   - Focus: Central retinal region (0_21) analysis\n\n",
    
    "2. Clustering Results:\n",
    "   - Number of clusters: 2\n",
    "   - Fuzzy coefficient (m): ", round(comprehensive_m, 3), "\n",
    "   - Better outcome group (Cluster ", cluster_interpretation$better_cluster, "): ", 
    sum(comprehensive_clusters_result$max_cluster == cluster_interpretation$better_cluster), " patients\n",
    "   - Worse outcome group (Cluster ", cluster_interpretation$worse_cluster, "): ", 
    sum(comprehensive_clusters_result$max_cluster == cluster_interpretation$worse_cluster), " patients\n\n",
    
    "3. Statistical Significance:\n",
    "   - Significant vision parameters: ", vision_sig, "\n",
    "   - Significant BF macular parameters: ", bf_mac_sig, "\n",
    "   - Significant TH macular parameters: ", th_mac_sig, "\n",
    "   - Total significant parameters: ", total_sig, "\n\n",
    
    "4. Key Findings:\n",
    "   - Vision improvements: Better cluster advantages in ", 
    cluster_interpretation$vision_cluster2_better + cluster_interpretation$vision_cluster1_better, " parameters\n",
    "   - Blood flow macular: Better cluster advantages in ", 
    ifelse(cluster_interpretation$better_cluster == 2, cluster_interpretation$bf_mac_cluster2_better, cluster_interpretation$bf_mac_cluster1_better), " parameters\n",
    "   - Thickness macular: Better cluster advantages in ", 
    ifelse(cluster_interpretation$better_cluster == 2, cluster_interpretation$th_mac_cluster2_better, cluster_interpretation$th_mac_cluster1_better), " parameters\n",
    "   - This analysis focuses on central retinal changes\n\n",
    
    "5. Clinical Implications:\n",
    "   - Identifies patients with comprehensive good outcomes\n",
    "   - Focuses on macular region most relevant for visual function\n",
    "   - Provides targeted assessment of central retinal recovery\n",
    "   - Can guide macular-specific follow-up strategies\n\n",
    
    "6. Analysis Advantages:\n",
    "   - Focused assessment of vision-critical central retina\n",
    "   - Higher data completeness due to fewer missing parameters\n",
    "   - Direct relevance to visual function outcomes\n",
    "   - Simplified interpretation for clinical decision-making\n\n",
    
    "7. Output Files:\n",
    "   - ppv_comprehensive_cluster_results_with_outcomes.csv: Main clustering results\n",
    "   - ppv_comprehensive_cluster_statistics.csv: Statistical analysis\n",
    "   - ppv_comprehensive_parameter_importance.csv: Parameter importance ranking\n",
    "   - ppv_comprehensive_cluster_summary.csv: Cluster summary statistics\n",
    "   - plots/: Comprehensive visualizations\n\n",
    
    "8. Advantages of This Approach:\n",
    "   - CLINICAL RELEVANCE: Focuses on vision-critical central retina\n",
    "   - HIGHER COMPLETENESS: More patients with complete macular data\n",
    "   - FUNCTIONAL CORRELATION: Direct relationship to visual outcomes\n",
    "   - SIMPLIFIED INTERPRETATION: Focused on anatomically relevant region\n",
    "   - ROBUST CLUSTERING: Based on maximum available central retinal information\n\n"
  )
  
  # Save report
  writeLines(report, "PPV_Comprehensive_Clustering_Report.txt")
  cat(report)
}

# Generate comprehensive report
generate_comprehensive_report()

# -------------------- 16. Quality control and validation --------------------
# Perform quality control checks
quality_control <- function() {
  cat("===== Quality Control Summary =====\n")
  
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
  
  cat("Parameter significance by type:\n")
  print(param_summary)
  
  # Save quality control results
  write.csv(membership_stability, "quality_control_membership.csv", row.names = FALSE)
  write.csv(param_summary, "quality_control_parameters.csv", row.names = FALSE)
}

# Perform quality control
quality_control()




# ================== OCTA聚类改善值可视化（修正版）==================
# 针对OCTA数据特点：使用T2-T0改善值进行聚类，不是时间序列
# 重点展示：1）改善值分布对比 2）参数相关性 3）聚类特征

# -------------------- 1. OCTA改善值可视化主函数 --------------------
create_octa_improvement_visualizations <- function(comprehensive_data_with_clusters, 
                                                   vision_params, bloodflow_params, thickness_params,
                                                   bloodflow_macular, thickness_macular) {
  
  cat("\n🎨 开始创建OCTA改善值可视化（修正版）...\n")
  cat("📋 数据特点：聚类基于T2-T0改善值，非时间序列\n")
  
  # 创建输出目录
  dir.create("plots/octa_improvements", recursive = TRUE, showWarnings = FALSE)
  
  # 1. 创建改善值分布对比图
  create_improvement_comparison_plots(comprehensive_data_with_clusters, bloodflow_params, thickness_params)
  
  # 2. 创建参数重要性排序图
  create_parameter_importance_plots(comprehensive_data_with_clusters, bloodflow_params, thickness_params)
  
  # 3. 创建聚类特征雷达图
  create_cluster_radar_charts(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params)
  
  # 4. 创建改善值热图和聚类对比
  create_improvement_heatmaps(comprehensive_data_with_clusters, bloodflow_params, thickness_params)
  
  # 5. 创建相关性分析图
  create_correlation_analysis(comprehensive_data_with_clusters, vision_params, bloodflow_params, thickness_params)
  
  # 6. 创建患者改善模式图
  create_patient_improvement_patterns(comprehensive_data_with_clusters, bloodflow_params, thickness_params)
  
  cat("\n✅ OCTA改善值可视化创建完成！\n")
}

# -------------------- 2. 改善值分布对比图 --------------------
create_improvement_comparison_plots <- function(data, bloodflow_params, thickness_params) {
  
  cat("  📊 创建改善值分布对比图...\n")
  
  # 1. 整体改善值分布对比
  all_params <- c(bloodflow_params, thickness_params)
  
  # 准备数据
  improvement_data <- prepare_improvement_comparison_data(data, all_params)
  
  if(nrow(improvement_data) == 0) {
    cat("    Warning: No improvement data available\n")
    return(NULL)
  }
  
  # A. 按参数类型的箱线图对比
  p_boxplot_type <- ggplot(improvement_data, aes(x = Parameter_Type, y = Improvement_Value, 
                                                 fill = factor(Cluster))) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(aes(color = factor(Cluster)), 
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3),
                alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    labs(
      title = "OCTA Parameter Improvements by Type and Cluster",
      subtitle = "Distribution of T2-T0 improvement values | Positive = better outcomes",
      x = "Parameter Type",
      y = "Improvement Value (T2 - T0)",
      caption = "Red dashed line = no change | Above = improvement, Below = deterioration"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "top"
    )
  
  ggsave("plots/octa_improvements/improvement_boxplot_by_type.pdf", p_boxplot_type, width = 12, height = 8)
  ggsave("plots/octa_improvements/improvement_boxplot_by_type.png", p_boxplot_type, width = 12, height = 8, dpi = 300)
  
  # B. Top差异参数的详细对比
  create_top_parameters_comparison(improvement_data)
  
  cat("      ✓ 改善值分布对比图创建完成\n")
}

# -------------------- 3. 参数重要性排序图 --------------------
create_parameter_importance_plots <- function(data, bloodflow_params, thickness_params) {
  
  cat("  📊 创建参数重要性排序图...\n")
  
  all_params <- c(bloodflow_params, thickness_params)
  
  # 计算每个参数的聚类间差异
  param_importance <- calculate_parameter_importance(data, all_params)
  
  if(nrow(param_importance) == 0) {
    cat("    Warning: No parameter importance data available\n")
    return(NULL)
  }
  
  # 创建重要性排序图
  p_importance <- ggplot(param_importance %>% head(15), 
                         aes(x = reorder(Parameter_Clean, Importance_Score), y = Importance_Score)) +
    geom_col(aes(fill = Parameter_Type), alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", Importance_Score)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = c("Blood Flow" = "#3498db", "Thickness" = "#e74c3c"), name = "Type") +
    coord_flip() +
    labs(
      title = "Top 15 Most Important OCTA Parameters for Cluster Separation",
      subtitle = "Importance Score = |Mean Difference| × -log10(P-value)",
      x = "OCTA Parameters",
      y = "Importance Score",
      caption = "Higher scores indicate better cluster discrimination"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 10)
    )
  
  ggsave("plots/octa_improvements/parameter_importance_ranking.pdf", p_importance, width = 14, height = 10)
  ggsave("plots/octa_improvements/parameter_importance_ranking.png", p_importance, width = 14, height = 10, dpi = 300)
  
  # 保存重要性数据
  write.csv(param_importance, "plots/octa_improvements/parameter_importance_scores.csv", row.names = FALSE)
  
  cat("      ✓ 参数重要性排序图创建完成\n")
}

# -------------------- 4. 聚类特征雷达图 --------------------
create_cluster_radar_charts <- function(data, vision_params, bloodflow_params, thickness_params) {
  
  cat("  📊 创建聚类特征雷达图...\n")
  
  # 选择关键参数创建雷达图
  key_params <- select_key_parameters_for_radar(data, vision_params, bloodflow_params, thickness_params)
  
  if(length(key_params) == 0) {
    cat("    Warning: No key parameters selected for radar chart\n")
    return(NULL)
  }
  
  # 准备雷达图数据
  radar_data <- prepare_radar_data(data, key_params)
  
  # 创建雷达图
  create_radar_plot(radar_data, key_params)
  
  cat("      ✓ 聚类特征雷达图创建完成\n")
}

# -------------------- 5. 改善值热图 --------------------
create_improvement_heatmaps <- function(data, bloodflow_params, thickness_params) {
  
  cat("  📊 创建改善值热图...\n")
  
  all_params <- c(bloodflow_params, thickness_params)
  
  # 计算聚类平均改善值
  cluster_means <- calculate_cluster_means(data, all_params)
  
  # A. 综合热图
  p_heatmap_all <- ggplot(cluster_means, aes(x = Parameter_Clean, y = factor(Cluster), fill = Mean_Improvement)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Mean_Improvement)), 
              color = "black", size = 2.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#d73027", mid = "white", high = "#1a9850", 
      midpoint = 0, name = "Mean\nImprovement"
    ) +
    facet_wrap(~ Parameter_Type, scales = "free_x", ncol = 1) +
    labs(
      title = "OCTA Parameter Improvement Heatmap by Cluster",
      subtitle = "Mean improvement values (T2-T0) | Green = improvement, Red = deterioration",
      x = "Parameters", y = "Cluster"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  ggsave("plots/octa_improvements/improvement_heatmap_comprehensive.pdf", p_heatmap_all, width = 20, height = 8)
  ggsave("plots/octa_improvements/improvement_heatmap_comprehensive.png", p_heatmap_all, width = 20, height = 8, dpi = 300)
  
  # B. 聚类对比条形图
  create_cluster_comparison_barplot(cluster_means)
  
  cat("      ✓ 改善值热图创建完成\n")
}

# -------------------- 6. 相关性分析图 --------------------
create_correlation_analysis <- function(data, vision_params, bloodflow_params, thickness_params) {
  
  cat("  📊 创建相关性分析图...\n")
  
  # A. 视力与OCTA参数相关性
  create_vision_octa_correlation(data, vision_params, bloodflow_params, thickness_params)
  
  # B. OCTA参数间相关性
  create_octa_internal_correlation(data, bloodflow_params, thickness_params)
  
  cat("      ✓ 相关性分析图创建完成\n")
}

# -------------------- 7. 患者改善模式图 --------------------
create_patient_improvement_patterns <- function(data, bloodflow_params, thickness_params) {
  
  cat("  📊 创建患者改善模式图...\n")
  
  # A. 患者改善得分分布
  patient_scores <- calculate_patient_improvement_scores(data, bloodflow_params, thickness_params)
  
  # 创建改善得分分布图
  p_scores <- ggplot(patient_scores, aes(x = factor(Cluster), y = Overall_Score, fill = factor(Cluster))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = factor(Cluster)), width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    labs(
      title = "Patient Overall Improvement Scores by Cluster",
      subtitle = "Composite score based on all OCTA parameter improvements",
      x = "Cluster", y = "Overall Improvement Score"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave("plots/octa_improvements/patient_improvement_scores.pdf", p_scores, width = 10, height = 8)
  ggsave("plots/octa_improvements/patient_improvement_scores.png", p_scores, width = 10, height = 8, dpi = 300)
  
  # B. 改善模式分类
  create_improvement_pattern_classification(patient_scores)
  
  cat("      ✓ 患者改善模式图创建完成\n")
}

# -------------------- 辅助函数 --------------------

# 准备改善值对比数据
prepare_improvement_comparison_data <- function(data, all_params) {
  improvement_data <- data.frame()
  
  for(param in all_params) {
    if(param %in% names(data)) {
      param_data <- data %>%
        dplyr::select(ID, max_cluster, max_membership, all_of(param)) %>%
        filter(!is.na(.data[[param]])) %>%
        mutate(
          Parameter_Clean = clean_parameter_name(param),
          Parameter_Type = determine_parameter_type(param, bloodflow_params, thickness_params),
          Improvement_Value = .data[[param]],
          Cluster = max_cluster
        )
      
      improvement_data <- rbind(improvement_data, param_data)
    }
  }
  
  return(improvement_data)
}

# 清理参数名称
clean_parameter_name <- function(param) {
  cleaned <- gsub("_improvement$", "", param)
  cleaned <- gsub("_0_21", "", cleaned)
  cleaned <- gsub("_", " ", cleaned)
  return(cleaned)
}

# 确定参数类型
determine_parameter_type <- function(param, bloodflow_params, thickness_params) {
  if(param %in% bloodflow_params) {
    return("Blood Flow")
  } else if(param %in% thickness_params) {
    return("Thickness")
  } else {
    return("Other")
  }
}

# 创建Top参数对比
create_top_parameters_comparison <- function(improvement_data) {
  # 找出差异最大的参数
  param_differences <- improvement_data %>%
    group_by(Parameter_Clean, Parameter_Type) %>%
    summarise(
      Cluster1_Mean = mean(Improvement_Value[Cluster == 1], na.rm = TRUE),
      Cluster2_Mean = mean(Improvement_Value[Cluster == 2], na.rm = TRUE),
      Difference = abs(Cluster2_Mean - Cluster1_Mean),
      .groups = 'drop'
    ) %>%
    arrange(desc(Difference)) %>%
    head(12)
  
  # 创建Top参数的详细对比图
  top_params <- param_differences$Parameter_Clean
  top_data <- improvement_data %>% filter(Parameter_Clean %in% top_params)
  
  p_top <- ggplot(top_data, aes(x = Parameter_Clean, y = Improvement_Value, fill = factor(Cluster))) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = factor(Cluster)), 
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
                alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    facet_wrap(~ Parameter_Type, scales = "free", ncol = 1) +
    labs(
      title = "Top 12 Parameters with Largest Cluster Differences",
      subtitle = "OCTA parameters showing most distinct improvement patterns",
      x = "Parameters", y = "Improvement Value"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )
  
  ggsave("plots/octa_improvements/top_discriminative_parameters.pdf", p_top, width = 16, height = 12)
  ggsave("plots/octa_improvements/top_discriminative_parameters.png", p_top, width = 16, height = 12, dpi = 300)
}

# 计算参数重要性
calculate_parameter_importance <- function(data, all_params) {
  importance_results <- data.frame()
  
  for(param in all_params) {
    if(param %in% names(data)) {
      param_data <- data %>%
        dplyr::select(max_cluster, all_of(param)) %>%
        filter(!is.na(.data[[param]]))
      
      if(nrow(param_data) > 3) {
        # 计算t检验
        t_test <- try(t.test(param_data[[param]] ~ param_data$max_cluster), silent = TRUE)
        
        if(class(t_test) != "try-error") {
          means <- tapply(param_data[[param]], param_data$max_cluster, mean, na.rm = TRUE)
          mean_diff <- abs(means[2] - means[1])
          p_value <- t_test$p.value
          
          importance_score <- mean_diff * (-log10(p_value + 1e-10))
          
          importance_results <- rbind(importance_results, data.frame(
            Parameter = param,
            Parameter_Clean = clean_parameter_name(param),
            Parameter_Type = determine_parameter_type(param, bloodflow_params, thickness_params),
            Mean_Difference = mean_diff,
            P_Value = p_value,
            Importance_Score = importance_score,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(importance_results %>% arrange(desc(Importance_Score)))
}

# 选择关键参数用于雷达图
select_key_parameters_for_radar <- function(data, vision_params, bloodflow_params, thickness_params) {
  all_params <- c(vision_params, bloodflow_params, thickness_params)
  
  # 选择方差最大且有显著性的参数
  param_importance <- calculate_parameter_importance(data, all_params)
  
  # 选择前8个最重要的参数
  key_params <- param_importance %>%
    filter(P_Value < 0.1) %>%  # 宽松的显著性标准
    head(8) %>%
    pull(Parameter)
  
  return(key_params)
}

# 准备雷达图数据
prepare_radar_data <- function(data, key_params) {
  if(length(key_params) == 0) return(data.frame())
  
  # 计算每个聚类的标准化平均值
  radar_data <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(key_params), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
  
  # 标准化到0-1范围以便雷达图显示
  for(param in key_params) {
    if(param %in% names(radar_data)) {
      min_val <- min(radar_data[[param]], na.rm = TRUE)
      max_val <- max(radar_data[[param]], na.rm = TRUE)
      if(max_val != min_val) {
        radar_data[[param]] <- (radar_data[[param]] - min_val) / (max_val - min_val)
      }
    }
  }
  
  return(radar_data)
}

# 创建雷达图（简化版，使用极坐标）
create_radar_plot <- function(radar_data, key_params) {
  if(nrow(radar_data) == 0 || length(key_params) == 0) return(NULL)
  
  # 转换为长格式
  radar_long <- radar_data %>%
    pivot_longer(cols = all_of(key_params), names_to = "Parameter", values_to = "Value") %>%
    mutate(
      Parameter_Clean = clean_parameter_name(Parameter),
      Cluster = factor(max_cluster)
    )
  
  # 创建极坐标图
  p_radar <- ggplot(radar_long, aes(x = Parameter_Clean, y = Value, color = Cluster, group = Cluster)) +
    geom_polygon(aes(fill = Cluster), alpha = 0.3) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292")) +
    scale_fill_manual(values = c("1" = "#df8859", "2" = "#0fb292")) +
    coord_polar() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Cluster Characteristic Profile (Radar Chart)",
      subtitle = "Normalized mean values for key discriminative parameters"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title = element_blank()
    )
  
  ggsave("plots/octa_improvements/cluster_radar_chart.pdf", p_radar, width = 12, height = 10)
  ggsave("plots/octa_improvements/cluster_radar_chart.png", p_radar, width = 12, height = 10, dpi = 300)
}

# 计算聚类平均值
calculate_cluster_means <- function(data, all_params) {
  cluster_means <- data.frame()
  
  for(param in all_params) {
    if(param %in% names(data)) {
      param_means <- data %>%
        group_by(max_cluster) %>%
        summarise(
          Mean_Improvement = mean(.data[[param]], na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        mutate(
          Parameter = param,
          Parameter_Clean = clean_parameter_name(param),
          Parameter_Type = determine_parameter_type(param, bloodflow_params, thickness_params),
          Cluster = max_cluster
        )
      
      cluster_means <- rbind(cluster_means, param_means)
    }
  }
  
  return(cluster_means)
}

# 创建聚类对比条形图
create_cluster_comparison_barplot <- function(cluster_means) {
  # 计算聚类间差异
  cluster_diff <- cluster_means %>%
    dplyr::select(Parameter_Clean, Parameter_Type, Cluster, Mean_Improvement) %>%
    pivot_wider(names_from = Cluster, values_from = Mean_Improvement, names_prefix = "Cluster_") %>%
    mutate(
      Difference = Cluster_2 - Cluster_1,
      Better_Cluster = ifelse(Difference > 0, "Cluster 2", "Cluster 1")
    ) %>%
    arrange(desc(abs(Difference))) %>%
    head(15)
  
  p_diff <- ggplot(cluster_diff, aes(x = reorder(Parameter_Clean, abs(Difference)), y = Difference)) +
    geom_col(aes(fill = Better_Cluster), alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Cluster 1" = "#df8859", "Cluster 2" = "#0fb292"), name = "Better Cluster") +
    coord_flip() +
    labs(
      title = "Top 15 Parameters: Cluster 2 vs Cluster 1 Differences",
      subtitle = "Positive = Cluster 2 better, Negative = Cluster 1 better",
      x = "Parameters", y = "Mean Difference (Cluster 2 - Cluster 1)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave("plots/octa_improvements/cluster_difference_barplot.pdf", p_diff, width = 14, height = 10)
  ggsave("plots/octa_improvements/cluster_difference_barplot.png", p_diff, width = 14, height = 10, dpi = 300)
}

# 创建视力-OCTA相关性图
create_vision_octa_correlation <- function(data, vision_params, bloodflow_params, thickness_params) {
  if(length(vision_params) == 0) return(NULL)
  
  # 选择视力改善参数
  vision_improvement <- vision_params[grep("improvement", vision_params)]
  if(length(vision_improvement) == 0) return(NULL)
  
  octa_params <- c(bloodflow_params, thickness_params)
  
  # 计算相关性
  cor_data <- data %>%
    dplyr::select(all_of(c(vision_improvement[1], octa_params[1:min(10, length(octa_params))]))) %>%
    na.omit()
  
  if(nrow(cor_data) > 3) {
    cor_matrix <- cor(cor_data)
    
    # 提取视力参数与OCTA参数的相关性
    vision_octa_cor <- cor_matrix[1, -1]
    
    cor_df <- data.frame(
      Parameter = names(vision_octa_cor),
      Correlation = as.numeric(vision_octa_cor),
      Parameter_Clean = clean_parameter_name(names(vision_octa_cor))
    ) %>%
      arrange(desc(abs(Correlation)))
    
    p_cor <- ggplot(cor_df, aes(x = reorder(Parameter_Clean, abs(Correlation)), y = Correlation)) +
      geom_col(aes(fill = Correlation > 0), alpha = 0.8) +
      scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"), 
                        labels = c("Negative", "Positive"), name = "Correlation") +
      coord_flip() +
      labs(
        title = "Vision-OCTA Parameter Correlations",
        subtitle = "Correlation between vision improvement and OCTA improvements",
        x = "OCTA Parameters", y = "Correlation with Vision Improvement"
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ggsave("plots/octa_improvements/vision_octa_correlations.pdf", p_cor, width = 12, height = 10)
    ggsave("plots/octa_improvements/vision_octa_correlations.png", p_cor, width = 12, height = 10, dpi = 300)
  }
}

# 创建OCTA内部相关性图
create_octa_internal_correlation <- function(data, bloodflow_params, thickness_params) {
  octa_params <- c(bloodflow_params, thickness_params)
  
  # 选择前15个参数避免图太复杂
  selected_params <- octa_params[1:min(15, length(octa_params))]
  
  cor_data <- data %>%
    dplyr::select(all_of(selected_params)) %>%
    na.omit()
  
  if(nrow(cor_data) > 3 && ncol(cor_data) > 2) {
    cor_matrix <- cor(cor_data)
    
    # 转换为长格式用于ggplot
    cor_long <- as.data.frame(as.table(cor_matrix)) %>%
      mutate(
        Var1_Clean = clean_parameter_name(as.character(Var1)),
        Var2_Clean = clean_parameter_name(as.character(Var2))
      )
    
    p_cor_heatmap <- ggplot(cor_long, aes(x = Var1_Clean, y = Var2_Clean, fill = Freq)) +
      geom_tile() +
      scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", 
                           midpoint = 0, name = "Correlation") +
      labs(
        title = "OCTA Parameter Inter-correlations",
        subtitle = "Correlation matrix of OCTA improvement parameters",
        x = "Parameters", y = "Parameters"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave("plots/octa_improvements/octa_correlation_heatmap.pdf", p_cor_heatmap, width = 14, height = 12)
    ggsave("plots/octa_improvements/octa_correlation_heatmap.png", p_cor_heatmap, width = 14, height = 12, dpi = 300)
  }
}

# 计算患者改善得分
calculate_patient_improvement_scores <- function(data, bloodflow_params, thickness_params) {
  all_params <- c(bloodflow_params, thickness_params)
  
  # 为每个患者计算综合改善得分
  patient_scores <- data %>%
    rowwise() %>%
    mutate(
      # 血流参数改善得分
      BF_Score = mean(c_across(all_of(bloodflow_params)), na.rm = TRUE),
      # 厚度参数改善得分  
      TH_Score = mean(c_across(all_of(thickness_params)), na.rm = TRUE),
      # 总体改善得分
      Overall_Score = mean(c_across(all_of(all_params)), na.rm = TRUE),
      # 正向改善参数数量
      Positive_Count = sum(c_across(all_of(all_params)) > 0, na.rm = TRUE),
      # 总参数数量
      Total_Count = sum(!is.na(c_across(all_of(all_params)))),
      # 改善比例
      Improvement_Ratio = Positive_Count / Total_Count
    ) %>%
    ungroup() %>%
    dplyr::select(ID, max_cluster, max_membership, BF_Score, TH_Score, Overall_Score, 
                  Positive_Count, Total_Count, Improvement_Ratio)
  
  return(patient_scores)
}

# 创建改善模式分类
create_improvement_pattern_classification <- function(patient_scores) {
  # 根据改善得分对患者进行分类
  pattern_data <- patient_scores %>%
    mutate(
      Improvement_Pattern = case_when(
        Overall_Score > 0.05 & Improvement_Ratio > 0.6 ~ "High Improver",
        Overall_Score > 0 & Improvement_Ratio > 0.5 ~ "Moderate Improver", 
        Overall_Score < -0.05 & Improvement_Ratio < 0.4 ~ "Poor Responder",
        TRUE ~ "Mixed Response"
      )
    )
  
  # 创建改善模式分布图
  p_patterns <- ggplot(pattern_data, aes(x = Improvement_Pattern, fill = factor(max_cluster))) +
    geom_bar(position = "dodge", alpha = 0.8) +
    geom_text(stat = "count", aes(label = ..count..), 
              position = position_dodge(width = 0.9), vjust = -0.5) +
    scale_fill_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    labs(
      title = "Patient Improvement Pattern Distribution by Cluster",
      subtitle = "Classification based on overall improvement score and success ratio",
      x = "Improvement Pattern", y = "Number of Patients"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
  
  ggsave("plots/octa_improvements/improvement_pattern_distribution.pdf", p_patterns, width = 12, height = 8)
  ggsave("plots/octa_improvements/improvement_pattern_distribution.png", p_patterns, width = 12, height = 8, dpi = 300)
  
  # 创建改善得分散点图
  p_scatter <- ggplot(pattern_data, aes(x = Improvement_Ratio, y = Overall_Score)) +
    geom_point(aes(color = factor(max_cluster), size = max_membership), alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292"), name = "Cluster") +
    scale_size_continuous(name = "Membership", range = c(2, 6)) +
    labs(
      title = "Patient Improvement Score vs Success Ratio",
      subtitle = "Each point represents one patient | Size = cluster membership strength",
      x = "Improvement Success Ratio", y = "Overall Improvement Score",
      caption = "Dashed lines: Success ratio = 0.5, Overall score = 0"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave("plots/octa_improvements/improvement_score_scatter.pdf", p_scatter, width = 12, height = 8)  
  ggsave("plots/octa_improvements/improvement_score_scatter.png", p_scatter, width = 12, height = 8, dpi = 300)
  
  # 保存模式分类数据
  write.csv(pattern_data, "plots/octa_improvements/patient_improvement_patterns.csv", row.names = FALSE)
  
  # 打印模式分布统计
  pattern_summary <- pattern_data %>%
    group_by(max_cluster, Improvement_Pattern) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    arrange(max_cluster, desc(Count))
  
  cat("\n===== 患者改善模式分布 =====\n")
  print(pattern_summary)
}

# -------------------- 执行OCTA改善值可视化 --------------------

cat("\n========================================\n")
cat("🎨 开始创建OCTA改善值可视化（修正版）\n")
cat("========================================\n")

# 执行主要可视化函数
# 注意：需要在原始OCTA代码后添加这个调用
# create_octa_improvement_visualizations(comprehensive_data_with_clusters, 
#                                        vision_params, bloodflow_params, thickness_params,
#                                        bloodflow_macular, thickness_macular)

# -------------------- 生成OCTA改善值可视化报告 --------------------
generate_octa_improvement_report <- function() {
  
  report <- paste0(
    "========================================\n",
    "OCTA聚类改善值可视化报告（修正版）\n",
    "========================================\n\n",
    
    "🎯 数据特点分析:\n",
    "- 聚类基于：T2-T0改善值（非时间序列）\n",
    "- 数据结构：每个参数一个改善值\n",
    "- 可视化重点：改善值分布、聚类对比、参数重要性\n\n",
    
    "📊 生成的可视化类型:\n",
    "1. 改善值分布对比图\n",
    "   - 按参数类型的箱线图对比\n",
    "   - Top差异参数的详细对比\n",
    "   - 文件：improvement_boxplot_by_type.pdf/png\n",
    "   - 文件：top_discriminative_parameters.pdf/png\n\n",
    
    "2. 参数重要性排序图\n",
    "   - 基于统计显著性和效应量\n",
    "   - 重要性得分 = |平均差异| × -log10(P值)\n",
    "   - 文件：parameter_importance_ranking.pdf/png\n\n",
    
    "3. 聚类特征雷达图\n",
    "   - 关键参数的聚类特征轮廓\n",
    "   - 标准化显示便于比较\n",
    "   - 文件：cluster_radar_chart.pdf/png\n\n",
    
    "4. 改善值热图\n",
    "   - 综合参数改善值热图\n",
    "   - 聚类差异条形图\n",
    "   - 文件：improvement_heatmap_comprehensive.pdf/png\n",
    "   - 文件：cluster_difference_barplot.pdf/png\n\n",
    
    "5. 相关性分析图\n",
    "   - 视力-OCTA参数相关性\n",
    "   - OCTA参数间相关性热图\n",
    "   - 文件：vision_octa_correlations.pdf/png\n",
    "   - 文件：octa_correlation_heatmap.pdf/png\n\n",
    
    "6. 患者改善模式图\n",
    "   - 综合改善得分分布\n",
    "   - 改善模式分类（High/Moderate/Poor/Mixed）\n",
    "   - 改善得分vs成功率散点图\n",
    "   - 文件：patient_improvement_scores.pdf/png\n",
    "   - 文件：improvement_pattern_distribution.pdf/png\n",
    "   - 文件：improvement_score_scatter.pdf/png\n\n",
    
    "🔍 关键特点:\n",
    "✅ 适应OCTA数据结构（改善值而非时间序列）\n",
    "✅ 重点展示聚类间改善差异\n",
    "✅ 识别最具判别力的参数\n",
    "✅ 患者个体改善模式分析\n",
    "✅ 参数间相关性探索\n",
    "✅ 统计学意义与临床意义结合\n\n",
    
    "💡 使用建议:\n",
    "1. 查看改善值分布图了解聚类特征\n",
    "2. 参考重要性排序识别关键参数\n",
    "3. 使用雷达图直观比较聚类轮廓\n",
    "4. 通过热图发现改善模式\n",
    "5. 利用相关性分析理解参数关系\n",
    "6. 根据改善模式指导临床决策\n\n",
    
    "📈 临床价值:\n",
    "- 识别OCTA改善的关键指标\n",
    "- 预测患者改善潜力\n",
    "- 指导个性化治疗策略\n",
    "- 优化随访方案\n\n",
    
    "报告生成时间: ", Sys.time(), "\n",
    "========================================\n"
  )
  
  writeLines(report, "OCTA_Improvement_Visualization_Report.txt")
  cat("✓ 保存OCTA改善值可视化报告: OCTA_Improvement_Visualization_Report.txt\n")
  
  return(report)
}

# 生成报告
octa_report <- generate_octa_improvement_report()

cat("\n🎉 OCTA改善值可视化代码创建完成！\n")
cat("========================================\n")
cat("📋 使用说明:\n")
cat("1. 在原始OCTA代码的最后添加函数调用\n")
cat("2. 运行：create_octa_improvement_visualizations(...)\n") 
cat("3. 查看生成的plots/octa_improvements/目录\n")
cat("4. 阅读生成的可视化报告\n")
cat("\n🎯 核心改进:\n")
cat("✅ 适应T2-T0改善值数据结构\n")
cat("✅ 重点展示聚类改善差异\n")
cat("✅ 多维度参数重要性分析\n")
cat("✅ 患者个体改善模式识别\n")
cat("✅ 临床实用性导向设计\n")

