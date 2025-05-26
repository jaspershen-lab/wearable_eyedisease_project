# ----------------------------------------------------
# Vision Improvement Clustering Analysis for PPV Group
# Using Mfuzz Fuzzy Clustering to Identify Patient Subgroups
# Focusing on 1-month vision improvement post-surgery
# ----------------------------------------------------

# Load required libraries
library(tidyverse)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(factoextra)
library(corrplot)
library(r4projects)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# -------------------- 1. Load data --------------------
# Load baseline information
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/mfuzz/vision_cluster", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/vision_cluster")

# -------------------- 2. Process vision data --------------------
# Create vision dataset with pre-surgery and post-surgery data
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,  # Right eye surgery
      surgery_eye_1 == 1 ~ os_corrected_bas,  # Left eye surgery
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # Both eyes (average)
      TRUE ~ NA_real_
    ),
    post_vision_1m = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1m,   # Right eye post-surgery 1 month
      surgery_eye_1 == 1 ~ os_corrected_1m,   # Left eye post-surgery 1 month
      surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,   # Both eyes post-surgery 1 month
      TRUE ~ NA_real_
    ),
    vision_improvement_1m = post_vision_1m - pre_vision
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision_1m, vision_improvement_1m, age, gender)

# -------------------- 3. Extract PPV group data --------------------
# Get PPV patient IDs (surgery type 1)
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

# Filter vision data for PPV patients
ppv_vision <- vision_data %>%
  filter(ID %in% ppv_patients)

# -------------------- 4. Prepare clustering data --------------------
# Save original ID list for later matching
original_ids <- ppv_vision$ID

# Create dataset for clustering with relevant features
ppv_features_for_clustering <- ppv_vision %>%
  dplyr::select(ID, vision_improvement_1m, pre_vision, age)

# Detailed analysis of missing values
na_count_by_row <- rowSums(is.na(ppv_features_for_clustering[,-1]))  # Exclude ID column
na_count_table <- table(na_count_by_row)

# Calculate missing values by column
na_count_by_col <- colSums(is.na(ppv_features_for_clustering[,-1]))

# Print missing values analysis
cat("\n===== Vision data missing values analysis =====\n")
cat("Total rows (patients):", nrow(ppv_features_for_clustering), "\n")
cat("Total columns (parameters):", ncol(ppv_features_for_clustering) - 1, "\n")  # Exclude ID

cat("\nPatient missing value distribution:\n")
for(na_count in names(na_count_table)) {
  cat("Patients with", na_count, "missing values:", na_count_table[na_count], "\n")
}

cat("\nMissing values by parameter:\n")
for(i in 2:length(ppv_features_for_clustering)) {
  col_name <- names(ppv_features_for_clustering)[i]
  col_na_count <- sum(is.na(ppv_features_for_clustering[[col_name]]))
  cat(col_name, ":", col_na_count, "missing values (", 
      round(col_na_count/nrow(ppv_features_for_clustering)*100, 1), "%)\n")
}

# Filter for complete data (no missing values)
complete_rows <- complete.cases(ppv_features_for_clustering[,-1])  # Exclude ID column
complete_data <- ppv_features_for_clustering[complete_rows, -1]  # Exclude ID for clustering
complete_ids <- original_ids[complete_rows]

# Print complete data statistics
cat("\n===== Complete data statistics =====\n")
cat("Original patient count:", nrow(ppv_features_for_clustering), "\n")
cat("Complete data patient count:", nrow(complete_data), 
    "(", round(nrow(complete_data)/nrow(ppv_features_for_clustering)*100, 1), "%)\n")
cat("Removed", sum(!complete_rows), "patients with missing values\n\n")

# Check if final dataset is suitable for clustering
if(nrow(complete_data) < 5) {
  stop("Not enough patients with complete data (less than 5). Consider using mean imputation or other strategies.")
} else if(nrow(complete_data) < 10) {
  cat("⚠️ Warning: Small number of patients with complete data (", nrow(complete_data), 
      "). Clustering reliability may be affected.\n\n")
}

# Standardize data for clustering
ppv_vision_std <- as.data.frame(scale(complete_data))

# Print processing summary
cat("===== Clustering data preparation complete =====\n")
cat("Number of parameters used:", ncol(ppv_vision_std), "\n")
cat("Parameters used:", paste(names(ppv_vision_std), collapse=", "), "\n")
cat("Number of patients used:", nrow(ppv_vision_std), "\n")
cat("Patient IDs used:", paste(complete_ids, collapse=", "), "\n\n")

# -------------------- 5. Perform Mfuzz clustering --------------------
# Prepare data for Mfuzz
prep_data_for_mfuzz <- function(data, row_ids) {
  data_matrix <- as.matrix(data)
  rownames(data_matrix) <- row_ids
  
  eset <- ExpressionSet(assayData = data_matrix)
  return(eset)
}

# Create ExpressionSet using complete data
ppv_vision_eset <- prep_data_for_mfuzz(ppv_vision_std, complete_ids)

# Standardize for Mfuzz
ppv_vision_eset_std <- standardise(ppv_vision_eset)

# Estimate optimal fuzzy coefficient m
ppv_m <- mestimate(ppv_vision_eset_std)
cat("Estimated optimal fuzzy coefficient m:", ppv_m, "\n")

# Perform clustering with 2 clusters
set.seed(123)
ppv_clustering <- mfuzz(ppv_vision_eset_std, c = 2, m = ppv_m)

# Get membership values and primary cluster assignments
ppv_membership <- ppv_clustering$membership
ppv_main_clusters <- apply(ppv_membership, 1, which.max)

# Create results dataframe
ppv_clusters_result <- data.frame(
  subject_id = complete_ids,
  max_cluster = ppv_main_clusters,
  max_membership = apply(ppv_membership, 1, max)
)

# Save clustering results
write.csv(ppv_clusters_result, "ppv_vision_cluster_results.csv", row.names = FALSE)

# -------------------- 6. Visualize clustering results --------------------
# Create combined dataset with original data and cluster assignments
ppv_vision_with_clusters <- ppv_vision %>%
  inner_join(ppv_clusters_result, by = c("ID" = "subject_id"))

# Function to visualize clusters
visualize_clusters <- function(data, params) {
  # Create plots directory
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  # Filter out any rows with NA in max_cluster
  data <- data %>% filter(!is.na(max_cluster))
  
  # Convert cluster to factor
  data$max_cluster <- as.factor(data$max_cluster)
  
  # Calculate mean values for each cluster
  cluster_means <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(params), mean, na.rm = TRUE))
  
  # Print cluster means
  cat("\nCluster means for key parameters:\n")
  print(cluster_means)
  
  # Create boxplot for vision improvement by cluster
  p_vision_improvement <- ggplot(data, aes(x = max_cluster, y = vision_improvement_1m, fill = max_cluster)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, width = 0.2) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      title = "1-Month Vision Improvement by Cluster",
      x = "Cluster",
      y = "Vision Improvement (Post - Pre)"
    )
  
  # Save vision improvement boxplot
  ggsave("plots/vision_improvement_boxplot.pdf", p_vision_improvement, width = 8, height = 6)
  ggsave("plots/vision_improvement_boxplot.png", p_vision_improvement, width = 8, height = 6, dpi = 300)
  
  # Create boxplot for pre-vision by cluster
  p_pre_vision <- ggplot(data, aes(x = max_cluster, y = pre_vision, fill = max_cluster)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, width = 0.2) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      title = "Pre-Surgery Vision by Cluster",
      x = "Cluster",
      y = "Pre-Surgery Vision"
    )
  
  # Save pre-vision boxplot
  ggsave("plots/pre_vision_boxplot.pdf", p_pre_vision, width = 8, height = 6)
  ggsave("plots/pre_vision_boxplot.png", p_pre_vision, width = 8, height = 6, dpi = 300)
  
  # Create boxplot for post-vision by cluster
  p_post_vision <- ggplot(data, aes(x = max_cluster, y = post_vision_1m, fill = max_cluster)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, width = 0.2) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      title = "Post-Surgery Vision (1 Month) by Cluster",
      x = "Cluster",
      y = "Post-Surgery Vision (1 Month)"
    )
  
  # Save post-vision boxplot
  ggsave("plots/post_vision_boxplot.pdf", p_post_vision, width = 8, height = 6)
  ggsave("plots/post_vision_boxplot.png", p_post_vision, width = 8, height = 6, dpi = 300)
  
  # Create scatter plot of pre vs. post vision colored by cluster
  p_scatter <- ggplot(data, aes(x = pre_vision, y = post_vision_1m, color = max_cluster)) +
    geom_point(aes(size = max_membership), alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      title = "Pre vs. Post Vision by Cluster",
      x = "Pre-Surgery Vision",
      y = "Post-Surgery Vision (1 Month)",
      color = "Cluster",
      size = "Membership"
    )
  
  # Save scatter plot
  ggsave("plots/pre_post_vision_scatter.pdf", p_scatter, width = 10, height = 8)
  ggsave("plots/pre_post_vision_scatter.png", p_scatter, width = 10, height = 8, dpi = 300)
  
  # Create PCA visualization if there are enough parameters
  if(length(params) > 1) {
    # Perform PCA
    pca_data <- data %>%
      dplyr::select(all_of(params))
    
    pca_result <- prcomp(pca_data, scale. = TRUE)
    
    # Prepare plot data
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
        title = "PCA of Vision Parameters",
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    # Save PCA plot
    ggsave("plots/vision_pca_plot.pdf", p_pca, width = 10, height = 8)
    ggsave("plots/vision_pca_plot.png", p_pca, width = 10, height = 8, dpi = 300)
  }
}

# Get parameters for visualization
vision_params <- c("vision_improvement_1m", "pre_vision", "post_vision_1m", "age")
vision_params <- vision_params[vision_params %in% names(ppv_vision_with_clusters)]

# Create visualizations
visualize_clusters(ppv_vision_with_clusters, vision_params)

# -------------------- 7. Statistical analysis of cluster differences --------------------
# Function to perform statistical analysis between clusters
analyze_cluster_differences <- function(data, params) {
  results <- data.frame(
    Parameter = character(),
    Cluster1_Mean = numeric(),
    Cluster2_Mean = numeric(),
    Mean_Difference = numeric(),
    P_Value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )
  
  for(param in params) {
    # Extract data for this parameter
    param_data <- data[, c("max_cluster", param)]
    
    # Calculate mean for each cluster
    means <- tapply(param_data[[param]], param_data$max_cluster, mean, na.rm = TRUE)
    
    # Perform t-test
    test_result <- t.test(reformulate("max_cluster", param), data = param_data)
    
    # Add to results
    results <- rbind(results, data.frame(
      Parameter = param,
      Cluster1_Mean = means["1"],
      Cluster2_Mean = means["2"],
      Mean_Difference = means["2"] - means["1"],
      P_Value = test_result$p.value,
      Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Add adjusted p-values
  results$P_Adjusted <- p.adjust(results$P_Value, method = "fdr")
  results$Significant_Adjusted <- ifelse(results$P_Adjusted < 0.05, "Yes", "No")
  
  # Sort by p-value
  results <- results %>% arrange(P_Value)
  
  return(results)
}

# Perform statistical analysis
ppv_stats <- analyze_cluster_differences(ppv_vision_with_clusters, vision_params)

# Save statistical results
write.csv(ppv_stats, "ppv_vision_cluster_statistics.csv", row.names = FALSE)

# Print significant parameters
cat("\nSignificant parameters differentiating clusters:\n")
print(ppv_stats[ppv_stats$Significant == "Yes", ])

# -------------------- 8. Interpret clusters --------------------
# Function to interpret which cluster shows better outcomes
interpret_clusters <- function(stats_results) {
  # Focus on vision improvement specifically
  vision_improvement_result <- stats_results %>%
    filter(Parameter == "vision_improvement_1m")
  
  if(nrow(vision_improvement_result) == 0) {
    cat("\nError: No vision improvement data found in statistical results\n")
    return(NULL)
  }
  
  # Determine which cluster has better vision improvement
  if(vision_improvement_result$Mean_Difference > 0) {
    better_cluster <- 2
    worse_cluster <- 1
  } else {
    better_cluster <- 1
    worse_cluster <- 2
  }
  
  # Check if the difference is significant
  is_significant <- vision_improvement_result$Significant == "Yes"
  
  cat("\nCluster interpretation:\n")
  cat("Cluster", better_cluster, "shows better vision improvement\n")
  cat("Cluster", worse_cluster, "shows worse vision improvement\n")
  cat("The difference is", ifelse(is_significant, "statistically significant", "not statistically significant"), "\n")
  
  return(list(
    better_cluster = better_cluster,
    worse_cluster = worse_cluster,
    is_significant = is_significant,
    mean_difference = abs(vision_improvement_result$Mean_Difference)
  ))
}

# Interpret clusters
cluster_interpretation <- interpret_clusters(ppv_stats)

# Add outcome labels to cluster results
if(!is.null(cluster_interpretation)) {
  ppv_clusters_result$outcome_quality <- ifelse(
    ppv_clusters_result$max_cluster == cluster_interpretation$better_cluster,
    "Better",
    "Worse"
  )
  
  # Save final results with outcome labels
  write.csv(ppv_clusters_result, "ppv_vision_cluster_results_with_outcomes.csv", row.names = FALSE)
}

# -------------------- 9. Cluster summary statistics --------------------
# Calculate summary statistics for each cluster
cluster_summary <- ppv_vision_with_clusters %>%
  group_by(max_cluster) %>%
  summarise(
    Count = n(),
    Mean_Membership = mean(max_membership, na.rm = TRUE),
    Min_Membership = min(max_membership, na.rm = TRUE),
    Max_Membership = max(max_membership, na.rm = TRUE),
    Mean_Vision_Improvement = mean(vision_improvement_1m, na.rm = TRUE),
    Mean_Pre_Vision = mean(pre_vision, na.rm = TRUE),
    Mean_Post_Vision = mean(post_vision_1m, na.rm = TRUE),
    Mean_Age = mean(age, na.rm = TRUE)
  ) %>%
  mutate(
    Outcome_Quality = ifelse(max_cluster == cluster_interpretation$better_cluster, "Better", "Worse")
  )

# Print cluster summary
cat("\nCluster summary:\n")
print(cluster_summary)

# Save cluster summary
write.csv(cluster_summary, "ppv_vision_cluster_summary.csv", row.names = FALSE)

# Print final message
cat("\nPPV group vision clustering analysis complete.\n")
cat("Clustering results saved to 'ppv_vision_cluster_results_with_outcomes.csv'\n")
cat("Cluster statistics saved to 'ppv_vision_cluster_statistics.csv'\n")
cat("Visualizations saved to 'plots' directory\n")
