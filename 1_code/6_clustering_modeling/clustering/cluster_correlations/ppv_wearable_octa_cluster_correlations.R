# Association Analysis Between Wearable Device Clusters and OCTA Improvement Clusters
# Based on the approach described in the referenced paper using contingency tables and Fisher's Exact test

# Load required libraries
library(tidyverse)
library(vcd)       # For additional association statistics
library(ggplot2)
library(gridExtra) # For arranging multiple plots
library(ggmosaic)
library(r4projects)
# Set working directory - update this to your project directory
# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# Step 1: Load the clustering results from both analyses
wearable_clusters  <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/ppv_cluster_results_time_windows.csv", check.names = FALSE)
# wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/ppv_cluster_results_all_metrics.csv")
# octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_octa_0_21_only_cluster_results.csv")
octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results.csv")
# octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/fcm_octa/WF_only_cluster/fcm_all_parameters_results.csv")
 
dir.create("3_data_analysis/6_clustering_modeling/cluster_association_analysis/ppv_octa/timepoint",
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/cluster_association_analysis/ppv_octa/timepoint")

# dir.create("3_data_analysis/6_clustering_modeling/cluster_association_analysis/ppv_octa/all_days",
#            recursive = TRUE, showWarnings = FALSE)
# setwd("3_data_analysis/6_clustering_modeling/cluster_association_analysis/ppv_octa/all_days")

# Print summary of the loaded data
cat("Wearable device clusters summary:\n")
print(summary(wearable_clusters))
cat("\nOCTA improvement clusters summary:\n")
print(summary(octa_clusters))

# Step 2: Merge the datasets by subject ID
# Note: There might be naming convention differences, adjust as needed
combined_clusters <- wearable_clusters %>%
  left_join(octa_clusters, by = "subject_id")

# Remove rows where either clustering is missing
combined_clusters_complete <- combined_clusters[!is.na(combined_clusters$max_cluster.x) & 
                                                  !is.na(combined_clusters$max_cluster.y), ]
combined_clusters_complete$wearable_cluster <- combined_clusters_complete$max_cluster.x
combined_clusters_complete$octa_cluster <- combined_clusters_complete$max_cluster.y

# Print information about matched patients
cat("\nTotal patients with wearable device data:", nrow(wearable_clusters))
cat("\nTotal patients with OCTA data:", nrow(octa_clusters))
cat("\nPatients with both types of data:", nrow(combined_clusters_complete))

# Step 3: Create a contingency table
contingency_table <- table(
  Wearable_Cluster = combined_clusters_complete$wearable_cluster,
  OCTA_Cluster = combined_clusters_complete$octa_cluster
)

# Step 4: Calculate expected frequencies and add row/column totals
row_totals <- rowSums(contingency_table)
col_totals <- colSums(contingency_table)
total_count <- sum(contingency_table)

expected_table <- outer(row_totals, col_totals) / total_count

# Step 5: Perform Fisher's Exact test for association
fisher_test_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)

# Step 6: Calculate additional association statistics
assoc_stats <- assocstats(contingency_table)

# Step 7: Print results
cat("\n===== ASSOCIATION ANALYSIS RESULTS =====\n\n")
cat("Observed Contingency Table:\n")
print(contingency_table)
cat("\nExpected Frequencies:\n")
print(expected_table)
cat("\nFisher's Exact Test Result:\n")
print(fisher_test_result)
cat("\nAssociation Statistics:\n")
print(assoc_stats)

# Step 8: Create visualizations

# Create a heatmap of the contingency table
contingency_df <- as.data.frame(as.table(contingency_table))
names(contingency_df) <- c("Wearable_Cluster", "OCTA_Cluster", "Frequency")

# Calculate percentages for cell labeling
contingency_df$Percentage <- contingency_df$Frequency / sum(contingency_df$Frequency) * 100

# Create the heatmap
heatmap_plot <- ggplot(contingency_df, aes(x = OCTA_Cluster, y = Wearable_Cluster, fill = Frequency)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Frequency, Percentage)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "#a488bf") +
  labs(title = "Association Between Wearable Device Clusters and OCTA Improvement Clusters",
       x = "OCTA Improvement Cluster",
       y = "Wearable Device Cluster",
       fill = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"))
heatmap_plot
# Create a mosaic plot for visualizing the association
mosaic_plot <- ggplot(contingency_df) +
  geom_mosaic(aes(weight = Frequency, x = product(Wearable_Cluster), fill = OCTA_Cluster)) +
  labs(title = "Mosaic Plot of Cluster Associations",
       fill = "OCTA Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Save the plots
ggsave("wearable_octa_association_heatmap.pdf", heatmap_plot, width = 10, height = 8)
ggsave("wearable_octa_association_mosaic.pdf", mosaic_plot, width = 10, height = 8)





# ===== SENSITIVITY ANALYSIS: EXCLUDING SH035 =====
# 敏感性分析：排除SH035样本后重新分析

cat("\n\n===== SENSITIVITY ANALYSIS (EXCLUDING SH035) =====\n\n")

# Step 1: Create dataset excluding SH035
combined_clusters_sensitivity <- combined_clusters_complete[combined_clusters_complete$subject_id != "SH035", ]

# Print information about the sensitivity analysis dataset
cat("Original dataset with both types of data:", nrow(combined_clusters_complete))
cat("\nSensitivity analysis dataset (excluding SH035):", nrow(combined_clusters_sensitivity))
cat("\nExcluded patients:", nrow(combined_clusters_complete) - nrow(combined_clusters_sensitivity))

# Check if SH035 was actually in the dataset
if ("SH035" %in% combined_clusters_complete$subject_id) {
  cat("\nSH035 was successfully excluded from the analysis")
} else {
  cat("\nSH035 was not found in the original dataset")
}

# Step 2: Create new contingency table for sensitivity analysis
contingency_table_sensitivity <- table(
  Wearable_Cluster = combined_clusters_sensitivity$wearable_cluster,
  OCTA_Cluster = combined_clusters_sensitivity$octa_cluster
)

# Step 3: Calculate expected frequencies for sensitivity analysis
row_totals_sensitivity <- rowSums(contingency_table_sensitivity)
col_totals_sensitivity <- colSums(contingency_table_sensitivity)
total_count_sensitivity <- sum(contingency_table_sensitivity)

expected_table_sensitivity <- outer(row_totals_sensitivity, col_totals_sensitivity) / total_count_sensitivity

# Step 4: Perform Fisher's Exact test for sensitivity analysis
fisher_test_sensitivity <- fisher.test(contingency_table_sensitivity, simulate.p.value = TRUE, B = 10000)

# Step 5: Calculate additional association statistics for sensitivity analysis
assoc_stats_sensitivity <- assocstats(contingency_table_sensitivity)

# Step 6: Print sensitivity analysis results
cat("\n===== SENSITIVITY ANALYSIS RESULTS =====\n\n")
cat("Observed Contingency Table (excluding SH035):\n")
print(contingency_table_sensitivity)
cat("\nExpected Frequencies (excluding SH035):\n")
print(expected_table_sensitivity)
cat("\nFisher's Exact Test Result (excluding SH035):\n")
print(fisher_test_sensitivity)
cat("\nAssociation Statistics (excluding SH035):\n")
print(assoc_stats_sensitivity)

# Step 7: Compare results between original and sensitivity analysis
cat("\n===== COMPARISON OF RESULTS =====\n\n")
cat("Original Analysis:\n")
cat("- Sample size:", nrow(combined_clusters_complete))
cat("\n- Fisher's Exact Test p-value:", format(fisher_test_result$p.value, scientific = TRUE, digits = 4))
cat("\n- Cramer's V:", format(assoc_stats$cramer, digits = 4))

cat("\n\nSensitivity Analysis (excluding SH035):\n")
cat("- Sample size:", nrow(combined_clusters_sensitivity))
cat("\n- Fisher's Exact Test p-value:", format(fisher_test_sensitivity$p.value, scientific = TRUE, digits = 4))
cat("\n- Cramer's V:", format(assoc_stats_sensitivity$cramer, digits = 4))

# Calculate the difference in p-values and effect sizes
p_value_diff <- abs(fisher_test_result$p.value - fisher_test_sensitivity$p.value)
cramer_v_diff <- abs(assoc_stats$cramer - assoc_stats_sensitivity$cramer)

cat("\n\nDifferences:\n")
cat("- P-value difference:", format(p_value_diff, scientific = TRUE, digits = 4))
cat("\n- Cramer's V difference:", format(cramer_v_diff, digits = 4))

# Step 8: Create visualizations for sensitivity analysis

# Create heatmap for sensitivity analysis
contingency_df_sensitivity <- as.data.frame(as.table(contingency_table_sensitivity))
names(contingency_df_sensitivity) <- c("Wearable_Cluster", "OCTA_Cluster", "Frequency")

# Calculate percentages for cell labeling
contingency_df_sensitivity$Percentage <- contingency_df_sensitivity$Frequency / sum(contingency_df_sensitivity$Frequency) * 100

# Create the heatmap for sensitivity analysis
heatmap_plot_sensitivity <- ggplot(contingency_df_sensitivity, aes(x = OCTA_Cluster, y = Wearable_Cluster, fill = Frequency)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Frequency, Percentage)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "#a488bf") +
  labs(title = "Sensitivity Analysis: Association Between Clusters (Excluding SH035)",
       x = "OCTA Improvement Cluster",
       y = "Wearable Device Cluster",
       fill = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"))

# Create mosaic plot for sensitivity analysis
mosaic_plot_sensitivity <- ggplot(contingency_df_sensitivity) +
  geom_mosaic(aes(weight = Frequency, x = product(Wearable_Cluster), fill = OCTA_Cluster)) +
  labs(title = "Sensitivity Analysis: Mosaic Plot of Cluster Associations (Excluding SH035)",
       fill = "OCTA Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Create comparison plots side by side
comparison_heatmap <- grid.arrange(
  heatmap_plot + labs(title = "Original Analysis\n(All Patients)"),
  heatmap_plot_sensitivity + labs(title = "Sensitivity Analysis\n(Excluding SH035)"),
  ncol = 2
)

# Step 9: Save sensitivity analysis results and plots
ggsave("sensitivity_analysis_heatmap.pdf", heatmap_plot_sensitivity, width = 10, height = 8)
ggsave("sensitivity_analysis_mosaic.pdf", mosaic_plot_sensitivity, width = 10, height = 8)
ggsave("comparison_heatmaps.pdf", comparison_heatmap, width = 16, height = 8)

# Step 10: Generate summary conclusion
cat("\n===== SENSITIVITY ANALYSIS CONCLUSION =====\n\n")

# Determine if results are robust
p_value_threshold <- 0.001  # Consider results robust if p-value difference < 0.001
effect_size_threshold <- 0.1  # Consider results robust if Cramer's V difference < 0.1

if (p_value_diff < p_value_threshold && cramer_v_diff < effect_size_threshold) {
  cat("CONCLUSION: The results are ROBUST to the exclusion of SH035.\n")
  cat("Both the statistical significance and effect size remain consistent.\n")
} else if (p_value_diff >= p_value_threshold) {
  cat("CONCLUSION: The statistical significance shows notable change after excluding SH035.\n")
  cat("This suggests that SH035 may have some influence on the overall association.\n")
} else if (cramer_v_diff >= effect_size_threshold) {
  cat("CONCLUSION: The effect size shows notable change after excluding SH035.\n")
  cat("This suggests that SH035 may influence the strength of association.\n")
}

cat("\nThis sensitivity analysis demonstrates whether the association between\n")
cat("wearable device clusters and OCTA improvement clusters depends on the\n")
cat("inclusion of the single PPV+IVI patient (SH035) versus the IVI-only patients.\n")

# Save all results to a text file
sink("sensitivity_analysis_results.txt")
cat("===== SENSITIVITY ANALYSIS RESULTS =====\n\n")
cat("Dataset Information:\n")
cat("- Original sample size:", nrow(combined_clusters_complete), "\n")
cat("- Sensitivity analysis sample size:", nrow(combined_clusters_sensitivity), "\n")
cat("- Excluded subject: SH035 (PPV+IVI patient)\n\n")

cat("Original Analysis Results:\n")
print(contingency_table)
print(fisher_test_result)
print(assoc_stats)

cat("\n\nSensitivity Analysis Results (excluding SH035):\n")
print(contingency_table_sensitivity)
print(fisher_test_sensitivity)
print(assoc_stats_sensitivity)

cat("\n\nComparison:\n")
cat("P-value difference:", format(p_value_diff, scientific = TRUE, digits = 4), "\n")
cat("Cramer's V difference:", format(cramer_v_diff, digits = 4), "\n")

if (p_value_diff < p_value_threshold && cramer_v_diff < effect_size_threshold) {
  cat("\nCONCLUSION: Results are robust to exclusion of SH035.\n")
} else {
  cat("\nCONCLUSION: Results show some sensitivity to exclusion of SH035.\n")
}
sink()

