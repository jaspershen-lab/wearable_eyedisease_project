# Association Analysis Between Wearable Device Clusters and OCTA Improvement Clusters
# Based on the approach described in the referenced paper using contingency tables and Fisher's Exact test
# 修正版本：使用正确的列名

# Load required libraries
library(tidyverse)
library(vcd)       # For additional association statistics
library(ggplot2)
library(gridExtra) # For arranging multiple plots
library(ggmosaic)

# Set working directory - update this to your project directory
# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# Step 1: Load the clustering results from both analyses
wearable_clusters  <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/ppv_cluster_results_time_windows.csv", check.names = FALSE)
# wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/ppv_cluster_results_all_metrics.csv")
# octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_comprehensive_cluster_results.csv")
octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/fcm_octa/WF_only_cluster/fcm_all_parameters_results.csv")

dir.create("3_data_analysis/6_clustering_modeling/cluster_association_analysis/ppv_octa/timepoint",
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/cluster_association_analysis/ppv_octa/timepoint")

# Print summary of the loaded data
cat("Wearable device clusters summary:\n")
print(summary(wearable_clusters))
cat("\nOCTA improvement clusters summary:\n")
print(summary(octa_clusters))

# Check column names to verify structure
cat("\nWearable clusters column names:", paste(colnames(wearable_clusters), collapse = ", "), "\n")
cat("OCTA clusters column names:", paste(colnames(octa_clusters), collapse = ", "), "\n")

# Step 1.5: Standardize column names for merging
# Rename Patient_ID to subject_id in octa_clusters for consistency
octa_clusters_renamed <- octa_clusters %>%
  rename(subject_id = Patient_ID,
         max_cluster = Hard_Cluster,
         max_membership = Max_Membership)

cat("\nAfter renaming - OCTA clusters column names:", paste(colnames(octa_clusters_renamed), collapse = ", "), "\n")

# Step 2: Merge the datasets by subject ID
# Using the standardized column names
combined_clusters <- wearable_clusters %>%
  inner_join(octa_clusters_renamed, by = "subject_id", suffix = c("_wearable", "_octa"))

cat("\nCombined dataset column names:", paste(colnames(combined_clusters), collapse = ", "), "\n")

# Remove rows where either clustering is missing
combined_clusters_complete <- combined_clusters[!is.na(combined_clusters$max_cluster_wearable) & 
                                                  !is.na(combined_clusters$max_cluster_octa), ]

# Create cleaner column names for analysis
combined_clusters_complete$wearable_cluster <- combined_clusters_complete$max_cluster_wearable
combined_clusters_complete$octa_cluster <- combined_clusters_complete$max_cluster_octa

# Print information about matched patients
cat("\nTotal patients with wearable device data:", nrow(wearable_clusters))
cat("\nTotal patients with OCTA data:", nrow(octa_clusters))
cat("\nPatients with both types of data:", nrow(combined_clusters_complete))

# Check the distribution of clusters
cat("\nWearable cluster distribution:\n")
print(table(combined_clusters_complete$wearable_cluster))
cat("\nOCTA cluster distribution:\n")
print(table(combined_clusters_complete$octa_cluster))

# Step 3: Create a contingency table
contingency_table <- table(
  Wearable_Cluster = combined_clusters_complete$wearable_cluster,
  OCTA_Cluster = combined_clusters_complete$octa_cluster
)

# Check if contingency table has sufficient data for analysis
if(any(dim(contingency_table) < 2)) {
  stop("Contingency table has insufficient dimensions for Fisher's test. Check your clustering results.")
}

if(sum(contingency_table) < 4) {
  stop("Insufficient total observations for reliable statistical analysis.")
}

# Step 4: Calculate expected frequencies and add row/column totals
row_totals <- rowSums(contingency_table)
col_totals <- colSums(contingency_table)
total_count <- sum(contingency_table)

expected_table <- outer(row_totals, col_totals) / total_count

# Step 5: Perform Fisher's Exact test for association
# Use appropriate method based on table size
if(min(dim(contingency_table)) <= 2 && sum(contingency_table) <= 100) {
  # Use exact method for small tables
  fisher_test_result <- fisher.test(contingency_table)
  cat("\nUsing exact Fisher's test (small table)\n")
} else {
  # Use simulation for larger tables
  fisher_test_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
  cat("\nUsing simulated Fisher's test (larger table)\n")
}

# Step 6: Calculate additional association statistics
tryCatch({
  assoc_stats <- assocstats(contingency_table)
}, error = function(e) {
  cat("Warning: Could not calculate association statistics:", e$message, "\n")
  assoc_stats <- NULL
})

# Calculate Cramér's V manually if assocstats fails
chi_test <- chisq.test(contingency_table, correct = FALSE)
cramers_v <- sqrt(chi_test$statistic / (total_count * (min(dim(contingency_table)) - 1)))

# Step 7: Print results
cat("\n===== ASSOCIATION ANALYSIS RESULTS =====\n\n")
cat("Observed Contingency Table:\n")
print(contingency_table)
cat("\nRow totals:", row_totals, "\n")
cat("Column totals:", col_totals, "\n")
cat("Total observations:", total_count, "\n")

cat("\nExpected Frequencies:\n")
print(round(expected_table, 2))

cat("\nFisher's Exact Test Result:\n")
print(fisher_test_result)

if(!is.null(assoc_stats)) {
  cat("\nAssociation Statistics:\n")
  print(assoc_stats)
} else {
  cat("\nCramér's V (manually calculated):", round(as.numeric(cramers_v), 4), "\n")
}

# Interpretation
cat("\n===== INTERPRETATION =====\n")
cat("Sample size:", total_count, "patients\n")
cat("Fisher's exact test p-value:", round(fisher_test_result$p.value, 4), "\n")

if(fisher_test_result$p.value < 0.05) {
  cat("Result: Significant association detected (p < 0.05)\n")
} else {
  cat("Result: No significant association detected (p >= 0.05)\n")
}

# Odds ratio interpretation (if 2x2 table)
if(all(dim(contingency_table) == 2)) {
  cat("Odds ratio:", round(fisher_test_result$estimate, 3), "\n")
  cat("95% CI for odds ratio: [", round(fisher_test_result$conf.int[1], 3), ",", 
      round(fisher_test_result$conf.int[2], 3), "]\n")
}

cat("Cramér's V (effect size):", round(as.numeric(cramers_v), 3), "\n")
if(cramers_v < 0.1) {
  cat("Effect size: Negligible\n")
} else if(cramers_v < 0.3) {
  cat("Effect size: Small\n")
} else if(cramers_v < 0.5) {
  cat("Effect size: Medium\n")
} else {
  cat("Effect size: Large\n")
}

# Step 8: Create visualizations

# Create a heatmap of the contingency table
contingency_df <- as.data.frame(as.table(contingency_table))
names(contingency_df) <- c("Wearable_Cluster", "OCTA_Cluster", "Frequency")

# Calculate percentages for cell labeling
contingency_df$Percentage <- contingency_df$Frequency / sum(contingency_df$Frequency) * 100

# Create the heatmap with improved styling
heatmap_plot <- ggplot(contingency_df, aes(x = OCTA_Cluster, y = Wearable_Cluster, fill = Frequency)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Frequency, Percentage)), 
            color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Count") +
  labs(title = "Association Between Wearable Device Clusters and OCTA Improvement Clusters",
       subtitle = paste("Fisher's exact test p-value:", round(fisher_test_result$p.value, 4),
                        "| Cramér's V:", round(as.numeric(cramers_v), 3),
                        "| n =", total_count),
       x = "OCTA Improvement Cluster",
       y = "Wearable Device Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"))

print(heatmap_plot)

# Create a mosaic plot for visualizing the association
tryCatch({
  mosaic_plot <- ggplot(contingency_df) +
    geom_mosaic(aes(weight = Frequency, x = product(Wearable_Cluster), fill = OCTA_Cluster)) +
    labs(title = "Mosaic Plot of Cluster Associations",
         subtitle = paste("Area represents frequency,", 
                          ifelse(fisher_test_result$p.value < 0.05, "Significant", "Non-significant"), 
                          "association"),
         fill = "OCTA Cluster") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11))
  
  print(mosaic_plot)
}, error = function(e) {
  cat("Could not create mosaic plot:", e$message, "\n")
  mosaic_plot <- NULL
})

# Save the plots
ggsave("wearable_octa_association_heatmap.pdf", heatmap_plot, width = 12, height = 8, dpi = 300)
cat("✓ Heatmap saved as: wearable_octa_association_heatmap.pdf\n")

if(!is.null(mosaic_plot)) {
  ggsave("wearable_octa_association_mosaic.pdf", mosaic_plot, width = 12, height = 8, dpi = 300)
  cat("✓ Mosaic plot saved as: wearable_octa_association_mosaic.pdf\n")
}

# Save results to CSV
results_summary <- data.frame(
  Metric = c("Total_Patients", "Fisher_P_Value", "Cramers_V", "Significant_Association"),
  Value = c(total_count, fisher_test_result$p.value, as.numeric(cramers_v), 
            fisher_test_result$p.value < 0.05),
  stringsAsFactors = FALSE
)

write.csv(contingency_table, "contingency_table.csv")
write.csv(results_summary, "association_analysis_summary.csv", row.names = FALSE)
write.csv(combined_clusters_complete, "combined_cluster_data.csv", row.names = FALSE)

cat("✓ Contingency table saved as: contingency_table.csv\n")
cat("✓ Analysis summary saved as: association_analysis_summary.csv\n")
cat("✓ Combined data saved as: combined_cluster_data.csv\n")

cat("\n===== ANALYSIS COMPLETE =====\n")