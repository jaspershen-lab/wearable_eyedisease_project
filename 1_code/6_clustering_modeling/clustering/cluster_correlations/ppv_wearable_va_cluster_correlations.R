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
# wearable_clusters  <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/ppv_cluster_results_time_windows.csv", check.names = FALSE)
wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/ppv_cluster_results_all_metrics.csv")
va_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/vision_cluster/ppv_vision_cluster_results.csv")

# Print summary of the loaded data
cat("Wearable device clusters summary:\n")
print(summary(wearable_clusters))
cat("\nVA improvement clusters summary:\n")
print(summary(va_clusters))

# Step 2: Merge the datasets by subject ID
# Note: There might be naming convention differences, adjust as needed
combined_clusters <- wearable_clusters %>%
  left_join(va_clusters, by = "subject_id")

# Remove rows where either clustering is missing
combined_clusters_complete <- combined_clusters[!is.na(combined_clusters$max_cluster.x) & 
                                                  !is.na(combined_clusters$max_cluster.y), ]
combined_clusters_complete$wearable_cluster <- combined_clusters_complete$max_cluster.x
combined_clusters_complete$va_cluster <- combined_clusters_complete$max_cluster.y

# Print information about matched patients
cat("\nTotal patients with wearable device data:", nrow(wearable_clusters))
cat("\nTotal patients with VA data:", nrow(va_clusters))
cat("\nPatients with both types of data:", nrow(combined_clusters_complete))

# Step 3: Create a contingency table
contingency_table <- table(
  Wearable_Cluster = combined_clusters_complete$wearable_cluster,
  VA_cluster = combined_clusters_complete$va_cluster
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
