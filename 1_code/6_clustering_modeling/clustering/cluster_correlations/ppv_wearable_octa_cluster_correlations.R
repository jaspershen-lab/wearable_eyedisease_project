# Association Analysis Between Wearable Device Clusters and OCTA Improvement Clusters
# Based on the approach described in the referenced paper using contingency tables and Fisher's Exact test

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
octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_comprehensive_cluster_results.csv")

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
  scale_fill_gradient(low = "white", high = "steelblue") +
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

# Step 9: Analyze which specific metrics from wearable data might be associated with OCTA outcomes

# Prepare data for analysis
# Extract the original wearable metrics data for patients in each cluster
wearable_data <- read.csv("mfuzz_D_Surg1_8h_filtered.csv")  # Adjust filename as needed

# Join with cluster information
wearable_data_with_clusters <- wearable_data %>%
  inner_join(combined_clusters_complete, by = c("subject_id" = "subject_id_wearable"))

# Analyze differences in wearable metrics between different OCTA outcome groups
metrics_analysis <- data.frame()

for (metric in c("mean_rhr_1", "cv_rhr_1", "cv_bo", "steps_max", "steps_total", "deep_sleep", "total_sleep")) {
  for (day in 0:7) {  # Analyze first week post-surgery
    col_name <- paste0("day_", day, "_", metric)
    
    if (col_name %in% colnames(wearable_data_with_clusters)) {
      # Extract data for this metric and day
      data_for_test <- wearable_data_with_clusters %>%
        dplyr::select(all_of(col_name), octa_cluster) %>%
        filter(!is.na(!!sym(col_name)))
      
      # Skip if insufficient data
      if (nrow(data_for_test) < 5) next
      
      # Perform t-test between OCTA cluster groups
      t_test_result <- t.test(reformulate("octa_cluster", col_name), data = data_for_test)
      
      # Add results to metrics_analysis dataframe
      metrics_analysis <- rbind(metrics_analysis, data.frame(
        Metric = metric,
        Day = day,
        P_Value = t_test_result$p.value,
        Mean_Diff = diff(t_test_result$estimate),
        Significant = ifelse(t_test_result$p.value < 0.05, "Yes", "No")
      ))
    }
  }
}

# Apply FDR correction for multiple testing
metrics_analysis$P_Adjusted <- p.adjust(metrics_analysis$P_Value, method = "fdr")
metrics_analysis$Significant_Adjusted <- ifelse(metrics_analysis$P_Adjusted < 0.05, "Yes", "No")

# Sort by significance
metrics_analysis <- metrics_analysis %>% arrange(P_Value)

# Print significant associations
cat("\n===== SIGNIFICANT ASSOCIATIONS BETWEEN SPECIFIC METRICS AND OCTA OUTCOMES =====\n\n")
print(metrics_analysis %>% filter(Significant == "Yes"))

# Save all results
write.csv(metrics_analysis, "wearable_octa_metric_associations.csv", row.names = FALSE)
write.csv(combined_clusters_complete, "combined_cluster_assignments.csv", row.names = FALSE)

cat("\nAssociation analysis complete. Results and visualizations have been saved.\n")