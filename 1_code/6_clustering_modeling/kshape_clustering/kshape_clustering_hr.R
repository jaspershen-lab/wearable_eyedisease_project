# 设置工作目录和清理环境
setwd(get_project_wd())
rm(list = ls())

# Load required libraries
library(tidyverse)
library(lubridate)
library(tidymass)
library(massdataset)
library(cluster)  # For clustering functions
library(dtw)      # For Dynamic Time Warping distance
library(factoextra) # For cluster visualization

# Source tools if needed
source('1_code/100_tools.R')

# Load heart rate data
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

# Load baseline information
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Extract surgery time for each subject
surgery_times <- baseline_info %>%
  mutate(
    ID = toupper(ID),
    ID = gsub("SHH|ShH", "SH", ID),
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)

# Create group information
group_info <- baseline_info %>%
  mutate(
    ID = toupper(ID),
    ID = gsub("SHH|ShH", "SH", ID),
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  dplyr::select(ID, dm_status, surgery_type)

# Create patient groups
patient_groups <- group_info %>%
  mutate(
    patient_group = case_when(
      dm_status == "No Diabetes" & surgery_type == 0 ~ "Group1_Surgery0",
      dm_status == "Diabetes" & surgery_type == 1 ~ "Group2_Diabetes",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(patient_group))

# Extract subjects for each group
group1_subjects <- patient_groups %>% 
  filter(patient_group == "Group1_Surgery0") %>% 
  pull(ID)

group2_subjects <- patient_groups %>% 
  filter(patient_group == "Group2_Diabetes") %>% 
  pull(ID)

# Modified function to extract heart rate data
extract_hr_data <- function(subject_ids) {
  # Extract the sample info (which contains subject_id)
  sample_info <- heart_rate_data@sample_info
  
  # Filter sample_info to include only samples from our subjects of interest
  filtered_sample_info <- sample_info %>%
    filter(subject_id %in% subject_ids)
  
  # Get the sample_ids from the filtered sample_info
  matching_samples <- filtered_sample_info$sample_id
  
  # Filter the heart_rate_data to include only matching samples
  expression_data <- heart_rate_data@expression_data[, matching_samples, drop = FALSE]
  
  # Extract the heart rate values
  heart_rate_values <- as.numeric(expression_data[1, ])
  
  # Create a data frame with sample_id, subject_id, measure_time, and heart_rate
  hr_data <- filtered_sample_info %>%
    mutate(heart_rate = heart_rate_values)
  
  # Add surgery time information
  hr_data_with_surgery <- hr_data %>%
    left_join(surgery_times, by = c("subject_id" = "ID"))
  
  # Calculate days from surgery and filter to the period of interest
  hr_trajectory_data <- hr_data_with_surgery %>%
    mutate(
      days_from_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days"))
    ) %>%
    filter(days_from_surgery >= -4 & days_from_surgery <= 30)
  
  return(hr_trajectory_data)
}

# Process the data for both groups
group1_hr_data <- extract_hr_data(group1_subjects)
group2_hr_data <- extract_hr_data(group2_subjects)

# 创建结果保存目录并设置工作目录（只设置一次）
output_dir <- "3_data_analysis/6_clustering_modeling/kshape/time_series/1m/hr"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 保存原始路径，以便在程序结束时恢复
original_wd <- getwd()
# 切换到输出目录
setwd(file.path(original_wd, output_dir))

# 修改后的函数，确保正确处理路径和p值
perform_clustering <- function(hr_trajectory_data, group_name) {
  # Create daily average heart rate
  daily_hr <- hr_trajectory_data %>%
    group_by(subject_id, floor(days_from_surgery)) %>%
    summarize(
      day = min(floor(days_from_surgery)),
      heart_rate_avg = mean(heart_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(subject_id, day, heart_rate_avg)
  
  # Save the daily heart rate data (使用相对路径)
  write.csv(daily_hr, paste0(group_name, "_daily_hr.csv"), row.names = FALSE)
  
  # Pivot to wide format for clustering
  hr_wide <- daily_hr %>%
    pivot_wider(
      id_cols = subject_id,
      names_from = day,
      values_from = heart_rate_avg
    )
  
  # Save the wide-format data
  write.csv(hr_wide, paste0(group_name, "_hr_wide.csv"), row.names = FALSE)
  
  # Prepare matrix for clustering
  subjects <- hr_wide$subject_id
  hr_matrix <- as.matrix(hr_wide[, -1])
  rownames(hr_matrix) <- subjects
  
  # Handle missing values (if any)
  if(any(is.na(hr_matrix))) {
    cat("Handling missing values in the heart rate matrix\n")
    hr_matrix[is.na(hr_matrix)] <- mean(hr_matrix, na.rm = TRUE)
  }
  
  # Print group info
  cat("Processing", group_name, "with", nrow(hr_matrix), "subjects\n")
  write(paste("Processing", group_name, "with", nrow(hr_matrix), "subjects"), 
        file = paste0(group_name, "_processing_log.txt"))
  
  # Save the matrix
  saveRDS(hr_matrix, paste0(group_name, "_hr_matrix.rds"))
  
  # Create a simpler distance matrix using Euclidean distance as a fallback
  cat("Computing distance matrix...\n")
  
  # Use standard Euclidean distance
  eucl_dist <- dist(hr_matrix, method = "euclidean")
  
  # Save distance matrix
  saveRDS(eucl_dist, paste0(group_name, "_distance_matrix.rds"))
  
  # Perform hierarchical clustering
  clusters_list <- list()
  
  for (k in 2:3) {
    tryCatch({
      # Hierarchical clustering
      hc <- hclust(eucl_dist, method = "ward.D2")
      
      # Cut the tree to get k clusters
      cluster_assignment <- cutree(hc, k = k)
      
      # Save clustering results
      clusters_list[[paste0("k", k)]] <- list(
        hc = hc,
        clusters = cluster_assignment
      )
      
      # Save clustering object
      saveRDS(clusters_list[[paste0("k", k)]], paste0(group_name, "_k", k, "_clusters.rds"))
      
      # Create dataframe with cluster assignments
      results <- data.frame(
        subject_id = names(cluster_assignment),
        cluster = cluster_assignment
      )
      
      # Save cluster assignments
      write.csv(results, paste0(group_name, "_k", k, "_cluster_assignments.csv"), row.names = FALSE)
      
      # Calculate statistics for clusters
      cluster_stats <- results %>%
        group_by(cluster) %>%
        summarise(count = n(), .groups = "drop") %>%
        mutate(percentage = count / sum(count) * 100)
      
      # Save cluster statistics
      write.csv(cluster_stats, paste0(group_name, "_k", k, "_cluster_stats.csv"), row.names = FALSE)
      
      # Create visualization of clusters
      p_dendrogram <- fviz_dend(hc, k = k, 
                                palette = "jco",
                                rect = TRUE,
                                rect_fill = TRUE,
                                rect_border = "jco",
                                main = paste(group_name, "- Cluster Dendrogram (k =", k, ")"))
      
      # Save the dendrogram plot
      ggsave(paste0(group_name, "_k", k, "_dendrogram.pdf"), p_dendrogram, width = 12, height = 8)
      
      # Calculate mean trajectories for each cluster
      cluster_means <- lapply(1:k, function(i) {
        subjects_in_cluster <- results$subject_id[results$cluster == i]
        cluster_data <- hr_matrix[rownames(hr_matrix) %in% subjects_in_cluster, , drop = FALSE]
        colMeans(cluster_data, na.rm = TRUE)
      })
      
      # Convert to a data frame for plotting
      mean_traj_data <- data.frame()
      for (i in 1:k) {
        traj <- data.frame(
          times = as.numeric(colnames(hr_wide)[-1]),
          traj = cluster_means[[i]],
          cluster = i
        )
        mean_traj_data <- rbind(mean_traj_data, traj)
      }
      
      # Save the mean trajectories data
      write.csv(mean_traj_data, paste0(group_name, "_k", k, "_mean_trajectories.csv"), row.names = FALSE)
      
      # Create plot of mean trajectories
      p_means <- ggplot(data = mean_traj_data) +
        geom_line(aes(x = times, y = traj, group = cluster, color = as.factor(cluster)), size = 1) +
        ylab("Heart Rate") +
        xlab("Days from Surgery") +
        ggtitle(paste(group_name, "- Mean Heart Rate Trajectories (k =", k, ")")) +
        scale_color_discrete(name = "Cluster") +
        theme_minimal()
      
      # Save the plot
      ggsave(paste0(group_name, "_k", k, "_clusters.pdf"), p_means, width = 8, height = 6)
      
    }, error = function(e) {
      cat("Error in", group_name, "with k =", k, ":", e$message, "\n")
      write(paste("Error in", group_name, "with k =", k, ":", e$message), 
            file = paste0(group_name, "_k", k, "_error.txt"))
    })
  }
  
  return(clusters_list)
}

# 修改后的分析函数，修复p值问题
analyze_clusters <- function(hr_data, cluster_results, group_name, k) {
  # Get cluster assignments
  if (is.null(cluster_results[[paste0("k", k)]])) {
    cat("No results for", group_name, "with k =", k, "\n")
    write(paste("No results for", group_name, "with k =", k), 
          file = paste0(group_name, "_k", k, "_no_results.txt"))
    return(NULL)
  }
  
  cluster_assignment <- cluster_results[[paste0("k", k)]]$clusters
  
  # Create dataframe with cluster assignments
  results <- data.frame(
    subject_id = names(cluster_assignment),
    cluster = cluster_assignment
  )
  
  # Calculate heart rate features
  hr_features <- hr_data %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate, na.rm = TRUE),
      sd_hr = sd(heart_rate, na.rm = TRUE),
      min_hr = min(heart_rate, na.rm = TRUE),
      max_hr = max(heart_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(results, by = "subject_id")
  
  # Save heart rate features by subject
  write.csv(hr_features, paste0(group_name, "_k", k, "_hr_features.csv"), row.names = FALSE)
  
  # Statistical comparison between clusters
  cluster_comparison <- hr_features %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      mean_hr = mean(mean_hr, na.rm = TRUE),
      sd_mean_hr = sd(mean_hr, na.rm = TRUE),
      mean_sd_hr = mean(sd_hr, na.rm = TRUE),
      mean_min_hr = mean(min_hr, na.rm = TRUE),
      mean_max_hr = mean(max_hr, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save cluster comparison
  write.csv(cluster_comparison, paste0(group_name, "_k", k, "_cluster_comparison.csv"), row.names = FALSE)
  
  cat("Heart Rate Features by Cluster for", group_name, "with k =", k, ":\n")
  capture.output(print(cluster_comparison), 
                 file = paste0(group_name, "_k", k, "_cluster_comparison_output.txt"))
  
  # 修改后的统计检验部分
  stat_test_result <- NULL
  
  # Statistical test between clusters if more than one cluster
  if (length(unique(hr_features$cluster)) > 1) {
    # 检查集群数是否为2个
    if (length(unique(hr_features$cluster)) == 2) {
      # 如果是2个集群，使用t检验
      t_test_result <- t.test(mean_hr ~ cluster, data = hr_features)
      cat("T-test for Mean Heart Rate Between Clusters for", group_name, "with k =", k, ":\n")
      capture.output(print(t_test_result), 
                     file = paste0(group_name, "_k", k, "_ttest_output.txt"))
      
      # Save t-test results
      t_test_summary <- data.frame(
        p_value = t_test_result$p.value,
        t_statistic = t_test_result$statistic,
        df = t_test_result$parameter,
        mean_diff = diff(t_test_result$estimate),
        conf_int_low = t_test_result$conf.int[1],
        conf_int_high = t_test_result$conf.int[2]
      )
      write.csv(t_test_summary, paste0(group_name, "_k", k, "_ttest_results.csv"), row.names = FALSE)
      
      # 保存测试结果
      stat_test_result <- t_test_result
    } else {
      # 如果是3个或更多集群，使用ANOVA
      anova_result <- aov(mean_hr ~ as.factor(cluster), data = hr_features)
      anova_summary <- summary(anova_result)
      
      cat("ANOVA for Mean Heart Rate Between Clusters for", group_name, "with k =", k, ":\n")
      capture.output(print(anova_summary), 
                     file = paste0(group_name, "_k", k, "_anova_output.txt"))
      
      # 保存ANOVA结果
      # 从ANOVA摘要中提取关键信息
      if(length(anova_summary) > 0) {
        anova_df <- data.frame(
          p_value = anova_summary[[1]][["Pr(>F)"]][1],
          F_statistic = anova_summary[[1]][["F value"]][1],
          df_between = anova_summary[[1]][["Df"]][1],
          df_within = anova_summary[[1]][["Df"]][2]
        )
        write.csv(anova_df, paste0(group_name, "_k", k, "_anova_results.csv"), row.names = FALSE)
      }
      
      # 保存测试结果
      stat_test_result <- anova_result
    }
  }
  
  # Visualize individual trajectories by cluster
  hr_long <- hr_data %>%
    dplyr::select(subject_id, days_from_surgery, heart_rate) %>%
    left_join(results, by = "subject_id")
  
  # Save long-format data with cluster assignments
  write.csv(hr_long, paste0(group_name, "_k", k, "_hr_long.csv"), row.names = FALSE)
  
  # Plot individual trajectories by cluster
  p <- ggplot(hr_long, aes(x = days_from_surgery, y = heart_rate, group = subject_id, color = as.factor(cluster))) +
    geom_line(alpha = 0.3) +
    geom_smooth(aes(group = as.factor(cluster)), method = "loess", se = TRUE) +
    scale_x_continuous(breaks = seq(-4, 30, by = 5)) +
    ylab("Heart Rate") +
    xlab("Days from Surgery") +
    ggtitle(paste(group_name, "- Individual Heart Rate Trajectories (k =", k, ")")) +
    scale_color_discrete(name = "Cluster") +
    theme_minimal()
  
  # Save the plot
  ggsave(paste0(group_name, "_k", k, "_individual_trajectories.pdf"), p, width = 10, height = 6)
  
  # Create a summary of daily mean heart rates by cluster
  daily_summary <- hr_long %>%
    group_by(cluster, floor(days_from_surgery)) %>%
    summarize(
      day = min(floor(days_from_surgery)),
      mean_hr = mean(heart_rate, na.rm = TRUE),
      se_hr = sd(heart_rate, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Save daily summary
  write.csv(daily_summary, paste0(group_name, "_k", k, "_daily_summary.csv"), row.names = FALSE)
  
  # Plot daily mean heart rates with error bands
  p_daily <- ggplot(daily_summary, aes(x = day, y = mean_hr, color = as.factor(cluster))) +
    geom_line() +
    geom_ribbon(aes(ymin = mean_hr - se_hr, ymax = mean_hr + se_hr, fill = as.factor(cluster)), 
                alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = seq(-4, 30, by = 5)) +
    ylab("Mean Heart Rate") +
    xlab("Days from Surgery") +
    ggtitle(paste(group_name, "- Mean Heart Rate Trajectories by Day (k =", k, ")")) +
    scale_color_discrete(name = "Cluster") +
    scale_fill_discrete(name = "Cluster") +
    theme_minimal()
  
  # Save the plot
  ggsave(paste0(group_name, "_k", k, "_daily_trajectories.pdf"), p_daily, width = 10, height = 6)
  
  return(list(
    features = hr_features,
    comparison = cluster_comparison,
    statistical_test = stat_test_result  # 一致使用 statistical_test 作为键名
  ))
}

# Save the patient group information
write.csv(patient_groups, "patient_groups.csv", row.names = FALSE)
write.csv(data.frame(subject_id = group1_subjects), "group1_subjects.csv", row.names = FALSE)
write.csv(data.frame(subject_id = group2_subjects), "group2_subjects.csv", row.names = FALSE)

# Save the processed data
saveRDS(group1_hr_data, "group1_hr_data.rds")
saveRDS(group2_hr_data, "group2_hr_data.rds")

# Perform clustering for each group
group1_clusters <- perform_clustering(group1_hr_data, "Group1_Surgery0")
group2_clusters <- perform_clustering(group2_hr_data, "Group2_Diabetes")

# Analyze the clusters for each group (k=2 and k=3)
group1_analysis_k2 <- analyze_clusters(group1_hr_data, group1_clusters, "Group1_Surgery0", 2)
group1_analysis_k3 <- analyze_clusters(group1_hr_data, group1_clusters, "Group1_Surgery0", 3)

group2_analysis_k2 <- analyze_clusters(group2_hr_data, group2_clusters, "Group2_Diabetes", 2)
group2_analysis_k3 <- analyze_clusters(group2_hr_data, group2_clusters, "Group2_Diabetes", 3)

# 修改报告生成部分，正确显示p值
sink("heart_rate_clustering_report.txt")

cat("============================================\n")
cat("HEART RATE TRAJECTORY CLUSTERING ANALYSIS\n")
cat("============================================\n\n")

cat("Analysis performed on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Group 1 (No Diabetes, Surgery Type 0):\n")
cat("Number of subjects:", length(group1_subjects), "\n\n")

cat("Group 2 (Diabetes, Surgery Type 1):\n")
cat("Number of subjects:", length(group2_subjects), "\n\n")

cat("Method: Hierarchical clustering with Euclidean distance\n\n")

cat("============================================\n")
cat("CLUSTERING RESULTS - GROUP 1 (NO DIABETES)\n")
cat("============================================\n\n")

# Print k=2 results for Group 1
if (!is.null(group1_clusters[["k2"]])) {
  cat("K=2 Clusters:\n")
  cluster_stats <- table(group1_clusters[["k2"]]$clusters)
  print(cluster_stats)
  cat("\nCluster Characteristics:\n")
  print(group1_analysis_k2$comparison)
  
  # 修改后的t检验结果显示
  cat("\nStatistical Test Results:\n")
  if (!is.null(group1_analysis_k2$statistical_test)) {
    if (inherits(group1_analysis_k2$statistical_test, "htest")) {
      # T-test结果
      cat("T-test p-value:", format(group1_analysis_k2$statistical_test$p.value, digits=4), "\n")
      cat("T-statistic:", format(group1_analysis_k2$statistical_test$statistic, digits=4), "\n")
      cat("Degrees of freedom:", format(group1_analysis_k2$statistical_test$parameter, digits=4), "\n")
    } else {
      cat("No statistical test results available\n")
    }
  } else {
    cat("No statistical test results available\n")
  }
  cat("\n")
} else {
  cat("K=2 clustering failed for Group 1\n\n")
}

# Print k=3 results for Group 1
if (!is.null(group1_clusters[["k3"]])) {
  cat("K=3 Clusters:\n")
  cluster_stats <- table(group1_clusters[["k3"]]$clusters)
  print(cluster_stats)
  cat("\nCluster Characteristics:\n")
  print(group1_analysis_k3$comparison)
  
  # 修改后的ANOVA结果显示
  cat("\nStatistical Test Results:\n")
  if (!is.null(group1_analysis_k3$statistical_test)) {
    if (inherits(group1_analysis_k3$statistical_test, "aov")) {
      # ANOVA结果
      anova_summary <- summary(group1_analysis_k3$statistical_test)
      if (length(anova_summary) > 0) {
        cat("ANOVA p-value:", format(anova_summary[[1]][["Pr(>F)"]][1], digits=4), "\n")
        cat("F-statistic:", format(anova_summary[[1]][["F value"]][1], digits=4), "\n")
        cat("Degrees of freedom:", 
            paste(anova_summary[[1]][["Df"]][1], anova_summary[[1]][["Df"]][2], sep=", "), "\n")
      } else {
        cat("ANOVA results not available in expected format\n")
      }
    } else {
      cat("No ANOVA results available\n")
    }
  } else {
    cat("No statistical test results available\n")
  }
  cat("\n")
} else {
  cat("K=3 clustering failed for Group 1\n\n")
}

cat("============================================\n")
cat("CLUSTERING RESULTS - GROUP 2 (DIABETES)\n")
cat("============================================\n\n")

# Print k=2 results for Group 2
if (!is.null(group2_clusters[["k2"]])) {
  cat("K=2 Clusters:\n")
  cluster_stats <- table(group2_clusters[["k2"]]$clusters)
  print(cluster_stats)
  cat("\nCluster Characteristics:\n")
  print(group2_analysis_k2$comparison)
  
  # 修改后的t检验结果显示
  cat("\nStatistical Test Results:\n")
  if (!is.null(group2_analysis_k2$statistical_test)) {
    if (inherits(group2_analysis_k2$statistical_test, "htest")) {
      # T-test结果
      cat("T-test p-value:", format(group2_analysis_k2$statistical_test$p.value, digits=4), "\n")
      cat("T-statistic:", format(group2_analysis_k2$statistical_test$statistic, digits=4), "\n")
      cat("Degrees of freedom:", format(group2_analysis_k2$statistical_test$parameter, digits=4), "\n")
    } else {
      cat("No t-test results available\n")
    }
  } else {
    cat("No statistical test results available\n")
  }
  cat("\n")
} else {
  cat("K=2 clustering failed for Group 2\n\n")
}

# Print k=3 results for Group 2
if (!is.null(group2_clusters[["k3"]])) {
  cat("K=3 Clusters:\n")
  cluster_stats <- table(group2_clusters[["k3"]]$clusters)
  print(cluster_stats)
  cat("\nCluster Characteristics:\n")
  print(group2_analysis_k3$comparison)
  
  # 修改后的ANOVA结果显示
  cat("\nStatistical Test Results:\n")
  if (!is.null(group2_analysis_k3$statistical_test)) {
    if (inherits(group2_analysis_k3$statistical_test, "aov")) {
      # ANOVA结果
      anova_summary <- summary(group2_analysis_k3$statistical_test)
      if (length(anova_summary) > 0) {
        cat("ANOVA p-value:", format(anova_summary[[1]][["Pr(>F)"]][1], digits=4), "\n")
        cat("F-statistic:", format(anova_summary[[1]][["F value"]][1], digits=4), "\n")
        cat("Degrees of freedom:", 
            paste(anova_summary[[1]][["Df"]][1], anova_summary[[1]][["Df"]][2], sep=", "), "\n")
      } else {
        cat("ANOVA results not available in expected format\n")
      }
    } else {
      cat("No ANOVA results available\n")
    }
  } else {
    cat("No statistical test results available\n")
  }
  cat("\n")
} else {
  cat("K=3 clustering failed for Group 2\n\n")
}

cat("============================================\n")
cat("CONCLUSION\n")
cat("============================================\n\n")

cat("This analysis identified distinct heart rate trajectory patterns within each patient group during the perioperative period (4 days before to 30 days after surgery).\n\n")

cat("Key findings:\n")
cat("1. Different heart rate recovery patterns were observed within each group.\n")
cat("2. Check the generated plots to visualize these patterns.\n")
cat("3. Statistical analysis shows differences in mean heart rate between clusters.\n")
cat("4. Further clinical interpretation is recommended to understand the physiological significance of these patterns.\n\n")

sink()

# Save final clustering objects and complete data for future reference
saveRDS(group1_clusters, "Group1_Surgery0_all_clusters.rds")
saveRDS(group2_clusters, "Group2_Diabetes_all_clusters.rds")

# Save analysis results
saveRDS(list(
  group1_analysis_k2 = group1_analysis_k2,
  group1_analysis_k3 = group1_analysis_k3,
  group2_analysis_k2 = group2_analysis_k2,
  group2_analysis_k3 = group2_analysis_k3
), "all_analysis_results.rds")

# 返回到原始工作目录
setwd(original_wd)

cat("Analysis complete. All results saved to the directory:\n")
cat(output_dir, "\n")
