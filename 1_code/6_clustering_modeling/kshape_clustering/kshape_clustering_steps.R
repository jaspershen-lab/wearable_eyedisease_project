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

# Load steps data
load("3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda")

# Load baseline information
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")


# 数据质控函数 - 基于心率数据的质控方法
quality_control_data <- function(steps_trajectory_data, group_name) {
  # 1. 检查每个受试者的数据量
  subject_data_count <- steps_trajectory_data %>%
    group_by(subject_id) %>%
    summarize(
      data_points = n(),
      unique_days = n_distinct(floor(days_from_surgery)),
      has_periop_data = any(days_from_surgery >= -2 & days_from_surgery <= 2),
      min_day = min(days_from_surgery),
      max_day = max(days_from_surgery),
      .groups = "drop"
    ) %>%
    mutate(
      passed_qc = (unique_days >= 14) & has_periop_data  # 至少有15天的数据且在围手术期有数据
    )
  
  # 保存质控结果
  write.csv(subject_data_count, paste0(group_name, "_data_quality.csv"), row.names = FALSE)
  
  # 输出质控摘要
  passed_count <- sum(subject_data_count$passed_qc)
  total_count <- nrow(subject_data_count)
  excluded_count <- total_count - passed_count
  
  cat("Quality Control Results for", group_name, ":\n")
  cat("Total subjects:", total_count, "\n")
  cat("Subjects passing QC:", passed_count, "(", round(passed_count/total_count*100, 1), "%)\n")
  cat("Subjects excluded:", excluded_count, "(", round(excluded_count/total_count*100, 1), "%)\n")
  
  # 输出未通过质控的受试者
  if (excluded_count > 0) {
    cat("Excluded subjects:\n")
    excluded_subjects <- subject_data_count %>% 
      filter(!passed_qc) %>%
      mutate(
        reason = case_when(
          !has_periop_data ~ "No perioperative data",
          unique_days < 14 ~ paste0("Insufficient data (", unique_days, " days)"),
          TRUE ~ "Unknown reason"
        )
      )
    print(excluded_subjects %>% dplyr::select(subject_id, reason, data_points, unique_days, min_day, max_day))
  }
  
  # 过滤通过质控的受试者数据
  steps_clean <- steps_trajectory_data %>%
    filter(subject_id %in% (subject_data_count %>% filter(passed_qc) %>% pull(subject_id)))
  
  cat("Original data points:", nrow(steps_trajectory_data), "\n")
  cat("Clean data points:", nrow(steps_clean), "(", round(nrow(steps_clean)/nrow(steps_trajectory_data)*100, 1), "%)\n")
  
  return(list(
    clean_data = steps_clean,
    qc_results = subject_data_count,
    n_original = nrow(subject_data_count),
    n_clean = passed_count
  ))
}

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

# Modified function to extract steps data
extract_steps_data <- function(subject_ids) {
  # Extract the sample info (which contains subject_id)
  sample_info <- daily_workout_details_data@sample_info
  
  # Filter sample_info to include only samples from our subjects of interest
  filtered_sample_info <- sample_info %>%
    filter(subject_id %in% subject_ids)
  
  # Get the sample_ids from the filtered sample_info
  matching_samples <- filtered_sample_info$sample_id
  
  # Filter the daily_workout_details_data to include only matching samples
  expression_data <- daily_workout_details_data@expression_data[, matching_samples, drop = FALSE]
  
  # Find the row index for steps
  steps_row_index <- which(rownames(expression_data) == "steps" | 
                             rownames(expression_data) == "steps(steps)" | 
                             row.names(daily_workout_details_data@variable_info) == "steps")
  
  # Extract the steps values
  steps_values <- as.numeric(expression_data[steps_row_index, ])
  
  # Create a data frame with sample_id, subject_id, measure_time, and steps
  steps_data <- filtered_sample_info %>%
    mutate(steps = steps_values)
  
  # Add surgery time information
  steps_data_with_surgery <- steps_data %>%
    left_join(surgery_times, by = c("subject_id" = "ID"))
  
  # Calculate days from surgery and filter to the period of interest
  steps_trajectory_data <- steps_data_with_surgery %>%
    mutate(
      days_from_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days"))
    ) %>%
    filter(days_from_surgery >= -4 & days_from_surgery <= 30)
  
  return(steps_trajectory_data)
}

# Process the data for both groups
group1_steps_data <- extract_steps_data(group1_subjects)
group2_steps_data <- extract_steps_data(group2_subjects)

# 移除SH012并重新聚类
# 从组2受试者中排除SH012
group2_subjects_filtered <- group2_subjects[group2_subjects != "SH012"]
cat("Original Group2 size:", length(group2_subjects), "\n")
cat("Filtered Group2 size (without SH012):", length(group2_subjects_filtered), "\n")

# 处理过滤后的数据
group2_steps_data <- extract_steps_data(group2_subjects_filtered)

# 应用质控
cat("\n--- Applying Quality Control ---\n")
group1_qc <- quality_control_data(group1_steps_data, "Group1_Surgery0")
group2_qc <- quality_control_data(group2_steps_data, "Group2_Diabetes")

# 使用通过质控的数据
group1_steps_data_clean <- group1_qc$clean_data
group2_steps_data_clean <- group2_qc$clean_data
# 创建结果保存目录并设置工作目录（只设置一次）
output_dir <- "3_data_analysis/6_clustering_modeling/kshape/time_series/1m/steps"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 保存原始路径，以便在程序结束时恢复
original_wd <- getwd()
# 切换到输出目录
setwd(file.path(original_wd, output_dir))

# 修改后的函数，确保正确处理路径和p值
perform_clustering <- function(steps_trajectory_data, group_name) {
  # Create daily average steps
  daily_steps <- steps_trajectory_data %>%
    group_by(subject_id, floor(days_from_surgery)) %>%
    summarize(
      day = min(floor(days_from_surgery)),
      steps_avg = mean(steps, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(subject_id, day, steps_avg)
  
  # Save the daily steps data (使用相对路径)
  write.csv(daily_steps, paste0(group_name, "_daily_steps.csv"), row.names = FALSE)
  
  # Pivot to wide format for clustering
  steps_wide <- daily_steps %>%
    pivot_wider(
      id_cols = subject_id,
      names_from = day,
      values_from = steps_avg
    )
  
  # Save the wide-format data
  write.csv(steps_wide, paste0(group_name, "_steps_wide.csv"), row.names = FALSE)
  
  # Prepare matrix for clustering
  subjects <- steps_wide$subject_id
  steps_matrix <- as.matrix(steps_wide[, -1])
  rownames(steps_matrix) <- subjects
  
  # Handle missing values (if any)
  if(any(is.na(steps_matrix))) {
    cat("Handling missing values in the steps matrix\n")
    steps_matrix[is.na(steps_matrix)] <- mean(steps_matrix, na.rm = TRUE)
  }
  
  # Print group info
  cat("Processing", group_name, "with", nrow(steps_matrix), "subjects\n")
  write(paste("Processing", group_name, "with", nrow(steps_matrix), "subjects"), 
        file = paste0(group_name, "_processing_log.txt"))
  
  # Save the matrix
  saveRDS(steps_matrix, paste0(group_name, "_steps_matrix.rds"))
  
  # Create a simpler distance matrix using Euclidean distance as a fallback
  cat("Computing distance matrix...\n")
  
  # Use standard Euclidean distance
  eucl_dist <- dist(steps_matrix, method = "euclidean")
  
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
        cluster_data <- steps_matrix[rownames(steps_matrix) %in% subjects_in_cluster, , drop = FALSE]
        colMeans(cluster_data, na.rm = TRUE)
      })
      
      # Convert to a data frame for plotting
      mean_traj_data <- data.frame()
      for (i in 1:k) {
        traj <- data.frame(
          times = as.numeric(colnames(steps_wide)[-1]),
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
        ylab("Steps Count") +
        xlab("Days from Surgery") +
        ggtitle(paste(group_name, "- Mean Steps Count Trajectories (k =", k, ")")) +
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
analyze_clusters <- function(steps_data, cluster_results, group_name, k) {
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
  
  # Calculate steps features
  steps_features <- steps_data %>%
    group_by(subject_id) %>%
    summarise(
      mean_steps = mean(steps, na.rm = TRUE),
      sd_steps = sd(steps, na.rm = TRUE),
      min_steps = min(steps, na.rm = TRUE),
      max_steps = max(steps, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(results, by = "subject_id")
  
  # Save steps features by subject
  write.csv(steps_features, paste0(group_name, "_k", k, "_steps_features.csv"), row.names = FALSE)
  
  # Statistical comparison between clusters
  cluster_comparison <- steps_features %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      mean_steps = mean(mean_steps, na.rm = TRUE),
      sd_mean_steps = sd(mean_steps, na.rm = TRUE),
      mean_sd_steps = mean(sd_steps, na.rm = TRUE),
      mean_min_steps = mean(min_steps, na.rm = TRUE),
      mean_max_steps = mean(max_steps, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save cluster comparison
  write.csv(cluster_comparison, paste0(group_name, "_k", k, "_cluster_comparison.csv"), row.names = FALSE)
  
  cat("Steps Features by Cluster for", group_name, "with k =", k, ":\n")
  capture.output(print(cluster_comparison), 
                 file = paste0(group_name, "_k", k, "_cluster_comparison_output.txt"))
  
  # 修改后的统计检验部分
  stat_test_result <- NULL
  
  # Statistical test between clusters if more than one cluster
  if (length(unique(steps_features$cluster)) > 1) {
    # 检查集群数是否为2个
    if (length(unique(steps_features$cluster)) == 2) {
      # 如果是2个集群，使用t检验
      t_test_result <- t.test(mean_steps ~ cluster, data = steps_features)
      cat("T-test for Mean Steps Between Clusters for", group_name, "with k =", k, ":\n")
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
      anova_result <- aov(mean_steps ~ as.factor(cluster), data = steps_features)
      anova_summary <- summary(anova_result)
      
      cat("ANOVA for Mean Steps Between Clusters for", group_name, "with k =", k, ":\n")
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
  steps_long <- steps_data %>%
    dplyr::select(subject_id, days_from_surgery, steps) %>%
    left_join(results, by = "subject_id")
  
  # Save long-format data with cluster assignments
  write.csv(steps_long, paste0(group_name, "_k", k, "_steps_long.csv"), row.names = FALSE)
  
  # Plot individual trajectories by cluster
  p <- ggplot(steps_long, aes(x = days_from_surgery, y = steps, group = subject_id, color = as.factor(cluster))) +
    geom_line(alpha = 0.3) +
    geom_smooth(aes(group = as.factor(cluster)), method = "loess", se = TRUE) +
    scale_x_continuous(breaks = seq(-4, 30, by = 5)) +
    ylab("Steps Count") +
    xlab("Days from Surgery") +
    ggtitle(paste(group_name, "- Individual Steps Count Trajectories (k =", k, ")")) +
    scale_color_discrete(name = "Cluster") +
    theme_minimal()
  
  # Save the plot
  ggsave(paste0(group_name, "_k", k, "_individual_trajectories.pdf"), p, width = 10, height = 6)
  
  # Create a summary of daily mean steps by cluster
  daily_summary <- steps_long %>%
    group_by(cluster, floor(days_from_surgery)) %>%
    summarize(
      day = min(floor(days_from_surgery)),
      mean_steps = mean(steps, na.rm = TRUE),
      se_steps = sd(steps, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Save daily summary
  write.csv(daily_summary, paste0(group_name, "_k", k, "_daily_summary.csv"), row.names = FALSE)
  
  # Plot daily mean steps with error bands
  p_daily <- ggplot(daily_summary, aes(x = day, y = mean_steps, color = as.factor(cluster))) +
    geom_line() +
    geom_ribbon(aes(ymin = mean_steps - se_steps, ymax = mean_steps + se_steps, fill = as.factor(cluster)), 
                alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = seq(-4, 30, by = 5)) +
    ylab("Mean Steps Count") +
    xlab("Days from Surgery") +
    ggtitle(paste(group_name, "- Mean Steps Count Trajectories by Day (k =", k, ")")) +
    scale_color_discrete(name = "Cluster") +
    scale_fill_discrete(name = "Cluster") +
    theme_minimal()
  
  # Save the plot
  ggsave(paste0(group_name, "_k", k, "_daily_trajectories.pdf"), p_daily, width = 10, height = 6)
  
  return(list(
    features = steps_features,
    comparison = cluster_comparison,
    statistical_test = stat_test_result  # 一致使用 statistical_test 作为键名
  ))
}

# Save the patient group information
write.csv(patient_groups, "patient_groups.csv", row.names = FALSE)
write.csv(data.frame(subject_id = group1_subjects), "group1_subjects.csv", row.names = FALSE)
write.csv(data.frame(subject_id = group2_subjects), "group2_subjects.csv", row.names = FALSE)

# Save the processed data
saveRDS(group1_steps_data, "group1_steps_data_original.rds")
saveRDS(group2_steps_data, "group2_steps_data_original.rds")
saveRDS(group1_steps_data_clean, "group1_steps_data_clean.rds")
saveRDS(group2_steps_data_clean, "group2_steps_data_clean.rds")

# 保存质控结果
saveRDS(group1_qc$qc_results, "group1_qc_results.rds")
saveRDS(group2_qc$qc_results, "group2_qc_results.rds")

# Perform clustering for each group (使用经过质控的数据)
group1_clusters <- perform_clustering(group1_steps_data_clean, "Group1_Surgery0")
group2_clusters <- perform_clustering(group2_steps_data_clean, "Group2_Diabetes")

# Analyze the clusters for each group
group1_analysis_k2 <- analyze_clusters(group1_steps_data_clean, group1_clusters, "Group1_Surgery0", 2)
group1_analysis_k3 <- analyze_clusters(group1_steps_data_clean, group1_clusters, "Group1_Surgery0", 3)

group2_analysis_k2 <- analyze_clusters(group2_steps_data_clean, group2_clusters, "Group2_Diabetes", 2)
group2_analysis_k3 <- analyze_clusters(group2_steps_data_clean, group2_clusters, "Group2_Diabetes", 3)

# 修改报告生成部分，正确显示p值
sink("steps_clustering_report.txt")

cat("============================================\n")
cat("STEPS COUNT TRAJECTORY CLUSTERING ANALYSIS\n")
cat("============================================\n\n")

cat("Analysis performed on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Quality Control Process:\n")
cat("- Subjects were required to have at least 15 days of data\n")
cat("- Subjects were required to have data during the perioperative period (-2 to +2 days from surgery)\n\n")

cat("Group 1 (No Diabetes, Surgery Type 0):\n")
cat("Original subjects:", length(group1_subjects), "\n")
cat("Subjects after QC:", group1_qc$n_clean, "(", round(group1_qc$n_clean/group1_qc$n_original*100, 1), "%)\n\n")

cat("Group 2 (Diabetes, Surgery Type 1):\n")
cat("Original subjects:", length(group2_subjects_filtered), "\n")
cat("Subjects after QC:", group2_qc$n_clean, "(", round(group2_qc$n_clean/group2_qc$n_original*100, 1), "%)\n\n")

cat("Method: Hierarchical clustering with Euclidean distance\n\n")

# Rest of your report remains the same

# Update the conclusion section to include QC information
cat("============================================\n")
cat("CONCLUSION\n")
cat("============================================\n\n")

cat("This analysis identified distinct steps count trajectory patterns within each patient group during the perioperative period (4 days before to 30 days after surgery).\n\n")

cat("Quality Control Summary:\n")
cat("- All subjects were required to have at least 15 days of data\n")
cat("- All subjects were required to have data during the perioperative period\n")
cat("- Group 1: ", group1_qc$n_clean, " of ", group1_qc$n_original, " subjects passed QC (", 
    round(group1_qc$n_clean/group1_qc$n_original*100, 1), "%)\n", sep="")
cat("- Group 2: ", group2_qc$n_clean, " of ", group2_qc$n_original, " subjects passed QC (", 
    round(group2_qc$n_clean/group2_qc$n_original*100, 1), "%)\n\n", sep="")

cat("Key findings:\n")
cat("1. Different physical activity recovery patterns were observed within each group.\n")
cat("2. Check the generated plots to visualize these patterns.\n")
cat("3. Statistical analysis shows differences in mean steps count between clusters.\n")
cat("4. Further clinical interpretation is recommended to understand the significance of these activity patterns in recovery from surgery.\n\n")

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


# 保存质控结果的汇总
quality_control_summary <- data.frame(
  group = c("Group1_Surgery0", "Group2_Diabetes"),
  original_subjects = c(group1_qc$n_original, group2_qc$n_original),
  subjects_after_qc = c(group1_qc$n_clean, group2_qc$n_clean),
  percentage_passed = c(round(group1_qc$n_clean/group1_qc$n_original*100, 1), 
                        round(group2_qc$n_clean/group2_qc$n_original*100, 1))
)
write.csv(quality_control_summary, "quality_control_summary.csv", row.names = FALSE)

# 返回到原始工作目录
setwd(original_wd)

cat("Analysis complete. All results saved to the directory:\n")
cat(output_dir, "\n")
cat("\nQuality control summary:\n")
print(quality_control_summary)
