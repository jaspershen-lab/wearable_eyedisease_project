# =====================================================================
# Integrated Complete Single Metric Correlation Analysis
# Includes data loading, analysis, and complete visualizations
# Output will be saved to specified directory
# =====================================================================

library(tidyverse)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(RColorBrewer)

# Store original working directory and target output directory
original_wd <- getwd()
target_output_dir <- "3_data_analysis/6_clustering_modeling/systematic_metric_selection"

cat("🎯 Integrated Complete Single Metric Correlation Analysis\n")
cat("=========================================================\n")
cat("📂 Current working directory:", original_wd, "\n")
cat("📁 Target output directory:", target_output_dir, "\n\n")

# Create target directory but don't change working directory yet
dir.create(target_output_dir, recursive = TRUE, showWarnings = FALSE)

# Clear environment except for our directory variables
all_objects <- ls()
objects_to_keep <- c("original_wd", "target_output_dir")
objects_to_remove <- setdiff(all_objects, objects_to_keep)
rm(list = objects_to_remove)

# -------------------- 1. Define Candidate Metrics --------------------

candidate_metrics <- list(
  rhr_metrics = c("cv_rhr_1", "mean_rhr_1", "max_rhr_1", "min_rhr_1"),
  steps_metrics = c("steps_max", "steps_mean", "steps_total")
)

time_windows <- c("baseline", "acute_recovery", "early_recovery", "mid_recovery", "late_recovery")

# -------------------- 2. Load and Process Wearable Data --------------------

load_wearable_data_fixed <- function() {
  
  cat("🔄 Loading and processing wearable data...\n")
  
  # Try multiple possible paths for the raw wearable data (from original working directory)
  possible_paths <- c(
    "3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv",
    "mfuzz_D_Surg1_8h_filtered.csv",
    "../mfuzz_D_Surg1_8h_filtered.csv",
    "data/mfuzz_D_Surg1_8h_filtered.csv"
  )
  
  wearable_raw_data <- NULL
  
  for(path in possible_paths) {
    if(file.exists(path)) {
      cat("✓ Found wearable data at:", path, "\n")
      wearable_raw_data <- read.csv(path, check.names = FALSE)
      break
    }
  }
  
  if(is.null(wearable_raw_data)) {
    cat("❌ Could not find raw wearable data file\n")
    cat("📝 Expected file patterns: mfuzz_D_Surg1_8h_filtered.csv\n")
    cat("🔍 Available CSV files:\n")
    csv_files <- list.files(".", pattern = "\\.csv$", recursive = TRUE)
    for(i in 1:min(10, length(csv_files))) {
      cat("  ", csv_files[i], "\n")
    }
    stop("Cannot proceed without wearable data")
  }
  
  cat("✓ Successfully loaded raw wearable data:", nrow(wearable_raw_data), "patients\n")
  
  # Process time windows
  time_windows_def <- list(
    baseline = list(days = -4:-1, name = "baseline"),
    acute_recovery = list(days = 0:3, name = "acute_recovery"),
    early_recovery = list(days = 4:7, name = "early_recovery"),
    mid_recovery = list(days = 8:15, name = "mid_recovery"),
    late_recovery = list(days = 16:30, name = "late_recovery")
  )
  
  # Initialize processed data
  processed_wearable_data <- data.frame(
    ID = wearable_raw_data$subject_id,
    stringsAsFactors = FALSE
  )
  
  # Process all metrics for all time windows
  all_metrics <- c("cv_rhr_1", "mean_rhr_1", "max_rhr_1", "min_rhr_1", 
                   "steps_max", "steps_mean", "steps_total")
  
  for(window_name in names(time_windows_def)) {
    window_info <- time_windows_def[[window_name]]
    window_days <- window_info$days
    
    cat(sprintf("  Processing %s time window (days %s)...\n", 
                window_name, paste(range(window_days), collapse = " to ")))
    
    for(metric in all_metrics) {
      # Find columns for this metric and time window
      window_cols <- c()
      for(day in window_days) {
        day_str <- paste0("day_", day, "_", metric)
        if(day_str %in% colnames(wearable_raw_data)) {
          window_cols <- c(window_cols, day_str)
        }
      }
      
      if(length(window_cols) > 0) {
        # Calculate mean for this time window
        metric_data <- wearable_raw_data %>%
          dplyr::select(subject_id, all_of(window_cols)) %>%
          mutate(
            valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
            metric_mean = ifelse(
              valid_count >= max(1, floor(length(window_cols)/2)),
              rowMeans(dplyr::select(., -subject_id, -valid_count), na.rm = TRUE),
              NA
            )
          ) %>%
          dplyr::select(subject_id, metric_mean)
        
        # Add to processed data
        col_name <- paste0(metric, "_", window_name)
        names(metric_data) <- c("ID", col_name)
        
        processed_wearable_data <- processed_wearable_data %>%
          left_join(metric_data, by = "ID")
        
        valid_count <- sum(!is.na(metric_data[[col_name]]))
        cat(sprintf("    ✓ %s: %d valid values\n", col_name, valid_count))
      }
    }
  }
  
  # Save processed data to target directory
  processed_data_path <- file.path(target_output_dir, "processed_wearable_correlation_data.csv")
  write.csv(processed_wearable_data, processed_data_path, row.names = FALSE)
  cat("✓ Saved processed data:", processed_data_path, "\n")
  
  cat("\n📊 Processing summary:\n")
  cat(sprintf("- Patients: %d\n", nrow(processed_wearable_data)))
  cat(sprintf("- Metrics: %d\n", ncol(processed_wearable_data) - 1))
  
  return(processed_wearable_data)
}

# -------------------- 3. Correlation Analysis Function --------------------

analyze_single_metrics_correlation_fixed <- function(wearable_data = NULL) {
  
  cat("========================================\n")
  cat("🎯 Single Metric vs OCTA Clustering Correlation Analysis\n")
  cat("========================================\n")
  
  # Load OCTA data
  cat("📂 Loading OCTA clustering data...\n")
  
  # Try multiple possible paths for OCTA data (from original working directory)
  octa_possible_paths <- c(
    "3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv",
    "../mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv",
    "ppv_comprehensive_cluster_results_with_outcomes.csv",
    "results/ppv_comprehensive_cluster_results_with_outcomes.csv"
  )
  
  octa_clusters <- NULL
  
  for(path in octa_possible_paths) {
    if(file.exists(path)) {
      cat("✓ Found OCTA data at:", path, "\n")
      octa_clusters <- read.csv(path)
      break
    }
  }
  
  if(is.null(octa_clusters)) {
    cat("❌ Could not find OCTA clustering data file\n")
    cat("📝 Expected file patterns: *cluster*, *ppv*, *octa*\n")
    cat("🔍 Available CSV files:\n")
    csv_files <- list.files(".", pattern = "\\.csv$", recursive = TRUE)
    for(i in 1:min(10, length(csv_files))) {
      cat("  ", csv_files[i], "\n")
    }
    stop("Cannot proceed without OCTA data")
  }
  
  # Use provided wearable data or load it
  if(is.null(wearable_data)) {
    cat("📂 Loading wearable data...\n")
    processed_data_path <- file.path(target_output_dir, "processed_wearable_correlation_data.csv")
    if(file.exists(processed_data_path)) {
      wearable_data <- read.csv(processed_data_path)
      cat("✓ Loaded existing processed wearable data from target directory\n")
    } else {
      cat("📊 Generating wearable data...\n")
      wearable_data <- load_wearable_data_fixed()
    }
  } else {
    cat("✓ Using provided wearable data\n")
  }
  
  # Merge data
  merged_data <- merge(octa_clusters, wearable_data, 
                       by.x = "subject_id", by.y = "ID", all.inner = TRUE)
  
  cat("✓ Successfully matched", nrow(merged_data), "patients\n\n")
  
  # -------------------- Analyze Each Candidate Metric --------------------
  
  all_results <- data.frame()
  
  cat("🔍 Analyzing candidate metrics vs OCTA clustering correlation:\n")
  cat("==========================================\n")
  
  # Analyze heart rate and steps metrics separately
  for(metric_type in names(candidate_metrics)) {
    metrics <- candidate_metrics[[metric_type]]
    
    cat(sprintf("\n📊 %s Analysis:\n", 
                ifelse(metric_type == "rhr_metrics", "Heart Rate Variability", "Step Count Activity")))
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    for(metric in metrics) {
      cat(sprintf("\n🔸 %s:\n", toupper(metric)))
      
      metric_results <- data.frame()
      
      for(window in time_windows) {
        col_name <- paste0(metric, "_", window)
        
        if(col_name %in% names(merged_data)) {
          
          # Get complete data
          complete_data <- merged_data[complete.cases(merged_data[, c(col_name, "max_cluster")]), ]
          
          if(nrow(complete_data) >= 5) {
            
            # Calculate correlation
            correlation <- cor(complete_data[[col_name]], complete_data$max_cluster, 
                               method = "spearman")
            
            # t-test
            t_test <- t.test(complete_data[[col_name]] ~ complete_data$max_cluster)
            
            # Descriptive statistics
            cluster_stats <- complete_data %>%
              group_by(max_cluster) %>%
              summarise(
                n = n(),
                mean = mean(.data[[col_name]], na.rm = TRUE),
                sd = sd(.data[[col_name]], na.rm = TRUE),
                .groups = 'drop'
              )
            
            # Effect size
            if(nrow(cluster_stats) == 2) {
              pooled_sd <- sqrt(((cluster_stats$n[1]-1)*cluster_stats$sd[1]^2 + 
                                   (cluster_stats$n[2]-1)*cluster_stats$sd[2]^2) / 
                                  (cluster_stats$n[1] + cluster_stats$n[2] - 2))
              cohens_d <- abs(cluster_stats$mean[2] - cluster_stats$mean[1]) / pooled_sd
            } else {
              cohens_d <- NA
            }
            
            # Save results
            metric_results <- rbind(metric_results, data.frame(
              metric_type = metric_type,
              metric = metric,
              time_window = window,
              n_patients = nrow(complete_data),
              correlation = correlation,
              p_value = t_test$p.value,
              cohens_d = cohens_d,
              cluster1_mean = cluster_stats$mean[cluster_stats$max_cluster == 1],
              cluster2_mean = cluster_stats$mean[cluster_stats$max_cluster == 2],
              cluster1_n = cluster_stats$n[cluster_stats$max_cluster == 1],
              cluster2_n = cluster_stats$n[cluster_stats$max_cluster == 2],
              is_significant = t_test$p.value < 0.05,
              correlation_strength = case_when(
                abs(correlation) >= 0.7 ~ "Strong",
                abs(correlation) >= 0.5 ~ "Moderate",
                abs(correlation) >= 0.3 ~ "Weak", 
                abs(correlation) >= 0.1 ~ "Very Weak",
                TRUE ~ "Negligible"
              ),
              direction = ifelse(correlation > 0, "Positive", "Negative"),
              stringsAsFactors = FALSE
            ))
            
            # Print results
            cat(sprintf("  %s: r=%.3f, p=%.3f, d=%.3f (%s%s)\n",
                        str_pad(window, 15), correlation, t_test$p.value, cohens_d,
                        ifelse(t_test$p.value < 0.05, "Significant", "Non-significant"),
                        ifelse(abs(correlation) >= 0.3, " ⭐", "")))
          } else {
            cat(sprintf("  %s: Insufficient data (n=%d)\n", str_pad(window, 15), nrow(complete_data)))
          }
        } else {
          cat(sprintf("  %s: Column not found\n", str_pad(window, 15)))
        }
      }
      
      # Find best time window for each metric
      if(nrow(metric_results) > 0) {
        best_performance <- metric_results[which.max(abs(metric_results$correlation)), ]
        cat(sprintf("  → Best performance: %s (r=%.3f)\n", 
                    best_performance$time_window, best_performance$correlation))
        
        all_results <- rbind(all_results, metric_results)
      } else {
        cat("  → No valid results for this metric\n")
      }
    }
  }
  
  return(list(
    all_results = all_results,
    merged_data = merged_data
  ))
}

# -------------------- 4. Results Summary Function --------------------

summarize_single_metrics_results <- function(all_results) {
  
  cat("\n📈 Single Metric Correlation Analysis Summary:\n")
  cat("==========================\n")
  
  # Overall statistics
  total_tests <- nrow(all_results)
  significant_tests <- sum(all_results$is_significant)
  meaningful_correlations <- sum(abs(all_results$correlation) >= 0.3)
  weak_correlations <- sum(abs(all_results$correlation) >= 0.1 & abs(all_results$correlation) < 0.3)
  
  cat("📊 Overall Statistics:\n")
  cat(sprintf("- Total tests: %d (metrics × time windows)\n", total_tests))
  cat(sprintf("- Significant correlations: %d (%.1f%%)\n", significant_tests, 100*significant_tests/total_tests))
  cat(sprintf("- Meaningful correlations (|r|≥0.3): %d (%.1f%%)\n", meaningful_correlations, 100*meaningful_correlations/total_tests))
  cat(sprintf("- Weak correlation trends (|r|≥0.1): %d (%.1f%%)\n", weak_correlations, 100*weak_correlations/total_tests))
  
  # Top performing metrics
  cat("\n🏆 Top 5 Best Performing Single Metrics:\n")
  
  top_metrics <- all_results %>%
    arrange(desc(abs(correlation))) %>%
    head(5)
  
  for(i in 1:nrow(top_metrics)) {
    row <- top_metrics[i, ]
    cat(sprintf("%2d. %s (%s): r=%.3f, p=%.3f %s\n",
                i, row$metric, row$time_window, row$correlation, row$p_value,
                ifelse(row$is_significant, "⭐", "")))
  }
  
  return(list(
    total_tests = total_tests,
    significant_tests = significant_tests,
    meaningful_correlations = meaningful_correlations,
    top_metrics = top_metrics
  ))
}

# -------------------- 5. Complete Visualization Function --------------------

create_complete_single_metrics_visualizations <- function(all_results, merged_data) {
  
  cat("\n🎨 Creating Complete Single Metric Correlation Visualizations...\n")
  cat("================================================================\n")
  
  # Create output directory in target location
  viz_output_dir <- file.path(target_output_dir, "single_metrics_correlation")
  dir.create(viz_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Check if we have results
  if(nrow(all_results) == 0) {
    cat("❌ No results to visualize\n")
    return(NULL)
  }
  
  cat("📊 Generating visualizations in:", viz_output_dir, "\n")
  
  # -------------------- 1. Correlation Heatmap --------------------
  cat("1. Creating correlation heatmap...\n")
  
  # Prepare correlation matrix
  correlation_matrix <- all_results %>%
    dplyr::select(metric, time_window, correlation) %>%
    pivot_wider(names_from = time_window, values_from = correlation, values_fill = 0) %>%
    column_to_rownames("metric") %>%
    as.matrix()
  
  # Create correlation heatmap (PDF)
  pdf(file.path(viz_output_dir, "correlation_heatmap.pdf"), width = 12, height = 8)
  corrplot(correlation_matrix,
           method = "color",
           type = "full",
           tl.cex = 0.8,
           tl.col = "black",
           title = "Single Metrics vs OCTA Clustering Correlation Heatmap",
           mar = c(0,0,3,0),
           col = colorRampPalette(c("#d73027", "white", "#1a9850"))(100),
           addCoef.col = "black",
           number.cex = 0.7)
  dev.off()
  
  # Create correlation heatmap (PNG)
  png(file.path(viz_output_dir, "correlation_heatmap.png"), width = 1200, height = 800, res = 100)
  corrplot(correlation_matrix,
           method = "color", 
           type = "full",
           tl.cex = 0.8,
           tl.col = "black",
           title = "Single Metrics vs OCTA Clustering Correlation Heatmap",
           mar = c(0,0,3,0),
           col = colorRampPalette(c("#d73027", "white", "#1a9850"))(100),
           addCoef.col = "black",
           number.cex = 0.7)
  dev.off()
  
  cat("   ✓ Saved: correlation_heatmap.pdf/png\n")
  
  # -------------------- 2. Top Metrics Bar Plot --------------------
  cat("2. Creating top metrics bar plot...\n")
  
  top_10 <- all_results %>%
    arrange(desc(abs(correlation))) %>%
    head(10) %>%
    mutate(
      metric_label = paste0(metric, "\n(", time_window, ")"),
      metric_type_name = ifelse(metric_type == "rhr_metrics", "Heart Rate Metrics", "Step Count Metrics")
    )
  
  p_top_metrics <- ggplot(top_10, aes(x = reorder(metric_label, abs(correlation)))) +
    geom_col(aes(y = correlation, fill = metric_type_name), alpha = 0.8, width = 0.7) +
    geom_text(aes(y = correlation, 
                  label = paste0("r=", sprintf("%.3f", correlation),
                                 ifelse(is_significant, " *", ""))),
              hjust = ifelse(top_10$correlation >= 0, -0.1, 1.1), 
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("Heart Rate Metrics" = "#3498db", "Step Count Metrics" = "#e74c3c"),
                      name = "Metric Type") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
    labs(
      title = "Top 10 Single Metrics vs OCTA Clustering Correlation",
      subtitle = "* indicates statistical significance (p<0.05)",
      x = "Metric (Time Window)",
      y = "Spearman Correlation Coefficient",
      caption = "Positive values = positive correlation, Negative values = negative correlation"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(viz_output_dir, "top_metrics_barplot.pdf"), 
         p_top_metrics, width = 12, height = 10)
  ggsave(file.path(viz_output_dir, "top_metrics_barplot.png"), 
         p_top_metrics, width = 12, height = 10, dpi = 300)
  
  cat("   ✓ Saved: top_metrics_barplot.pdf/png\n")
  
  # -------------------- 3. Correlation by Type Comparison --------------------
  cat("3. Creating correlation by type comparison...\n")
  
  p_by_type <- ggplot(all_results, aes(x = metric, y = correlation, fill = time_window)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.8) +
    facet_wrap(~ metric_type, scales = "free_x", 
               labeller = labeller(metric_type = c("rhr_metrics" = "Heart Rate Variability Metrics",
                                                   "steps_metrics" = "Step Count Activity Metrics"))) +
    scale_fill_brewer(palette = "Set3", name = "Time Window") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = c(-0.3, 0.3), linetype = "dotted", color = "red", alpha = 0.7) +
    labs(
      title = "All Candidate Metrics vs OCTA Clustering Correlation Comparison",
      subtitle = "Grouped by metric type and time window | Red dotted lines = ±0.3 (meaningful correlation threshold)",
      x = "Metric",
      y = "Spearman Correlation Coefficient"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(viz_output_dir, "correlation_by_type.pdf"), 
         p_by_type, width = 14, height = 10)
  ggsave(file.path(viz_output_dir, "correlation_by_type.png"), 
         p_by_type, width = 14, height = 10, dpi = 300)
  
  cat("   ✓ Saved: correlation_by_type.pdf/png\n")
  
  # -------------------- 4. Correlation Strength Distribution --------------------
  cat("4. Creating correlation strength distribution...\n")
  
  # Create ordered factor for strength categories
  all_results$correlation_strength_ordered <- factor(all_results$correlation_strength,
                                                     levels = c("Negligible", "Very Weak", "Weak", "Moderate", "Strong"),
                                                     ordered = TRUE)
  
  p_strength_dist <- ggplot(all_results, aes(x = correlation_strength_ordered, fill = metric_type)) +
    geom_bar(position = "dodge", alpha = 0.8, width = 0.7) +
    geom_text(stat = "count", aes(label = after_stat(count)), 
              position = position_dodge(width = 0.7), vjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("rhr_metrics" = "#3498db", "steps_metrics" = "#e74c3c"),
                      labels = c("rhr_metrics" = "Heart Rate Metrics", "steps_metrics" = "Step Count Metrics"),
                      name = "Metric Type") +
    labs(
      title = "Single Metric Correlation Strength Distribution",
      subtitle = "Number of metrics in each strength category",
      x = "Correlation Strength Category",
      y = "Number of Metrics"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(viz_output_dir, "correlation_strength_distribution.pdf"), 
         p_strength_dist, width = 10, height = 8)
  ggsave(file.path(viz_output_dir, "correlation_strength_distribution.png"), 
         p_strength_dist, width = 10, height = 8, dpi = 300)
  
  cat("   ✓ Saved: correlation_strength_distribution.pdf/png\n")
  
  cat("\n🎉 Core Visualizations Generated Successfully!\n")
  cat("===============================================\n")
  cat("📂 Location:", viz_output_dir, "\n")
  cat("📊 Files created: 8 (4 visualizations × 2 formats)\n")
  
  return(list(
    correlation_matrix = correlation_matrix,
    top_metrics = top_10
  ))
}

# -------------------- 6. Generate Report --------------------

generate_single_metrics_report <- function(all_results, summary_results) {
  
  # Key statistics
  total_tests <- nrow(all_results)
  significant_tests <- sum(all_results$is_significant)
  meaningful_correlations <- sum(abs(all_results$correlation) >= 0.3)
  weak_trends <- sum(abs(all_results$correlation) >= 0.1 & abs(all_results$correlation) < 0.3)
  
  top_metric <- all_results[which.max(abs(all_results$correlation)), ]
  
  report <- paste0(
    "========================================\n",
    "Single Metrics vs OCTA Clustering Correlation Analysis Report\n",
    "========================================\n\n",
    
    "🎯 Analysis Objective:\n",
    "Directly evaluate the correlation between each candidate wearable metric and OCTA clustering\n",
    "No composite metrics used, focusing on individual metric predictive potential\n\n",
    
    "📊 Candidate Metrics:\n",
    "- Heart Rate Variability: cv_rhr_1, mean_rhr_1, max_rhr_1, min_rhr_1\n",
    "- Step Count Activity: steps_max, steps_mean, steps_total\n",
    "- Time Windows: baseline, acute_recovery, early_recovery, mid_recovery, late_recovery\n\n",
    
    "📈 Overall Results:\n",
    sprintf("- Total tests: %d (7 metrics × 5 time windows)\n", total_tests),
    sprintf("- Significant correlations: %d (%.1f%%)\n", significant_tests, 100*significant_tests/total_tests),
    sprintf("- Meaningful correlations (|r|≥0.3): %d (%.1f%%)\n", meaningful_correlations, 100*meaningful_correlations/total_tests),
    sprintf("- Weak correlation trends (0.1≤|r|<0.3): %d (%.1f%%)\n", weak_trends, 100*weak_trends/total_tests),
    sprintf("- Strongest correlation: %.3f (%s in %s time window)\n", 
            top_metric$correlation, top_metric$metric, top_metric$time_window),
    "\n",
    
    "🏆 Best Performing Metrics:\n"
  )
  
  # Add Top 5 metrics
  top_5 <- all_results %>%
    arrange(desc(abs(correlation))) %>%
    head(5)
  
  for(i in 1:nrow(top_5)) {
    row <- top_5[i, ]
    report <- paste0(report,
                     sprintf("%d. %s (%s): r=%.3f, p=%.3f (%s)\n",
                             i, row$metric, row$time_window, row$correlation, row$p_value,
                             ifelse(row$is_significant, "Significant", "Non-significant")))
  }
  
  # Conclusions
  report <- paste0(report,
                   "\n💡 Key Conclusions:\n",
                   if(meaningful_correlations > 0) {
                     sprintf("✅ Found %d metrics showing meaningful correlation (|r|≥0.3)\n", meaningful_correlations)
                   } else if(weak_trends > 0) {
                     sprintf("⚠️ Found %d metrics showing weak correlation trends (|r|≥0.1)\n", weak_trends)
                   } else {
                     "❌ All candidate metrics show weak correlation with OCTA clustering\n"
                   },
                   
                   if(significant_tests > 0) {
                     sprintf("✅ %d metrics achieved statistical significance\n", significant_tests)
                   } else {
                     "⚠️ No metrics achieved statistical significance, but correlation trends still have reference value\n"
                   },
                   
                   sprintf("🎯 Recommended focus: %s (%s time window), strongest correlation (r=%.3f)\n",
                           top_metric$metric, top_metric$time_window, top_metric$correlation),
                   "\n",
                   
                   "📁 Output Files:\n",
                   "- correlation_heatmap.pdf/png: All metrics correlation heatmap\n",
                   "- top_metrics_barplot.pdf/png: Top 10 metrics bar chart\n",
                   "- correlation_by_type.pdf/png: Comparison grouped by type\n",
                   "- correlation_strength_distribution.pdf/png: Correlation strength distribution\n",
                   "- single_metrics_correlation_results.csv: Detailed data\n\n",
                   
                   "Report generated at: ", Sys.time(), "\n",
                   "========================================"
  )
  
  writeLines(report, file.path(target_output_dir, "single_metrics_correlation", "single_metrics_correlation_report.txt"))
  
  cat(report)
  
  return(report)
}

# -------------------- 7. Main Execution Function --------------------

run_complete_single_metrics_analysis <- function() {
  
  cat("🚀 Starting Complete Single Metrics vs OCTA Clustering Analysis\n")
  cat("===============================================================\n")
  cat("📂 Working directory:", getwd(), "\n")
  cat("📁 Output will be saved in:", file.path(target_output_dir, "single_metrics_correlation"), "\n\n")
  
  # Step 1: Load/generate wearable data
  cat("📊 Step 1: Processing wearable data...\n")
  wearable_data <- load_wearable_data_fixed()
  
  # Step 2: Run correlation analysis using the generated data
  cat("\n🔍 Step 2: Running correlation analysis...\n")
  analysis_results <- analyze_single_metrics_correlation_fixed(wearable_data)
  
  if(nrow(analysis_results$all_results) == 0) {
    cat("❌ No correlation results generated. Please check your data.\n")
    return(NULL)
  }
  
  # Step 3: Summarize results
  cat("\n📈 Step 3: Summarizing results...\n")
  summary_results <- summarize_single_metrics_results(analysis_results$all_results)
  
  # Step 4: Create visualizations
  cat("\n🎨 Step 4: Creating visualizations...\n")
  viz_results <- create_complete_single_metrics_visualizations(
    analysis_results$all_results, 
    analysis_results$merged_data
  )
  
  # Step 5: Save results
  cat("\n💾 Step 5: Saving results...\n")
  results_output_dir <- file.path(target_output_dir, "single_metrics_correlation")
  write.csv(analysis_results$all_results, 
            file.path(results_output_dir, "single_metrics_correlation_results.csv"), 
            row.names = FALSE)
  cat("✓ Saved detailed results:", file.path(results_output_dir, "single_metrics_correlation_results.csv"), "\n")
  
  # Step 6: Generate report
  cat("\n📄 Step 6: Generating report...\n")
  report <- generate_single_metrics_report(analysis_results$all_results, summary_results)
  
  cat("\n🎉 Complete Analysis Finished!\n")
  cat("===============================\n")
  
  # Output key findings
  top_metric <- analysis_results$all_results[which.max(abs(analysis_results$all_results$correlation)), ]
  meaningful_count <- sum(abs(analysis_results$all_results$correlation) >= 0.3)
  trend_count <- sum(abs(analysis_results$all_results$correlation) >= 0.1)
  
  cat("🎯 Key Findings:\n")
  cat(sprintf("✓ Strongest correlation metric: %s (%s) r=%.3f\n", 
              top_metric$metric, top_metric$time_window, top_metric$correlation))
  cat(sprintf("✓ Meaningful correlations: %d metrics\n", meaningful_count))
  cat(sprintf("✓ Showing trends: %d metrics\n", trend_count))
  
  if(trend_count > 0) {
    cat("✅ Found correlation trends between wearable metrics and OCTA clustering!\n")
  } else {
    cat("⚠️ Overall weak correlations, may need alternative strategies\n")
  }
  
  cat("\n📂 All results saved in:", file.path(target_output_dir, "single_metrics_correlation"), "\n")
  cat("📊 Generated files:\n")
  cat("   - 4 visualization types (PDF + PNG formats)\n")
  cat("   - Detailed CSV results\n")
  cat("   - Comprehensive text report\n")
  cat("\n🗂️ Full output path structure:\n")
  cat("📁", target_output_dir, "/\n")
  cat("  ├── processed_wearable_correlation_data.csv\n")
  cat("  └── single_metrics_correlation/\n")
  cat("      ├── correlation_heatmap.pdf/png\n")
  cat("      ├── top_metrics_barplot.pdf/png\n")
  cat("      ├── correlation_by_type.pdf/png\n")
  cat("      ├── correlation_strength_distribution.pdf/png\n")
  cat("      ├── single_metrics_correlation_results.csv\n")
  cat("      └── single_metrics_correlation_report.txt\n")
  
  return(list(
    wearable_data = wearable_data,
    all_results = analysis_results$all_results,
    merged_data = analysis_results$merged_data,
    summary_results = summary_results,
    visualization_results = viz_results,
    report = report
  ))
}

# -------------------- Execute Complete Analysis --------------------

cat("🎯 Integrated Complete Single Metric Correlation Analysis\n")
cat("=========================================================\n")
cat("📂 Current working directory:", getwd(), "\n")
cat("📁 All outputs will be saved in: single_metrics_correlation/\n\n")
cat("This integrated script includes:\n")
cat("✅ Data loading and processing\n")
cat("✅ Correlation analysis\n")
cat("✅ Results summarization\n")
cat("✅ Complete visualizations (4 types)\n")
cat("✅ Comprehensive reporting\n\n")

cat("🚀 To run the complete analysis:\n")
cat("results <- run_complete_single_metrics_analysis()\n\n")

cat("📊 Expected outputs in single_metrics_correlation/:\n")
cat("1. correlation_heatmap.pdf/png\n")
cat("2. top_metrics_barplot.pdf/png\n")
cat("3. correlation_by_type.pdf/png\n")
cat("4. correlation_strength_distribution.pdf/png\n")
cat("5. single_metrics_correlation_results.csv\n")
cat("6. single_metrics_correlation_report.txt\n\n")

cat("🔍 Analysis Features:\n")
cat("• Automatic file path detection\n")
cat("• Direct data usage (no re-reading processed files)\n")
cat("• Comprehensive statistical analysis\n")
cat("• Professional visualizations\n")
cat("• Detailed reporting\n\n")

# Automatically run the analysis
cat("🚀 Starting analysis automatically...\n")
cat("📂 Output location:", file.path(target_output_dir, "single_metrics_correlation"), "\n\n")
results <- run_complete_single_metrics_analysis()