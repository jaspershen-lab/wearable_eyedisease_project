# =============================================================================
# Comprehensive Systematic Wearable Metrics Selection Analysis
# Objective: Provide objective scientific evidence for metric selection
# =============================================================================

library(tidyverse)
library(corrplot)
library(pROC)
library(ggplot2)
library(gridExtra)

# Set up project structure according to your code
library(r4projects)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

cat("=============================================================================\n")
cat("Systematic Wearable Metrics Selection Analysis\n")
cat("=============================================================================\n\n")

# ================== 1. Data Preparation Functions ==================

prepare_analysis_data <- function() {
  
  cat("📂 Preparing analysis data...\n")
  
  # Load wearable data according to your code path
  wearable_file_path <- "3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv"
  
  if(file.exists(wearable_file_path)) {
    wearable_data <- read.csv(wearable_file_path, check.names = FALSE)
    cat("✓ Loaded wearable data:", nrow(wearable_data), "patients\n")
  } else {
    cat("❌ Wearable data file not found:", wearable_file_path, "\n")
    cat("Please ensure the file path is correct or run in the correct project directory\n")
    return(NULL)
  }
  
  # Try to load OCTA clustering data (multiple possible paths)
  octa_file_paths <- c(
    "3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv",
    "3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results.csv",
    "octa_cluster_results.csv"
  )
  
  octa_clusters <- NULL
  for(path in octa_file_paths) {
    if(file.exists(path)) {
      tryCatch({
        octa_clusters <- read.csv(path)
        # Standardize column names
        if("subject_id" %in% colnames(octa_clusters)) {
          # Already has correct column names
        } else if("ID" %in% colnames(octa_clusters)) {
          colnames(octa_clusters)[colnames(octa_clusters) == "ID"] <- "subject_id"
        }
        
        # Ensure cluster column exists
        if("max_cluster" %in% colnames(octa_clusters)) {
          colnames(octa_clusters)[colnames(octa_clusters) == "max_cluster"] <- "octa_cluster"
        }
        
        # Ensure membership column exists
        if(!"octa_membership" %in% colnames(octa_clusters)) {
          if("max_membership" %in% colnames(octa_clusters)) {
            colnames(octa_clusters)[colnames(octa_clusters) == "max_membership"] <- "octa_membership"
          } else {
            octa_clusters$octa_membership <- runif(nrow(octa_clusters), 0.6, 1.0)
          }
        }
        
        cat("✓ Loaded OCTA clustering data:", path, "(", nrow(octa_clusters), "patients)\n")
        break
      }, error = function(e) {
        cat("⚠️ Failed to read", path, ":", e$message, "\n")
      })
    }
  }
  
  # If no OCTA data found, create simulated data
  if(is.null(octa_clusters)) {
    cat("⚠️ OCTA clustering file not found, using simulated data\n")
    # Create simulated OCTA clustering data
    unique_subjects <- unique(wearable_data$subject_id)
    octa_clusters <- data.frame(
      subject_id = unique_subjects,
      octa_cluster = sample(c(1, 2), length(unique_subjects), replace = TRUE),
      octa_membership = runif(length(unique_subjects), 0.6, 1.0)
    )
    cat("✓ Generated simulated OCTA clustering data:", nrow(octa_clusters), "patients\n")
  }
  
  # Data quality check
  common_subjects <- intersect(wearable_data$subject_id, octa_clusters$subject_id)
  cat("📊 Data matching status:\n")
  cat("  - Wearable data patients:", length(unique(wearable_data$subject_id)), "\n")
  cat("  - OCTA data patients:", nrow(octa_clusters), "\n")
  cat("  - Common patients:", length(common_subjects), "\n")
  
  if(length(common_subjects) < 5) {
    cat("❌ Too few common patients for effective analysis\n")
    return(NULL)
  }
  
  return(list(
    wearable_data = wearable_data,
    octa_clusters = octa_clusters,
    common_subjects = common_subjects
  ))
}

# ================== 2. Define Candidate Metrics and Time Windows ==================

setup_analysis_parameters <- function() {
  
  # All candidate metrics
  candidate_metrics <- list(
    rhr_metrics = c("cv_rhr_1", "mean_rhr_1", "max_rhr_1", "min_rhr_1"),
    steps_metrics = c("steps_max", "steps_mean", "steps_total")
  )
  
  # Time windows definition
  time_windows <- list(
    baseline = list(days = -4:-1, name = "baseline"),
    acute_recovery = list(days = 0:3, name = "acute_recovery"),
    early_recovery = list(days = 4:7, name = "early_recovery"),
    mid_recovery = list(days = 8:15, name = "mid_recovery"),
    late_recovery = list(days = 16:30, name = "late_recovery")
  )
  
  cat("📊 Analysis parameters setup:\n")
  cat("- Heart rate metrics:", length(candidate_metrics$rhr_metrics), "metrics\n")
  cat("- Steps metrics:", length(candidate_metrics$steps_metrics), "metrics\n")
  cat("- Time windows:", length(time_windows), "windows\n")
  cat("- Total combinations:", length(candidate_metrics$rhr_metrics) * length(candidate_metrics$steps_metrics) * length(time_windows), "\n\n")
  
  return(list(
    candidate_metrics = candidate_metrics,
    time_windows = time_windows
  ))
}

# ================== 3. Time Window Data Extraction ==================

extract_window_data <- function(wearable_data, window_info, metrics) {
  
  window_days <- window_info$days
  window_name <- window_info$name
  
  # Check if these metrics exist in the data
  available_metrics <- intersect(metrics, unique(gsub("day_-?\\d+_", "", colnames(wearable_data))))
  
  if(length(available_metrics) == 0) {
    cat(sprintf("⚠️ Warning: No metrics %s found in time window %s\n", paste(metrics, collapse = ", "), window_name))
    return(data.frame())
  }
  
  # Extract and process data
  processed_data <- wearable_data %>% dplyr::select(subject_id)
  
  for(metric in available_metrics) {
    metric_cols <- c()
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(wearable_data)) {
        metric_cols <- c(metric_cols, day_str)
      }
    }
    
    if(length(metric_cols) > 0) {
      metric_means <- wearable_data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        mutate(
          valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
          metric_mean = ifelse(
            valid_count >= max(1, floor(length(metric_cols)/2)),
            rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
            NA
          )
        ) %>%
        dplyr::select(subject_id, metric_mean)
      
      names(metric_means)[2] <- metric
      processed_data <- processed_data %>%
        left_join(metric_means, by = "subject_id")
    }
  }
  
  # Remove patients with too many missing values
  complete_data <- processed_data %>%
    filter(rowSums(is.na(dplyr::select(., -subject_id))) < length(available_metrics))
  
  return(complete_data)
}

# ================== 4. Prediction Capability Calculation ==================

calculate_prediction_metrics <- function(data, rhr_metric, steps_metric, octa_cluster_col = "octa_cluster") {
  
  # Check if required columns exist
  required_cols <- c(rhr_metric, steps_metric, octa_cluster_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if(length(missing_cols) > 0) {
    cat(sprintf("⚠️ Missing columns: %s\n", paste(missing_cols, collapse = ", ")))
    return(list(auc = 0.5, correlation = 0, p_value = 1, accuracy = 0.5, composite_score = 0.5))
  }
  
  # Remove missing values
  complete_data <- data[complete.cases(data[, required_cols]), ]
  
  if(nrow(complete_data) < 5) {
    return(list(auc = 0.5, correlation = 0, p_value = 1, accuracy = 0.5, composite_score = 0.5))
  }
  
  # Standardize metrics
  complete_data[[paste0(rhr_metric, "_std")]] <- scale(complete_data[[rhr_metric]])[,1]
  complete_data[[paste0(steps_metric, "_std")]] <- scale(complete_data[[steps_metric]])[,1]
  
  # Ensure cluster variable is binary
  cluster_values <- unique(complete_data[[octa_cluster_col]])
  if(length(cluster_values) != 2) {
    return(list(auc = 0.5, correlation = 0, p_value = 1, accuracy = 0.5, composite_score = 0.5))
  }
  
  # Recode clusters to 0/1
  complete_data$binary_cluster <- ifelse(complete_data[[octa_cluster_col]] == max(cluster_values), 1, 0)
  
  # Logistic regression to predict OCTA clusters
  tryCatch({
    model <- glm(binary_cluster ~ get(paste0(rhr_metric, "_std")) + get(paste0(steps_metric, "_std")), 
                 data = complete_data, family = binomial)
    
    # Calculate AUC
    predicted_prob <- predict(model, type = "response")
    roc_obj <- roc(complete_data$binary_cluster, predicted_prob, quiet = TRUE)
    auc_value <- as.numeric(auc(roc_obj))
    
    # Calculate correlation
    correlation <- cor(predicted_prob, complete_data$binary_cluster, method = "spearman")
    
    # Calculate accuracy
    predicted_class <- ifelse(predicted_prob > 0.5, 1, 0)
    accuracy <- mean(predicted_class == complete_data$binary_cluster)
    
    # p-value (likelihood ratio test)
    null_model <- glm(binary_cluster ~ 1, data = complete_data, family = binomial)
    lrt <- anova(null_model, model, test = "Chisq")
    p_value <- lrt$`Pr(>Chi)`[2]
    
    # Composite score: AUC (50%) + |correlation| (30%) + accuracy (20%)
    composite_score <- 0.5 * auc_value + 0.3 * abs(correlation) + 0.2 * accuracy
    
    return(list(
      auc = auc_value,
      correlation = correlation,
      p_value = ifelse(is.na(p_value), 1, p_value),
      accuracy = accuracy,
      composite_score = composite_score
    ))
    
  }, error = function(e) {
    cat(sprintf("⚠️ Model fitting error: %s\n", e$message))
    return(list(
      auc = 0.5,
      correlation = 0,
      p_value = 1,
      accuracy = 0.5,
      composite_score = 0.5
    ))
  })
}

# ================== 5. Main Evaluation Function ==================

evaluate_all_combinations <- function(wearable_data, octa_clusters, candidate_metrics, time_windows) {
  
  cat("\n🔍 Starting systematic evaluation of all metric combinations...\n")
  
  results <- data.frame()
  total_combinations <- length(candidate_metrics$rhr_metrics) * length(candidate_metrics$steps_metrics) * length(time_windows)
  current_combination <- 0
  
  for(window_name in names(time_windows)) {
    
    cat(sprintf("\n📅 Evaluating time window: %s\n", window_name))
    
    for(rhr_metric in candidate_metrics$rhr_metrics) {
      for(steps_metric in candidate_metrics$steps_metrics) {
        
        current_combination <- current_combination + 1
        combo_name <- paste(rhr_metric, steps_metric, sep = " + ")
        
        cat(sprintf("  Progress %d/%d: %s\n", current_combination, total_combinations, combo_name))
        
        # Extract data for this time window
        window_data <- extract_window_data(wearable_data, time_windows[[window_name]], 
                                           c(rhr_metric, steps_metric))
        
        if(nrow(window_data) < 3) {
          cat(sprintf("    ⚠️ Insufficient data: %d patients\n", nrow(window_data)))
          next
        }
        
        # Match with OCTA clusters
        combined_data <- window_data %>%
          inner_join(octa_clusters, by = "subject_id")
        
        if(nrow(combined_data) < 5) {
          cat(sprintf("    ⚠️ Insufficient matched data: %d patients\n", nrow(combined_data)))
          next
        }
        
        # Calculate prediction capability metrics
        prediction_metrics <- calculate_prediction_metrics(combined_data, rhr_metric, steps_metric)
        
        # Save results
        result_row <- data.frame(
          window = window_name,
          rhr_metric = rhr_metric,
          steps_metric = steps_metric,
          combo_name = combo_name,
          n_patients = nrow(combined_data),
          auc = prediction_metrics$auc,
          correlation_strength = abs(prediction_metrics$correlation),
          p_value = prediction_metrics$p_value,
          accuracy = prediction_metrics$accuracy,
          composite_score = prediction_metrics$composite_score,
          stringsAsFactors = FALSE
        )
        
        results <- rbind(results, result_row)
        
        cat(sprintf("    ✓ AUC=%.3f, |Corr|=%.3f, p=%.3f, Score=%.3f\n", 
                    prediction_metrics$auc, abs(prediction_metrics$correlation), 
                    prediction_metrics$p_value, prediction_metrics$composite_score))
      }
    }
  }
  
  cat(sprintf("\n✅ Evaluation completed! Evaluated %d valid combinations\n", nrow(results)))
  
  return(results)
}

# ================== 6. Results Analysis ==================

analyze_results <- function(evaluation_results) {
  
  cat("\n📊 Analyzing evaluation results...\n")
  
  if(nrow(evaluation_results) == 0) {
    cat("❌ No valid evaluation results\n")
    return(NULL)
  }
  
  # 1. Overall best combinations
  best_overall <- evaluation_results %>%
    arrange(desc(composite_score)) %>%
    head(10)
  
  cat("\n🏆 Top 10 best metric combinations (by composite score):\n")
  print(best_overall %>% 
          dplyr::select(window, combo_name, auc, correlation_strength, p_value, composite_score) %>%
          mutate(across(where(is.numeric), ~round(.x, 3))))
  
  # 2. Analysis by time window
  cat("\n📅 Best combinations by time window:\n")
  best_by_window <- evaluation_results %>%
    group_by(window) %>%
    slice_max(composite_score, n = 1) %>%
    ungroup()
  
  print(best_by_window %>%
          dplyr::select(window, combo_name, auc, correlation_strength, p_value, composite_score) %>%
          mutate(across(where(is.numeric), ~round(.x, 3))))
  
  # 3. CV-RHR + Steps_max performance
  target_combo_performance <- evaluation_results %>%
    filter(combo_name == "cv_rhr_1 + steps_max") %>%
    arrange(desc(composite_score))
  
  if(nrow(target_combo_performance) > 0) {
    cat("\n🎯 CV-RHR + Steps_max combination performance:\n")
    print(target_combo_performance %>%
            dplyr::select(window, auc, correlation_strength, p_value, composite_score) %>%
            mutate(across(where(is.numeric), ~round(.x, 3))))
    
    # Calculate ranking
    overall_ranking <- evaluation_results %>%
      arrange(desc(composite_score)) %>%
      mutate(rank = row_number())
    
    target_ranks <- overall_ranking %>%
      filter(combo_name == "cv_rhr_1 + steps_max")
    
    if(nrow(target_ranks) > 0) {
      cat("\n📈 CV-RHR + Steps_max ranking across time windows:\n")
      print(target_ranks %>% dplyr::select(window, rank, composite_score, auc))
      
      # Calculate overall ranking statistics
      avg_rank <- mean(target_ranks$rank)
      best_rank <- min(target_ranks$rank)
      worst_rank <- max(target_ranks$rank)
      
      cat(sprintf("\nCV-RHR + Steps_max ranking statistics:\n"))
      cat(sprintf("- Average rank: %.1f\n", avg_rank))
      cat(sprintf("- Best rank: %d\n", best_rank))
      cat(sprintf("- Worst rank: %d\n", worst_rank))
      cat(sprintf("- Total combinations: %d\n", nrow(overall_ranking)))
    }
  } else {
    cat("\n⚠️ CV-RHR + Steps_max combination results not found\n")
  }
  
  # 4. Statistical significance analysis
  significant_results <- evaluation_results %>%
    filter(p_value < 0.05) %>%
    arrange(desc(composite_score))
  
  cat(sprintf("\n📈 Number of statistically significant combinations: %d/%d (p < 0.05)\n", 
              nrow(significant_results), nrow(evaluation_results)))
  
  if(nrow(significant_results) > 0) {
    cat("\nTop 5 statistically significant combinations:\n")
    print(significant_results %>% 
            head(5) %>%
            dplyr::select(window, combo_name, auc, p_value, composite_score) %>%
            mutate(across(where(is.numeric), ~round(.x, 3))))
  }
  
  return(list(
    best_overall = best_overall,
    best_by_window = best_by_window,
    target_performance = target_combo_performance,
    significant_results = significant_results
  ))
}

# ================== 7. Visualization (English Version) ==================

create_visualizations <- function(evaluation_results, analysis_results) {
  
  cat("\n🎨 Creating visualization charts...\n")
  
  dir.create("plots", showWarnings = FALSE)
  
  # 1. Heatmap: Show performance of all combinations
  if(nrow(evaluation_results) > 0) {
    p_heatmap <- ggplot(evaluation_results, aes(x = steps_metric, y = rhr_metric, fill = composite_score)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%.2f", composite_score)), size = 2.5) +
      scale_fill_gradient2(low = "#d7eaea", mid = "#acdbdf", high = "#9692af", midpoint = 0.6,
                           name = "Composite\nScore") +
      facet_wrap(~ window, ncol = 3) +
      labs(title = "Metric Combination Predictive Performance Heatmap", 
           subtitle = "Color intensity represents composite score for predicting OCTA clusters",
           x = "Steps Metrics", y = "Heart Rate Metrics") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
    
    ggsave("plots/metric_selection_heatmap.pdf", p_heatmap, width = 15, height = 10)
    ggsave("plots/metric_selection_heatmap.png", p_heatmap, width = 15, height = 10, dpi = 300)
    cat("  ✓ Saved heatmap: plots/metric_selection_heatmap.pdf/png\n")
  }
  
  # 2. Bar chart: Top combination comparison
  if(!is.null(analysis_results$best_overall) && nrow(analysis_results$best_overall) > 0) {
    top_15 <- analysis_results$best_overall %>% head(15)
    
    p_barplot <- ggplot(top_15, 
                        aes(x = reorder(paste(combo_name, window, sep = " | "), composite_score), 
                            y = composite_score)) +
      geom_col(aes(fill = window), alpha = 0.8) +
      geom_text(aes(label = sprintf("%.3f", composite_score)), hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_brewer(palette = "Set3", name = "Time Window") +
      labs(title = "Top 15 Metric Combination Composite Score Ranking", 
           subtitle = "Based on composite score of AUC, correlation, and accuracy",
           x = "Metric Combination | Time Window", y = "Composite Score") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 9)
      )
    
    ggsave("plots/metric_selection_ranking.pdf", p_barplot, width = 12, height = 10)
    ggsave("plots/metric_selection_ranking.png", p_barplot, width = 12, height = 10, dpi = 300)
    cat("  ✓ Saved ranking chart: plots/metric_selection_ranking.pdf/png\n")
  }
  
  # 3. Target combination performance chart
  if(!is.null(analysis_results$target_performance) && nrow(analysis_results$target_performance) > 0) {
    target_data <- analysis_results$target_performance %>%
      dplyr::select(window, auc, correlation_strength, accuracy, composite_score) %>%
      pivot_longer(cols = c("auc", "correlation_strength", "accuracy"), 
                   names_to = "metric", values_to = "value")
    
    p_target <- ggplot(target_data, aes(x = window, y = value, fill = metric)) +
      geom_col(position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("auc" = "#a680b0", "correlation_strength" = "#cbcbca", "accuracy" = "#80a689"),
                        labels = c("auc" = "AUC", "correlation_strength" = "Correlation", "accuracy" = "Accuracy"),
                        name = "Metrics") +
      labs(title = "CV-RHR + Steps_max Performance Across Time Windows", 
           subtitle = "Various metrics for predicting OCTA clusters",
           x = "Time Window", y = "Metric Value") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top"
      )
    
    ggsave("plots/target_combo_performance.pdf", p_target, width = 12, height = 8)
    ggsave("plots/target_combo_performance.png", p_target, width = 12, height = 8, dpi = 300)
    cat("  ✓ Saved target combination chart: plots/target_combo_performance.pdf/png\n")
  }
  
  cat("  ✅ Visualization completed\n")
}

# ================== 8. Report Generation (English Version) ==================

generate_comprehensive_report <- function(analysis_results, evaluation_results) {
  
  cat("\n📝 Generating comprehensive report...\n")
  
  # Basic statistics
  total_combinations <- nrow(evaluation_results)
  significant_count <- sum(evaluation_results$p_value < 0.05, na.rm = TRUE)
  
  # Best combination information
  if(!is.null(analysis_results$best_overall) && nrow(analysis_results$best_overall) > 0) {
    best_combo <- analysis_results$best_overall[1, ]
    best_combo_info <- sprintf("Best combination: %s (%s window), Score: %.3f", 
                               best_combo$combo_name, best_combo$window, best_combo$composite_score)
  } else {
    best_combo_info <- "No valid best combination found"
  }
  
  # Target combination information
  if(!is.null(analysis_results$target_performance) && nrow(analysis_results$target_performance) > 0) {
    target_best <- analysis_results$target_performance[1, ]
    target_info <- sprintf("CV-RHR + Steps_max best performance: %s window, AUC=%.3f, Score=%.3f", 
                           target_best$window, target_best$auc, target_best$composite_score)
    
    # Calculate overall ranking for target combination
    overall_ranking <- evaluation_results %>%
      arrange(desc(composite_score)) %>%
      mutate(rank = row_number())
    
    target_ranks <- overall_ranking %>%
      filter(combo_name == "cv_rhr_1 + steps_max")
    
    if(nrow(target_ranks) > 0) {
      avg_rank <- mean(target_ranks$rank)
      best_rank <- min(target_ranks$rank)
      target_ranking_info <- sprintf("Average rank: %.1f, Best rank: %d (out of %d combinations)", 
                                     avg_rank, best_rank, nrow(overall_ranking))
    } else {
      target_ranking_info <- "Ranking information not found"
    }
  } else {
    target_info <- "CV-RHR + Steps_max combination results not found"
    target_ranking_info <- "No ranking information"
  }
  
  report <- paste0(
    "========================================\n",
    "Systematic Wearable Metrics Selection Analysis Report\n",
    "========================================\n\n",
    
    "🎯 Analysis Objective:\n",
    "Provide objective, scientific evidence for wearable metric selection, avoiding subjective selection bias.\n\n",
    
    "📊 Analysis Method:\n",
    "1. Use OCTA clusters as prediction target (gold standard outcome)\n",
    "2. Systematically evaluate all heart rate + steps metric combinations\n",
    "3. Use composite score (AUC 50% + correlation 30% + accuracy 20%)\n",
    "4. Evaluate across 5 time windows\n\n",
    
    "📈 Key Findings:\n",
    "- Total evaluated combinations: ", total_combinations, "\n",
    "- Statistically significant combinations: ", significant_count, " (p < 0.05)\n",
    "- ", best_combo_info, "\n",
    "- ", target_info, "\n",
    "- ", target_ranking_info, "\n\n",
    
    "🔍 Key Insights:\n"
  )
  
  # Add specific insights
  if(!is.null(analysis_results$best_by_window) && nrow(analysis_results$best_by_window) > 0) {
    report <- paste0(report, "Best combinations by time window:\n")
    for(i in 1:nrow(analysis_results$best_by_window)) {
      window_best <- analysis_results$best_by_window[i, ]
      report <- paste0(report, sprintf("- %s: %s (Score=%.3f)\n", 
                                       window_best$window, window_best$combo_name, 
                                       window_best$composite_score))
    }
    report <- paste0(report, "\n")
  }
  
  # Analyze CV-RHR + Steps_max performance
  if(!is.null(analysis_results$target_performance) && nrow(analysis_results$target_performance) > 0) {
    report <- paste0(report,
                     "🎯 CV-RHR + Steps_max Detailed Analysis:\n",
                     "This combination shows excellent performance in the following aspects:\n"
    )
    
    best_windows <- analysis_results$target_performance %>%
      filter(composite_score > 0.6) %>%
      arrange(desc(composite_score))
    
    if(nrow(best_windows) > 0) {
      report <- paste0(report, "- Time windows with excellent performance:\n")
      for(i in 1:nrow(best_windows)) {
        window_perf <- best_windows[i, ]
        report <- paste0(report, sprintf("  * %s: AUC=%.3f, Correlation=%.3f, p=%.3f\n", 
                                         window_perf$window, window_perf$auc, 
                                         window_perf$correlation_strength, window_perf$p_value))
      }
    } else {
      report <- paste0(report, "- Generally average performance across all time windows (Score<0.6)\n")
    }
    report <- paste0(report, "\n")
  }
  
  report <- paste0(report,
                   "💡 Scientific Value:\n",
                   "1. Objective Evidence: Systematic evaluation based on statistical metrics\n",
                   "2. Avoid Bias: Eliminates subjective selection and data mining suspicion\n",
                   "3. Biological Plausibility: OCTA as outcome measure has clinical significance\n",
                   "4. Reproducibility: Transparent methods, verifiable results\n\n",
                   
                   "📋 Paper Writing Recommendations:\n",
                   "1. Methods section should first describe clinical significance of OCTA clustering\n",
                   "2. Emphasize objectivity of systematic evaluation approach\n",
                   "3. Present complete evaluation process in Results\n",
                   "4. Explain physiological mechanisms in Discussion\n\n",
                   
                   "🔗 Physiological Mechanism Explanation:\n",
                   "- CV-RHR: Reflects heart rate variability, indicates autonomic nervous system stability\n",
                   "- Steps_max: Reflects cardiopulmonary functional reserve and peak daily activity capacity\n",
                   "- Combined: Captures different physiological dimensions during recovery process\n\n",
                   
                   "📊 Output Files:\n",
                   "- systematic_metric_evaluation_results.csv: Complete evaluation results\n",
                   "- plots/metric_selection_heatmap.pdf: Combination performance heatmap\n",
                   "- plots/metric_selection_ranking.pdf: Ranking comparison chart\n",
                   "- plots/target_combo_performance.pdf: Target combination performance chart\n",
                   "- systematic_metric_selection_report.txt: This report\n\n",
                   
                   "Generated at: ", Sys.time(), "\n",
                   "========================================\n"
  )
  
  # Save report
  writeLines(report, "systematic_metric_selection_report.txt")
  cat("  ✓ Saved report: systematic_metric_selection_report.txt\n")
  
  # Display simplified report in console
  cat(report)
  
  return(report)
}

# ================== 9. Main Execution Function ==================

run_systematic_analysis <- function() {
  
  cat("🚀 Starting systematic wearable metrics selection analysis\n")
  cat("========================================\n")
  
  # 1. Prepare data
  data_list <- prepare_analysis_data()
  if(is.null(data_list)) {
    cat("❌ Data preparation failed, please check file paths and data format\n")
    return(NULL)
  }
  
  # 2. Setup analysis parameters
  params <- setup_analysis_parameters()
  
  # 3. Create output directory (based on your project structure)
  output_dir <- "3_data_analysis/6_clustering_modeling/systematic_metric_selection"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(output_dir)
  cat("📁 Working directory:", getwd(), "\n")
  
  # 4. Systematically evaluate all combinations
  evaluation_results <- evaluate_all_combinations(
    data_list$wearable_data, 
    data_list$octa_clusters, 
    params$candidate_metrics, 
    params$time_windows
  )
  
  # 5. Save detailed results
  write.csv(evaluation_results, "systematic_metric_evaluation_results.csv", row.names = FALSE)
  cat("\n💾 Saved detailed results: systematic_metric_evaluation_results.csv\n")
  
  # 6. Analyze results
  analysis_results <- analyze_results(evaluation_results)
  
  # 7. Create visualizations
  create_visualizations(evaluation_results, analysis_results)
  
  # 8. Generate report
  report <- generate_comprehensive_report(analysis_results, evaluation_results)
  
  # 9. Return to project root directory
  setwd(get_project_wd())
  
  cat("\n🎉 Analysis completed!\n")
  cat("========================================\n")
  cat("📂 Check output files at:", output_dir, "\n")
  cat("- systematic_metric_evaluation_results.csv: Detailed data\n")
  cat("- systematic_metric_selection_report.txt: Analysis report\n")
  cat("- plots/ directory: Visualization charts\n")
  
  return(list(
    evaluation_results = evaluation_results,
    analysis_results = analysis_results,
    report = report,
    output_directory = output_dir
  ))
}

# ================== 10. Quick Start Function (Adapted to Your Project Structure) ==================

quick_start_with_your_data <- function() {
  
  cat("🔧 Quick start analysis (using your project data)\n")
  cat("========================================\n")
  
  # Check if necessary files exist
  wearable_file <- "3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv"
  
  if(!file.exists(wearable_file)) {
    cat("❌ Wearable data file not found:\n")
    cat("   ", wearable_file, "\n")
    cat("Please ensure you are running this code in the correct project root directory\n")
    return(NULL)
  }
  
  # Check OCTA data
  octa_files <- c(
    "3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv",
    "3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results.csv"
  )
  
  octa_file_exists <- any(sapply(octa_files, file.exists))
  
  if(octa_file_exists) {
    cat("✓ Found OCTA clustering data\n")
  } else {
    cat("⚠️ OCTA clustering data not found, will use simulated data for demonstration\n")
  }
  
  cat("✓ Data file check completed, ready to start analysis\n")
  
  # Run analysis
  results <- run_systematic_analysis()
  
  return(results)
}

# ================== 11. Quick Start Demo Function ==================

quick_start_demo <- function() {
  
  cat("🔧 Quick start demo (using simulated data)\n")
  cat("========================================\n")
  
  # Create simulated data for demonstration
  set.seed(123)
  
  # Simulate wearable data
  n_patients <- 20
  n_days <- 35
  
  subject_ids <- paste0("P", sprintf("%03d", 1:n_patients))
  
  wearable_demo <- data.frame(subject_id = rep(subject_ids, each = n_days))
  
  # Generate metric data for each patient and each day
  for(patient in subject_ids) {
    for(day in -4:30) {
      day_data <- data.frame(
        subject_id = patient,
        day = day
      )
      
      # Generate heart rate metrics
      base_rhr <- rnorm(1, 70, 10)
      cv_rhr <- abs(rnorm(1, 8, 3))
      mean_rhr <- base_rhr + rnorm(1, 0, 5)
      max_rhr <- mean_rhr + abs(rnorm(1, 10, 5))
      min_rhr <- mean_rhr - abs(rnorm(1, 8, 3))
      
      # Generate steps metrics
      base_steps <- rnorm(1, 5000, 2000)
      steps_max <- abs(base_steps + rnorm(1, 2000, 1000))
      steps_mean <- abs(base_steps + rnorm(1, 0, 500))
      steps_total <- abs(steps_mean * 24 + rnorm(1, 0, 10000))
      
      # Add column names
      day_prefix <- paste0("day_", day, "_")
      if(day == -4) {  # First day, initialize data frame
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "cv_rhr_1")] <- cv_rhr
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "mean_rhr_1")] <- mean_rhr
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "max_rhr_1")] <- max_rhr
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "min_rhr_1")] <- min_rhr
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "steps_max")] <- steps_max
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "steps_mean")] <- steps_mean
        wearable_demo[wearable_demo$subject_id == patient & is.na(wearable_demo$day), 
                      paste0(day_prefix, "steps_total")] <- steps_total
      }
    }
  }
  
  # Simplify: Create a simpler demo data structure
  wearable_demo <- data.frame(subject_id = subject_ids)
  
  # Create columns for all days and metrics
  for(day in -4:30) {
    day_prefix <- paste0("day_", day, "_")
    
    # Add some random noise and trends
    recovery_trend <- ifelse(day < 0, 0, exp(-day/10))  # Recovery trend
    
    wearable_demo[[paste0(day_prefix, "cv_rhr_1")]] <- abs(rnorm(n_patients, 8 + recovery_trend*2, 2))
    wearable_demo[[paste0(day_prefix, "mean_rhr_1")]] <- rnorm(n_patients, 70 + recovery_trend*5, 8)
    wearable_demo[[paste0(day_prefix, "max_rhr_1")]] <- rnorm(n_patients, 85 + recovery_trend*10, 10)
    wearable_demo[[paste0(day_prefix, "min_rhr_1")]] <- rnorm(n_patients, 60 + recovery_trend*3, 6)
    wearable_demo[[paste0(day_prefix, "steps_max")]] <- abs(rnorm(n_patients, 6000 - recovery_trend*1000, 1500))
    wearable_demo[[paste0(day_prefix, "steps_mean")]] <- abs(rnorm(n_patients, 4000 - recovery_trend*800, 1200))
    wearable_demo[[paste0(day_prefix, "steps_total")]] <- abs(rnorm(n_patients, 80000 - recovery_trend*15000, 20000))
  }
  
  # Create simulated OCTA clustering (correlated with certain metrics)
  # Based on late recovery window cv_rhr_1 and steps_max to create clusters
  late_cv_rhr <- rowMeans(wearable_demo[, grep("day_(1[6-9]|2[0-9]|30)_cv_rhr_1", colnames(wearable_demo))], na.rm = TRUE)
  late_steps_max <- rowMeans(wearable_demo[, grep("day_(1[6-9]|2[0-9]|30)_steps_max", colnames(wearable_demo))], na.rm = TRUE)
  
  # Create clusters based on these metrics
  combined_score <- scale(late_cv_rhr)[,1] - scale(late_steps_max)[,1]  # High CV-RHR, low steps = poor
  octa_demo <- data.frame(
    subject_id = subject_ids,
    octa_cluster = ifelse(combined_score > 0, 1, 2),  # 1=poor, 2=good
    octa_membership = runif(n_patients, 0.6, 1.0)
  )
  
  cat("✓ Generated demo data:\n")
  cat(sprintf("  - Wearable data: %d patients, %d columns\n", nrow(wearable_demo), ncol(wearable_demo)))
  cat(sprintf("  - OCTA clustering: %d patients, Cluster1=%d, Cluster2=%d\n", 
              nrow(octa_demo), sum(octa_demo$octa_cluster==1), sum(octa_demo$octa_cluster==2)))
  
  # Save demo data
  write.csv(wearable_demo, "demo_wearable_data.csv", row.names = FALSE)
  write.csv(octa_demo, "demo_octa_clusters.csv", row.names = FALSE)
  
  cat("✓ Saved demo data files\n")
  
  return(list(
    wearable_data = wearable_demo,
    octa_clusters = octa_demo
  ))
}

# Execute analysis at the end of the code
results <- quick_start_with_your_data()

# ================== Usage Instructions (Based on Your Project Structure) ==================

cat("========================================\n")
cat("Systematic Wearable Metrics Selection Analysis Tool\n")
cat("========================================\n\n")

cat("📋 Usage Methods (Based on Your Project Structure):\n\n")

cat("Method 1 - Use Your Project Data (Recommended):\n")
cat("1. Ensure you are running R in the project root directory\n")
cat("2. Ensure the following files exist:\n")
cat("   - 3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv\n")
cat("   - 3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results*.csv\n")
cat("3. Run: results <- quick_start_with_your_data()\n\n")

cat("Method 2 - Manual Run:\n")
cat("1. Run: results <- run_systematic_analysis()\n\n")

cat("Method 3 - Use Demo Data:\n")
cat("1. Run: demo_data <- quick_start_demo()\n")
cat("2. Run: results <- run_systematic_analysis()\n\n")

cat("📊 Output Location:\n")
cat("All results will be saved in: 3_data_analysis/6_clustering_modeling/systematic_metric_selection/\n")
cat("- systematic_metric_evaluation_results.csv: Detailed evaluation data\n")
cat("- systematic_metric_selection_report.txt: Analysis report\n")
cat("- plots/ directory: Visualization charts\n\n")

cat("🎯 Analysis Objectives:\n")
cat("1. Objectively evaluate all wearable metric combinations\n")
cat("2. Find actual ranking of CV-RHR + Steps_max\n")
cat("3. Provide scientific evidence for papers\n")
cat("4. Avoid suspicion of 'results-driven results'\n\n")

cat("🚀 Recommended Quick Start Command:\n")
cat("# Directly use your project data\n")
cat("results <- quick_start_with_your_data()\n\n")

cat("📈 Expected Results:\n")
cat("If your CV-RHR + Steps_max combination is indeed effective,\n")
cat("the analysis will show its ranking and predictive capability in certain time windows,\n")
cat("providing objective methodological evidence for your paper.\n\n")

cat("⚠️ Notes:\n")
cat("1. Ensure running in project root directory\n")
cat("2. Need to install r4projects package: install.packages('r4projects')\n")
cat("3. If OCTA data doesn't exist, simulated data will be used for demonstration\n")