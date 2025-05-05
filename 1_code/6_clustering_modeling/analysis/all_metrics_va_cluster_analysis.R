# -----------------------------------------------------
# Analyzing Vision Improvements by Clusters for PPV and Cataract groups
# Focusing on vision_improvement_1w and vision_improvement_1m
# -----------------------------------------------------
library(r4projects)
library(tidyverse)
library(ggplot2)
library(rstatix)      # For statistical testing
library(ggpubr)       # For creating publication-ready plots
library(corrplot)     # For correlation visualizations
library(gtsummary)    # For creating summary tables
library(ggstatsplot)  # For statistical plots with significance indicators
library(rlang)        # For proper handling of non-standard evaluation
setwd(get_project_wd())
rm(list = ls())

# -----------------------------------------------------
# 1. Load the cluster results and vision data
# -----------------------------------------------------
# Assuming these files were created in your previous clustering analysis
ppv_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/ppv_cluster_results_all_metrics.csv", check.names = FALSE)
cat_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/cataract_cluster_results_all_metrics.csv", check.names = FALSE)
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Extract the subject_id and cluster information
ppv_cluster_info <- data.frame(
  ID = ppv_clusters$subject_id,
  cluster = ppv_clusters$max_cluster,
  membership = ppv_clusters$max_membership
)

cat_cluster_info <- data.frame(
  ID = cat_clusters$subject_id,
  cluster = cat_clusters$max_cluster,
  membership = cat_clusters$max_membership
)

# -----------------------------------------------------
# 2. Create vision dataset with both 1-week and 1-month post-surgery data
# -----------------------------------------------------
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
    post_vision_1w = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1w,   # Right eye post-surgery 1 week
      surgery_eye_1 == 1 ~ os_corrected_1w,   # Left eye post-surgery 1 week
      surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # Both eyes post-surgery 1 week
      TRUE ~ NA_real_
    ),
    post_vision_1m = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1m,   # Right eye post-surgery 1 month
      surgery_eye_1 == 1 ~ os_corrected_1m,   # Left eye post-surgery 1 month
      surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,   # Both eyes post-surgery 1 month
      TRUE ~ NA_real_
    ),
    vision_improvement_1w = post_vision_1w - pre_vision,
    vision_improvement_1m = post_vision_1m - pre_vision,
    vision_improved_1m = if_else(vision_improvement_1m >= 0, 1, 0),  # Binary improvement indicator
    vision_improved_factor_1m = factor(vision_improved_1m, 
                                       levels = c(0, 1), 
                                       labels = c("NoImprovement", "Improved"))  # Factor version
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision_1w, post_vision_1m,
                vision_improvement_1w, vision_improvement_1m, 
                vision_improved_1m, vision_improved_factor_1m,
                age, gender)

# -----------------------------------------------------
# 3. Merge cluster information with vision improvement data
# -----------------------------------------------------
# For PPV group
ppv_vision_analysis <- vision_data %>%
  inner_join(ppv_cluster_info, by = "ID") %>%
  mutate(surgery_type = "PPV")

# For Cataract group
cat_vision_analysis <- vision_data %>%
  inner_join(cat_cluster_info, by = "ID") %>%
  mutate(surgery_type = "Cataract")

# Combine data for overall analysis if needed
all_vision_analysis <- bind_rows(ppv_vision_analysis, cat_vision_analysis)

# -----------------------------------------------------
# 4. Create visualization using ggbetweenstats for vision improvement
# -----------------------------------------------------
# Function to create plots with ggbetweenstats
create_betweenstats_plot <- function(data, param_col, title_prefix, surgery_type) {
  # Convert cluster to factor for better plotting
  data$cluster <- as.factor(data$cluster)
  
  # Create a more readable parameter name
  param_name <- param_col
  
  # Create the plot with ggbetweenstats using formula approach
  p <- ggstatsplot::ggbetweenstats(
    data = data,
    x = cluster,
    y = !!sym(param_col), # Use !!sym() for proper evaluation
    type = "parametric",
    pairwise.display = "all",  # Show all pairwise comparisons, not just significant ones
    p.adjust.method = "holm",
    effsize.type = "unbiased",
    results.subtitle = TRUE,
    xlab = "Cluster",
    ylab = paste(param_name, ""),
    title = paste0(title_prefix, " - ", param_name),
    subtitle = paste("Surgery Type:", surgery_type),
    centrality.plotting = TRUE,
    centrality.type = "parametric",
    centrality.point.args = list(size = 5, color = "darkred"),
    point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                      alpha = 0.4, size = 3, stroke = 0),
    boxplot.args = list(width = 0.3, alpha = 0.2),
    violin.args = list(width = 0.5, alpha = 0.2),
    ggsignif.args = list(textsize = 3, tip_length = 0.01),
    ggtheme = ggplot2::theme_bw(),  # Use theme_bw() instead of theme_ggstatsplot
    package = "RColorBrewer",
    palette = "Dark2"
  )
  
  return(p)
}

# Create directory for plots
dir.create("3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/plots", 
           recursive = TRUE, showWarnings = FALSE)

# Generate plots for PPV vision improvement
vision_params <- c("vision_improvement_1w", "vision_improvement_1m")

for (param in vision_params) {
  p <- create_betweenstats_plot(
    ppv_vision_analysis, 
    param, 
    "Vision Improvement", 
    "PPV"
  )
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/vision_analysis/plots/ppv_", 
           param, 
           ".pdf"),
    p,
    width = 10,
    height = 8
  )
}

# Generate plots for cataract vision improvement
for (param in vision_params) {
  p <- create_betweenstats_plot(
    cat_vision_analysis, 
    param, 
    "Vision Improvement", 
    "Cataract"
  )
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/vision_analysis/plots/cat_", 
           param, 
           ".pdf"),
    p,
    width = 10,
    height = 8
  )
}

# -----------------------------------------------------
# 5. Statistical analysis of differences between clusters
# -----------------------------------------------------
# Function to test for differences between clusters
test_cluster_differences <- function(data, param_cols) {
  results_list <- list()
  
  for (param in param_cols) {
    # ANOVA to test for differences between clusters
    model_formula <- as.formula(paste(param, "~ cluster"))
    anova_result <- data %>%
      anova_test(model_formula)
    
    # Add results to the list
    results_list[[paste0(param, "_anova")]] <- anova_result
    
    # If ANOVA is significant, add post-hoc tests
    if (anova_result$p < 0.05) {
      posthoc <- data %>%
        pairwise_t_test(
          model_formula,
          p.adjust.method = "holm"
        )
      
      results_list[[paste0(param, "_posthoc")]] <- posthoc
    }
  }
  
  return(results_list)
}

# Test differences for PPV vision improvement parameters
ppv_vision_tests <- test_cluster_differences(
  ppv_vision_analysis,
  vision_params
)

# Test differences for cataract vision improvement parameters
cat_vision_tests <- test_cluster_differences(
  cat_vision_analysis,
  vision_params
)

# Create summary table of significant results
create_significant_summary <- function(test_results, params, output_file) {
  sig_results <- data.frame(
    Parameter = character(),
    P_Value = numeric(),
    Significant_Pairs = character(),
    stringsAsFactors = FALSE
  )
  
  for (param in params) {
    anova_result <- test_results[[paste0(param, "_anova")]]
    
    if (!is.null(anova_result) && anova_result$p < 0.05) {
      param_name <- param
      
      # Get post-hoc results if available
      posthoc <- test_results[[paste0(param, "_posthoc")]]
      sig_pairs <- ""
      
      if (!is.null(posthoc)) {
        sig_pairs_df <- posthoc %>% 
          filter(p.adj < 0.05) %>%
          mutate(pair = paste(group1, "vs", group2))
        
        if (nrow(sig_pairs_df) > 0) {
          sig_pairs <- paste(sig_pairs_df$pair, collapse = "; ")
        }
      }
      
      # Add to summary
      sig_results <- sig_results %>%
        rbind(data.frame(
          Parameter = param_name,
          P_Value = anova_result$p,
          Significant_Pairs = sig_pairs,
          stringsAsFactors = FALSE
        ))
    }
  }
  
  # Sort by p-value
  sig_results <- sig_results %>%
    arrange(P_Value)
  
  # Save to file
  write.csv(sig_results, output_file, row.names = FALSE)
  
  return(sig_results)
}

# Create summary tables
ppv_sig_summary <- create_significant_summary(
  ppv_vision_tests,
  vision_params,
  "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/ppv_vision_significant_results.csv"
)

cat_sig_summary <- create_significant_summary(
  cat_vision_tests,
  vision_params,
  "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/cat_vision_significant_results.csv"
)

# -----------------------------------------------------
# 6. Correlation between cluster membership and vision improvement
# -----------------------------------------------------
# Function to analyze correlation between membership and improvement
analyze_membership_correlation <- function(data, param_cols) {
  results <- data.frame(
    Parameter = character(),
    Rho = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (param in param_cols) {
    # Spearman correlation test
    corr_test <- cor.test(data$membership, data[[param]], 
                          method = "spearman", exact = FALSE)
    
    # Add results to dataframe
    results <- results %>%
      rbind(data.frame(
        Parameter = param,
        Rho = corr_test$estimate,
        P_value = corr_test$p.value,
        stringsAsFactors = FALSE
      ))
  }
  
  # Sort by absolute correlation strength
  results <- results %>%
    arrange(P_value)
  
  return(results)
}

# Analyze correlation for PPV vision improvement
ppv_membership_corr <- analyze_membership_correlation(
  ppv_vision_analysis,
  vision_params
)

# Analyze correlation for cataract vision improvement
cat_membership_corr <- analyze_membership_correlation(
  cat_vision_analysis,
  vision_params
)

# Save correlation results
write.csv(
  ppv_membership_corr,
  "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/ppv_membership_correlation_vision.csv",
  row.names = FALSE
)

write.csv(
  cat_membership_corr,
  "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/cat_membership_correlation_vision.csv",
  row.names = FALSE
)

# -----------------------------------------------------
# 7. Additional Analysis: Vision Improvement by Cluster (barplots)
# -----------------------------------------------------
# Function to create barplots showing mean vision improvement by cluster
create_vision_barplot <- function(data, param, title, surgery_type) {
  # Make sure cluster is a factor for proper handling in ggplot
  data$cluster <- as.factor(data$cluster)
  
  # Compute mean and standard error by cluster
  summary_data <- data %>%
    group_by(cluster) %>%
    summarise(
      mean_improvement = mean(!!sym(param), na.rm = TRUE),
      se_improvement = sd(!!sym(param), na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  
  # Create barplot
  p <- ggplot(summary_data, aes(x = cluster, y = mean_improvement, fill = cluster)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    geom_errorbar(aes(ymin = mean_improvement - se_improvement, 
                      ymax = mean_improvement + se_improvement),
                  width = 0.25, position = position_dodge(0.7)) +
    labs(
      title = paste0(title, " by Cluster"),
      subtitle = paste("Surgery Type:", surgery_type),
      x = "Cluster",
      y = paste("Mean", param),
      caption = "Error bars represent standard error of the mean"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_fill_brewer(palette = "Dark2")
  
  return(p)
}

# Create barplots for PPV
for (param in vision_params) {
  p <- create_vision_barplot(
    ppv_vision_analysis,
    param,
    param,
    "PPV"
  )
  
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/vision_analysis/plots/ppv_", 
           param, "_barplot.pdf"),
    p,
    width = 8,
    height = 6
  )
}

# Create barplots for Cataract
for (param in vision_params) {
  p <- create_vision_barplot(
    cat_vision_analysis,
    param,
    param,
    "Cataract"
  )
  
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/vision_analysis/plots/cat_", 
           param, "_barplot.pdf"),
    p,
    width = 8,
    height = 6
  )
}

# -----------------------------------------------------
# 8. Additional Analysis: Classification using logistic regression
# -----------------------------------------------------
# Function to build a logistic regression model for predicting vision improvement
build_logistic_model <- function(data, output_file_prefix) {
  # Prepare data - focus on 1-month improvement (binary outcome)
  model_data <- data %>%
    # Remove rows with missing values
    filter(!is.na(vision_improved_1m), !is.na(membership), !is.na(cluster))
  
  # Convert cluster to factor
  model_data$cluster <- as.factor(model_data$cluster)
  
  # Check if there's enough variation in the outcome variable
  outcome_table <- table(model_data$vision_improved_1m)
  
  # Save basic counts to a file
  sink(paste0(output_file_prefix, "_outcome_distribution.txt"))
  cat("Distribution of vision improvement outcomes:\n")
  print(outcome_table)
  cat("\nPercentage improved: ", 
      round(outcome_table["1"] / sum(outcome_table) * 100, 1), 
      "%\n", sep="")
  sink()
  
  # If there's not enough variation, return early with a warning
  if(length(outcome_table) < 2 || min(outcome_table) < 3) {
    warning("Not enough variation in outcome variable for logistic regression.")
    
    sink(paste0(output_file_prefix, "_logistic_model_summary.txt"))
    cat("WARNING: Not enough variation in outcome variable for logistic regression.\n")
    cat("At least 3 observations in each outcome category (improved/not improved) are needed.\n")
    cat("\nOutcome distribution:\n")
    print(outcome_table)
    sink()
    
    # Return NULL instead of a model
    return(list(model = NULL, odds_ratios = NULL, error = "Insufficient variation"))
  }
  
  # Try to fit logistic regression model
  tryCatch({
    # Fit model
    model <- glm(vision_improved_1m ~ cluster + membership, 
                 family = binomial(link = "logit"), 
                 data = model_data)
    
    # Create summary table
    model_summary <- summary(model)
    
    # Save model summary
    sink(paste0(output_file_prefix, "_logistic_model_summary.txt"))
    print(model_summary)
    sink()
    
    # Try to calculate confidence intervals
    ci_result <- tryCatch({
      confint(model)
    }, error = function(e) {
      # If confidence intervals fail, return NA matrix
      warning("Could not calculate confidence intervals: ", e$message)
      matrix(NA, nrow = length(coef(model)), ncol = 2,
             dimnames = list(names(coef(model)), c("2.5 %", "97.5 %")))
    })
    
    # Calculate odds ratios
    or_data <- data.frame(
      Variable = names(coef(model)),
      OR = exp(coef(model)),
      Lower_CI = exp(ci_result)[,1],
      Upper_CI = exp(ci_result)[,2],
      P_value = coef(summary(model))[,4]
    )
    
    # Save odds ratios
    write.csv(
      or_data,
      paste0(output_file_prefix, "_odds_ratios.csv"),
      row.names = FALSE
    )
    
    # Return both model and odds ratios
    return(list(model = model, odds_ratios = or_data, error = NULL))
  }, error = function(e) {
    # Handle any other errors in model fitting
    warning("Error in logistic regression: ", e$message)
    
    sink(paste0(output_file_prefix, "_logistic_model_summary.txt"))
    cat("ERROR in logistic regression: ", e$message, "\n")
    sink()
    
    return(list(model = NULL, odds_ratios = NULL, error = e$message))
  })
}

# Build logistic models
ppv_model <- build_logistic_model(
  ppv_vision_analysis,
  "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/ppv"
)

cat_model <- build_logistic_model(
  cat_vision_analysis,
  "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/cat"
)

# -----------------------------------------------------
# 9. Save processed data for further analysis
# -----------------------------------------------------
# Save the merged datasets
saveRDS(ppv_vision_analysis, "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/ppv_vision_analysis.rds")
saveRDS(cat_vision_analysis, "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/cat_vision_analysis.rds")
saveRDS(all_vision_analysis, "3_data_analysis/6_clustering_modeling/vision_analysis/all_metrics/all_vision_analysis.rds")

# Print a summary of the analysis
cat("\n=========================================================\n")
cat("Vision Improvement Analysis Complete\n")
cat("=========================================================\n")

# Print summary of significant parameters
cat("\nPPV - Significant Vision Improvement Parameters (p < 0.05):\n")
if(nrow(ppv_sig_summary) > 0) {
  print(ppv_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

cat("\nCataract - Significant Vision Improvement Parameters (p < 0.05):\n")
if(nrow(cat_sig_summary) > 0) {
  print(cat_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

# Print summary of top membership correlations
cat("\nPPV Membership Correlations with Vision Improvement:\n")
print(ppv_membership_corr)

cat("\nCataract Membership Correlations with Vision Improvement:\n")
print(cat_membership_corr)

# -----------------------------------------------------
# 10. Create combined visualization with multiple timepoints
# -----------------------------------------------------
# Function to create a combined plot showing improvement at both timepoints
create_combined_timepoint_plot <- function(data, surgery_type) {
  # Make sure cluster is a factor for proper handling in ggplot
  data$cluster <- as.factor(data$cluster)
  
  # Convert data to long format for easier plotting
  long_data <- data %>%
    pivot_longer(
      cols = c(vision_improvement_1w, vision_improvement_1m),
      names_to = "timepoint",
      values_to = "improvement"
    ) %>%
    mutate(
      timepoint = factor(
        timepoint,
        levels = c("vision_improvement_1w", "vision_improvement_1m"),
        labels = c("1 Week", "1 Month")
      )
    )
  
  # Create faceted plot
  p <- ggplot(long_data, aes(x = cluster, y = improvement, fill = timepoint)) +
    geom_boxplot(width = 0.7, position = position_dodge(0.8)) +
    # Fix the jitter issue by removing the width parameter
    geom_point(aes(color = timepoint), alpha = 0.4, size = 3,
               position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
    labs(
      title = paste(surgery_type, "- Vision Improvement Over Time by Cluster"),
      x = "Cluster",
      y = "Vision Improvement",
      fill = "Time Point",
      color = "Time Point"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1")
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/vision_analysis/plots/", 
           tolower(surgery_type), "_combined_timepoints.pdf"),
    p,
    width = 10,
    height = 8
  )
  
  return(p)
}

# Create combined timepoint plots
ppv_combined <- create_combined_timepoint_plot(ppv_vision_analysis, "PPV")
cat_combined <- create_combined_timepoint_plot(cat_vision_analysis, "Cataract")
