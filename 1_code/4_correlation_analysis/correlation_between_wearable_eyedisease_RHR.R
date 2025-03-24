library(tidyverse)
library(pROC)
library(ggplot2)
library(lme4)      # For mixed effects models
library(emmeans)   # For marginal means analysis
library(ggpubr)    # For arranging plots
library(tableone)  # For descriptive statistics
library(cosinor)
library(nlme)

setwd(get_project_wd())
rm(list = ls())

##########################
# 1. Load Data
##########################
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")


##########################
# 2. Data Preparation
##########################
# Get unique IDs from heart rate data
heart_rate_ids <- unique(daily_rhr_result$subject_id)

# Filter baseline info to include only participants with heart rate data
baseline_info_filtered <- baseline_info %>% 
  filter(ID %in% heart_rate_ids)
cat("Number of participants in filtered baseline data:", nrow(baseline_info_filtered), "\n")

# Create diabetes and cataract indicator variables
baseline_info_filtered <- baseline_info_filtered %>%
  mutate(
    cataract_2 = case_when(
      cataract == 1 ~ 0,  # No cataract
      cataract %in% c(2, 3, 4) ~ 1, # Has cataract
      TRUE ~ NA_real_
    ),
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  # Has diabetes history
      diabetes_history == 2 ~ 0,  # No diabetes history
      TRUE ~ NA_real_
    )
  )

# Print diabetes status distribution
cat("Number of participants with diabetes:", sum(baseline_info_filtered$dm_2 == 1, na.rm = TRUE), "\n")
cat("Number of participants without diabetes:", sum(baseline_info_filtered$dm_2 == 0, na.rm = TRUE), "\n")

# Create analysis results directory
dir.create("3_data_analysis/4_correlation_analysis/RHR/RHR_diabetes_comparison", recursive = TRUE)
setwd("3_data_analysis/4_correlation_analysis/RHR/RHR_diabetes_comparison")

##########################
# 3. Merge RHR and Diabetes Data
##########################
# Create simplified baseline info with only ID, diabetes status, age, and gender
dm_status <- data.frame(
  subject_id = baseline_info_filtered$ID,
  dm_2 = baseline_info_filtered$dm_2,
  age = baseline_info_filtered$age,
  gender = baseline_info_filtered$gender
)

# Merge RHR and diabetes status data
diabetes_rhr_data <- daily_rhr_result %>%
  left_join(dm_status, by = "subject_id") %>%
  filter(!is.na(dm_2))  # Remove participants with missing diabetes status

# Print information about the final dataset
cat("Final number of participants with complete data:", n_distinct(diabetes_rhr_data$subject_id), "\n")

##########################
# 4. Convert Wide Data to Long Format
##########################
cat("\nConverting data to long format...\n")

# Identify demographic columns
demo_cols <- c("subject_id", "surgery_date", "dm_2", "age", "gender")
actual_demo_cols <- names(diabetes_rhr_data)[names(diabetes_rhr_data) %in% demo_cols]

# Get all RHR columns (excluding demographic columns)
rhr_cols <- names(diabetes_rhr_data)[!names(diabetes_rhr_data) %in% actual_demo_cols]

# Create long format data
long_data <- diabetes_rhr_data %>%
  pivot_longer(
    cols = all_of(rhr_cols),
    names_to = "variable",
    values_to = "value"
  ) %>%
  # Extract day, statistic type, and RHR measurement type
  mutate(
    day = as.numeric(str_extract(variable, "-?\\d+")),
    stat_type = str_extract(variable, "(mean|min|max|median|sd|cv|iqr|skew|kurt)"),
    rhr_type = str_extract(variable, "rhr_\\d+")
  ) %>%
  # Remove missing values and filter for days -10 to 30
  filter(!is.na(value), day >= -7, day <= 30)

# Print long data summary
cat("Long data summary:\n")
cat("Total observations:", nrow(long_data), "\n")
cat("Unique days:", n_distinct(long_data$day), "\n")
cat("Statistical metrics:", paste(unique(long_data$stat_type), collapse = ", "), "\n")
cat("RHR types:", paste(unique(long_data$rhr_type), collapse = ", "), "\n")

##########################
# 5. Compare RHR Metrics by Diabetes Status
##########################

# Group data by day, stat_type, and rhr_type
comparison_results <- long_data %>%
  group_by(day, stat_type, rhr_type) %>%
  # Perform t-test for each group
  summarize(
    n_diabetes = sum(dm_2 == 1, na.rm = TRUE),
    n_non_diabetes = sum(dm_2 == 0, na.rm = TRUE),
    mean_diabetes = mean(value[dm_2 == 1], na.rm = TRUE),
    mean_non_diabetes = mean(value[dm_2 == 0], na.rm = TRUE),
    sd_diabetes = sd(value[dm_2 == 1], na.rm = TRUE),
    sd_non_diabetes = sd(value[dm_2 == 0], na.rm = TRUE),
    # Perform t-test only if there are at least 2 observations in each group
    p_value = ifelse(
      sum(dm_2 == 1, na.rm = TRUE) >= 2 & sum(dm_2 == 0, na.rm = TRUE) >= 2,
      t.test(value ~ dm_2)$p.value,
      NA
    ),
    # Calculate Cohen's d effect size
    cohens_d = ifelse(
      sum(dm_2 == 1, na.rm = TRUE) >= 2 & sum(dm_2 == 0, na.rm = TRUE) >= 2,
      (mean(value[dm_2 == 1], na.rm = TRUE) - mean(value[dm_2 == 0], na.rm = TRUE)) / 
        sqrt((sd(value[dm_2 == 1], na.rm = TRUE)^2 + sd(value[dm_2 == 0], na.rm = TRUE)^2) / 2),
      NA
    ),
    # Add significance indicator
    significant = p_value < 0.05,
    .groups = "drop"
  )

# Apply FDR correction for multiple testing
comparison_results <- comparison_results %>%
  group_by(stat_type, rhr_type) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "fdr"),
    significant_adjusted = p_adjusted < 0.05
  ) %>%
  ungroup()

# Create summary table of significant findings
significant_findings <- comparison_results %>%
  filter(significant_adjusted == TRUE) %>%
  arrange(stat_type, rhr_type, day)

cat("Found", nrow(significant_findings), "statistically significant differences (FDR-adjusted)\n")

# Save comparison results
write.csv(comparison_results, "diabetes_rhr_comparison_results.csv", row.names = FALSE)
write.csv(significant_findings, "diabetes_rhr_significant_findings.csv", row.names = FALSE)

##########################
# 6. Visualization: Plot RHR Mean by Day
##########################

# Mean RHR by day visualization
for(rhr_type_filter in c("rhr_1", "rhr_50")) {
  cat("Creating time series plot for", rhr_type_filter, "...\n")
  
  # Filter data
  plot_data <- long_data %>%
    filter(
      stat_type == "mean",
      rhr_type == rhr_type_filter
    )
  
  # Calculate averages by day and diabetes status
  avg_data <- plot_data %>%
    group_by(day, dm_2) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      dm_label = ifelse(dm_2 == 1, "DM", "No DM"),
      lower_ci = mean_value - 1.96 * se_value,
      upper_ci = mean_value + 1.96 * se_value
    )
  
  # Create plot
  p <- ggplot(avg_data, aes(x = day, y = mean_value, color = dm_label, group = dm_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = dm_label), 
                alpha = 0.2, color = NA) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("DM" = "#e6550d", "No DM" = "#3182bd")) +
    scale_fill_manual(values = c("DM" = "#e6550d", "No DM" = "#3182bd")) +
    labs(
      title = paste("Mean daily RHR by diabetes status"),
      subtitle = paste("RHR type:", rhr_type_filter),
      x = "Days relative to surgery",
      y = "RHR mean",
      color = "DM status",
      fill = "DM status"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  # Display and save plot
  print(p)
  file_name <- paste0("daily_rhr_mean_", rhr_type_filter, "_by_diabetes.png")
  ggsave(file_name, p, width = 10, height = 6, dpi = 300)
  cat("Plot saved as:", file_name, "\n")
}

##########################
# 7. Visualization: Statistical Significance
##########################
# Create p-value plots for different RHR metrics
for(stat_type_filter in c("mean", "sd", "cv")) {
  for(rhr_type_filter in c("rhr_1", "rhr_50")) {
    cat("Creating p-value plot for", stat_type_filter, rhr_type_filter, "...\n")
    
    # Filter data
    plot_data <- comparison_results %>%
      filter(
        stat_type == stat_type_filter,
        rhr_type == rhr_type_filter,
        !is.na(p_value)  # Remove entries with NA p-values
      )
    
    # Check if we have enough data
    if(nrow(plot_data) == 0) {
      cat("No data available for the specified filters.\n")
      next
    }
    
    # Create p-value plot
    p <- ggplot(plot_data, aes(x = day, y = -log10(p_value))) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
      geom_point(aes(color = significant_adjusted), size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgray")) +
      labs(
        title = paste("Statistical significance of diabetes effect on RHR", stat_type_filter),
        subtitle = paste("RHR type:", rhr_type_filter),
        x = "Days relative to surgery",
        y = "-log10(p-value)",
        color = "Significant (FDR adj.)"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      )
    
    # Display and save plot
    print(p)
    file_name <- paste0("pvalues_", stat_type_filter, "_", rhr_type_filter, "_diabetes.png")
    ggsave(file_name, p, width = 10, height = 6, dpi = 300)
    cat("Plot saved as:", file_name, "\n")
  }
}

##########################
# 8. Visualization: Period Comparison
##########################

# Create boxplots comparing periods for different RHR metrics
for(stat_type_filter in c("mean", "sd", "cv")) {
  for(rhr_type_filter in c("rhr_1", "rhr_50")) {
    cat("Creating period comparison plot for", stat_type_filter, rhr_type_filter, "...\n")
    
    # Filter data
    plot_data <- long_data %>%
      filter(
        stat_type == stat_type_filter,
        rhr_type == rhr_type_filter
      ) %>%
      mutate(
        period = case_when(
          day < 0 ~ "Pre-surgery",
          day >= 0 & day <= 30 ~ "Post-surgery (0-30d)"
        )
      )
    
    # Check if we have enough data
    if(nrow(plot_data) == 0) {
      cat("No data available for the specified filters.\n")
      next
    }
    
    # Create boxplot
    p <- ggplot(plot_data, 
                aes(x = period, y = value, fill = factor(dm_2))) +
      geom_boxplot(outlier.size = 1, width = 0.7, position = position_dodge(0.8)) +
      scale_fill_manual(
        values = c("0" = "#3182bd", "1" = "#e6550d"),
        labels = c("0" = "No DM", "1" = "DM")
      ) +
      labs(
        title = paste("RHR", stat_type_filter, "by surgical period and diabetes status"),
        x = "Surgical period",
        y = paste("RHR", stat_type_filter),
        fill = "DM status"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Add p-values
    p <- p + stat_compare_means(
      aes(group = factor(dm_2)),
      label = "p.format",
      method = "t.test",
      label.y.npc = 0.95
    )
    
    # Display and save plot
    print(p)
    file_name <- paste0("period_comparison_", stat_type_filter, "_", rhr_type_filter, "_diabetes.png")
    ggsave(file_name, p, width = 10, height = 6, dpi = 300)
    cat("Plot saved as:", file_name, "\n")
  }
}

##########################
# 9. Longitudinal Analysis
##########################
# Analyze different metrics
metrics <- list(
  list(stat = "mean", rhr = "rhr_1"),
  list(stat = "mean", rhr = "rhr_50"),
  list(stat = "sd", rhr = "rhr_1"),
  list(stat = "cv", rhr = "rhr_1")
)

# Run analysis for each metric
longitudinal_results <- list()
for(metric in metrics) {
  cat("Analyzing", metric$stat, metric$rhr, "...\n")
  
  # Filter data
  filtered_data <- long_data %>%
    filter(stat_type == metric$stat, rhr_type == metric$rhr) %>%
    mutate(
      time_period = case_when(
        day < 0 ~ "Pre-surgery",
        day >= 0 & day <= 30 ~ "Post-surgery",
        TRUE ~ NA_character_
      )
    )
  
  # Build model formula
  formula_str <- "value ~ dm_2 * day + age + gender + (1|subject_id)"
  
  # Try running mixed effects model
  tryCatch({
    model <- lmer(as.formula(formula_str), data = filtered_data)
    
    # Display model results summary
    cat("Fixed effects:\n")
    print(summary(model)$coefficients)
    
    # Get estimated marginal means for diabetes status
    emm <- emmeans(model, ~ dm_2)
    cat("Marginal means:\n")
    print(emm)
    
    contrasts <- contrast(emm, "pairwise")
    cat("Contrasts:\n")
    print(contrasts)
    
    # Save results
    longitudinal_results[[paste(metric$stat, metric$rhr, sep="_")]] <- list(
      model = model,
      summary = summary(model),
      emm = emm,
      contrasts = contrasts
    )
    
  }, error = function(e) {
    # If model fails, try simpler model
    cat("Full model failed, error:", e$message, "\n")
    cat("Trying simpler model...\n")
    
    tryCatch({
      # Simpler model without interaction
      model <- lmer(value ~ dm_2 + day + age + (1|subject_id), data = filtered_data)
      
      # Display model results summary
      cat("Simplified model fixed effects:\n")
      print(summary(model)$coefficients)
      
      # Get estimated marginal means for diabetes status
      emm <- emmeans(model, ~ dm_2)
      cat("Marginal means:\n")
      print(emm)
      
      contrasts <- contrast(emm, "pairwise")
      cat("Contrasts:\n")
      print(contrasts)
      
      # Save results
      longitudinal_results[[paste(metric$stat, metric$rhr, sep="_")]] <- list(
        model = model,
        summary = summary(model),
        emm = emm,
        contrasts = contrasts,
        note = "Used simplified model due to convergence issues with full model"
      )
      
    }, error = function(e2) {
      cat("Simplified model also failed, error:", e2$message, "\n")
    })
  })
}

##########################
# 10. Save All Results
##########################

# Save processed data
save(diabetes_rhr_data, long_data, comparison_results, longitudinal_results, significant_findings,
     file = "diabetes_rhr_analysis.rda")

##########################
# 11. Create Summary Report
##########################
cat("\nCreating summary report...\n")

# Analysis parameters summary
cat("Analysis Summary Report\n")
cat("====================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d"), "\n")
cat("Total participants analyzed:", n_distinct(diabetes_rhr_data$subject_id), "\n")
cat("Participants with diabetes:", sum(diabetes_rhr_data$dm_2 == 1, na.rm = TRUE), "\n")
cat("Participants without diabetes:", sum(diabetes_rhr_data$dm_2 == 0, na.rm = TRUE), "\n")
cat("Days analyzed:", min(long_data$day), "to", max(long_data$day), "(relative to surgery)\n")
cat("RHR metrics analyzed:", paste(unique(long_data$stat_type), collapse = ", "), "\n")
cat("RHR types analyzed:", paste(unique(long_data$rhr_type), collapse = ", "), "\n\n")

# Significant findings summary
cat("Significant Findings Summary\n")
cat("--------------------------\n\n")
cat("Total significant findings (FDR-adjusted):", nrow(significant_findings), "\n\n")

if(nrow(significant_findings) > 0) {
  # Group significant findings by metric
  sig_by_metric <- significant_findings %>%
    group_by(stat_type, rhr_type) %>%
    summarize(
      n_findings = n(),
      mean_effect_size = mean(abs(cohens_d), na.rm = TRUE),
      days = paste(sort(day), collapse = ", "),
      .groups = "drop"
    )
  
  # Print summary by metric
  for(i in 1:nrow(sig_by_metric)) {
    cat(sig_by_metric$stat_type[i], " (", sig_by_metric$rhr_type[i], "): ", 
        sig_by_metric$n_findings[i], " significant days\n", sep = "")
    cat("  Mean effect size (Cohen's d): ", round(sig_by_metric$mean_effect_size[i], 2), "\n")
    cat("  Days with significant differences: ", sig_by_metric$days[i], "\n\n")
  }
  
  # Highlight strongest findings
  top_findings <- significant_findings %>%
    arrange(desc(abs(cohens_d))) %>%
    head(5)
  
  cat("Top 5 Strongest Effects:\n")
  for(i in 1:nrow(top_findings)) {
    cat(i, ". ", top_findings$stat_type[i], " (", top_findings$rhr_type[i], "), Day ", 
        top_findings$day[i], ", Cohen's d = ", round(top_findings$cohens_d[i], 2), 
        ", p = ", format(top_findings$p_adjusted[i], digits = 3), "\n", sep = "")
  }
}

# Save summary report
sink("diabetes_rhr_analysis_summary.txt")
cat("Analysis Summary Report\n")
cat("====================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d"), "\n")
cat("Total participants analyzed:", n_distinct(diabetes_rhr_data$subject_id), "\n")
cat("Participants with diabetes:", sum(diabetes_rhr_data$dm_2 == 1, na.rm = TRUE), "\n")
cat("Participants without diabetes:", sum(diabetes_rhr_data$dm_2 == 0, na.rm = TRUE), "\n")
cat("Days analyzed:", min(long_data$day), "to", max(long_data$day), "(relative to surgery)\n")
cat("RHR metrics analyzed:", paste(unique(long_data$stat_type), collapse = ", "), "\n")
cat("RHR types analyzed:", paste(unique(long_data$rhr_type), collapse = ", "), "\n\n")

cat("Significant Findings Summary\n")
cat("--------------------------\n\n")
cat("Total significant findings (FDR-adjusted):", nrow(significant_findings), "\n\n")

if(nrow(significant_findings) > 0) {
  # Group significant findings by metric
  sig_by_metric <- significant_findings %>%
    group_by(stat_type, rhr_type) %>%
    summarize(
      n_findings = n(),
      mean_effect_size = mean(abs(cohens_d), na.rm = TRUE),
      days = paste(sort(day), collapse = ", "),
      .groups = "drop"
    )
  
  # Print summary by metric
  for(i in 1:nrow(sig_by_metric)) {
    cat(sig_by_metric$stat_type[i], " (", sig_by_metric$rhr_type[i], "): ", 
        sig_by_metric$n_findings[i], " significant days\n", sep = "")
    cat("  Mean effect size (Cohen's d): ", round(sig_by_metric$mean_effect_size[i], 2), "\n")
    cat("  Days with significant differences: ", sig_by_metric$days[i], "\n\n")
  }
  
  # Highlight strongest findings
  top_findings <- significant_findings %>%
    arrange(desc(abs(cohens_d))) %>%
    head(5)
  
  cat("Top 5 Strongest Effects:\n")
  for(i in 1:nrow(top_findings)) {
    cat(i, ". ", top_findings$stat_type[i], " (", top_findings$rhr_type[i], "), Day ", 
        top_findings$day[i], ", Cohen's d = ", round(top_findings$cohens_d[i], 2), 
        ", p = ", format(top_findings$p_adjusted[i], digits = 3), "\n", sep = "")
  }
}
sink()



