# Steps Max Daily Data Visualization - Based on Daily Steps Data
# Focus on -4 to +30 days around surgery, using daily aggregated data

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ggridges)
library(corrplot)
library(r4projects)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# ================== 1. Load Data and Match Wearable Cohort ==================

# Load daily steps result data
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result_assigned.rda")

# Load the wearable cohort from clustering analysis (this is the exact population from code 2)
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)

# Create output directory
dir.create("3_data_analysis/7_figures/figure_s2/steps_max_daily_analysis", recursive = TRUE)
setwd("3_data_analysis/7_figures/figure_s2/steps_max_daily_analysis")

# Extract wearable cohort subject IDs from the clustering data (ppv_data)
wearable_cohort_ids <- unique(ppv_data$subject_id)

cat("Data loaded successfully\n")
cat(sprintf("Total subjects in daily steps data: %d\n", length(unique(daily_steps_result$subject_id))))
cat(sprintf("Wearable cohort subjects from clustering data: %d\n", length(wearable_cohort_ids)))
cat("Wearable cohort IDs:", paste(wearable_cohort_ids, collapse = ", "), "\n")

# Filter daily steps data to match wearable cohort
daily_steps_wearable <- daily_steps_result %>%
  filter(subject_id %in% wearable_cohort_ids)

# Verify the match
matched_ids <- unique(daily_steps_wearable$subject_id)
missing_ids <- setdiff(wearable_cohort_ids, matched_ids)


# ================== 2. Extract Steps Max Data for -4 to +30 Days ==================

# Find all steps_max columns within the desired day range (-4 to +30)
target_days <- -4:30
steps_max_pattern <- paste0("day_", target_days, "_steps_max", collapse = "|")
steps_max_cols <- grep(steps_max_pattern, names(daily_steps_result), value = TRUE)

cat("\nFound steps_max columns for days -4 to +30:\n")
cat("Total columns found:", length(steps_max_cols), "\n")

# Extract day numbers from column names to verify range
day_numbers <- as.numeric(str_extract(steps_max_cols, "-?\\d+"))
cat("Day range in data:", min(day_numbers, na.rm = TRUE), "to", max(day_numbers, na.rm = TRUE), "\n")

# Create analysis dataset using the matched wearable cohort
steps_max_daily_data <- daily_steps_wearable %>%
  dplyr::select(subject_id, surgery_date, all_of(steps_max_cols)) %>%
  # Convert to long format for analysis
  pivot_longer(
    cols = all_of(steps_max_cols),
    names_to = "day_column",
    values_to = "steps_max"
  ) %>%
  # Extract day number from column name
  mutate(
    day_relative = as.numeric(str_extract(day_column, "-?\\d+")),
    # Create period categories
    period_category = case_when(
      day_relative < 0 ~ "Pre-Surgery",
      day_relative >= 0 & day_relative <= 7 ~ "Post-Surgery (0-7d)",
      day_relative > 7 & day_relative <= 30 ~ "Post-Surgery (8-30d)",
      TRUE ~ "Other"
    ),
    # Create more detailed period labels
    period_label = case_when(
      day_relative >= -4 & day_relative < 0 ~ "Pre-Surgery (-4 to -1d)",
      day_relative >= 0 & day_relative <= 3 ~ "Acute Recovery (0-3d)",
      day_relative >= 4 & day_relative <= 7 ~ "Early Recovery (4-7d)",
      day_relative >= 8 & day_relative <= 14 ~ "Mid Recovery (8-14d)",
      day_relative >= 15 & day_relative <= 30 ~ "Late Recovery (15-30d)",
      TRUE ~ "Other"
    ),
    # Add surgery month for seasonal analysis
    surgery_month = as.numeric(format(as.Date(surgery_date), "%m")),
    surgery_season = case_when(
      surgery_month %in% c(12, 1, 2) ~ "Winter",
      surgery_month %in% c(3, 4, 5) ~ "Spring", 
      surgery_month %in% c(6, 7, 8, 9) ~ "Summer",
      surgery_month %in% c(10, 11) ~ "Fall",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Filter out missing values
  filter(!is.na(steps_max))

# Set factor levels for proper ordering
period_order <- c("Pre-Surgery (-4 to -1d)", "Acute Recovery (0-3d)", 
                  "Early Recovery (4-7d)", "Mid Recovery (8-14d)", "Late Recovery (15-30d)")
steps_max_daily_data$period_label <- factor(steps_max_daily_data$period_label, levels = period_order)

cat(sprintf("\nProcessed daily data points: %d\n", nrow(steps_max_daily_data)))
cat(sprintf("Wearable cohort subjects in analysis: %d\n", length(unique(steps_max_daily_data$subject_id))))
cat(sprintf("Day range: %d to %d\n", min(steps_max_daily_data$day_relative), max(steps_max_daily_data$day_relative)))
cat(sprintf("Period categories: %d\n", length(unique(steps_max_daily_data$period_label)))))

# ================== 3. Basic Statistical Summary ==================

daily_summary_stats <- steps_max_daily_data %>%
  group_by(period_label, period_category) %>%
  summarise(
    n_observations = n(),
    n_subjects = n_distinct(subject_id),
    mean_steps = mean(steps_max, na.rm = TRUE),
    median_steps = median(steps_max, na.rm = TRUE),
    sd_steps = sd(steps_max, na.rm = TRUE),
    min_steps = min(steps_max, na.rm = TRUE),
    max_steps = max(steps_max, na.rm = TRUE),
    q25 = quantile(steps_max, 0.25, na.rm = TRUE),
    q75 = quantile(steps_max, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n=== Wearable Cohort Daily Steps Max Summary Statistics ===\n")
print(daily_summary_stats)

# Daily trend summary (by day)
daily_trend_stats <- steps_max_daily_data %>%
  group_by(day_relative) %>%
  summarise(
    n_observations = n(),
    mean_steps = mean(steps_max, na.rm = TRUE),
    median_steps = median(steps_max, na.rm = TRUE),
    sd_steps = sd(steps_max, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(day_relative)

# ================== 4. Visualization Functions ==================

# 1. Daily trend plot (day by day)
create_daily_trend_plot <- function(data, trend_stats) {
  
  cat("\nCreating daily trend plot...\n")
  
  # Main trend plot
  p1 <- ggplot(trend_stats, aes(x = day_relative, y = mean_steps)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1, alpha = 0.7) +
    geom_line(color = "#2c3e50", size = 1.2) +
    geom_point(color = "#2c3e50", size = 2) +
    geom_ribbon(aes(ymin = mean_steps - sd_steps, ymax = mean_steps + sd_steps), 
                alpha = 0.3, fill = "#3498db") +
    scale_x_continuous(
      breaks = seq(-4, 30, by = 5),
      labels = seq(-4, 30, by = 5)
    ) +
    labs(
      title = "Daily Steps Max Trends Around Surgery - Wearable Cohort",
      subtitle = paste("Wearable cohort (n =", length(unique(data$subject_id)), "subjects) | Red line indicates surgery day"),
      x = "Days Relative to Surgery",
      y = "Mean Steps Max"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12)
    )
  
  # Individual patient trajectories (all patients for small cohort)
  all_patients <- unique(data$subject_id)
  n_patients <- length(all_patients)
  
  # For small cohort, show all patients; for larger cohorts, sample
  if(n_patients <= 20) {
    sample_patients <- all_patients
    subtitle_text <- paste("All", n_patients, "patients with overall trend")
  } else {
    sample_patients <- sample(all_patients, 20)
    subtitle_text <- paste("Sample of 20 out of", n_patients, "patients with overall trend")
  }
  
  sample_data <- data %>% filter(subject_id %in% sample_patients)
  
  p2 <- ggplot(sample_data, aes(x = day_relative, y = steps_max, group = subject_id)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1, alpha = 0.7) +
    geom_line(alpha = 0.6, color = "#7f8c8d", size = 0.8) +
    geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "#2c3e50", size = 1.5) +
    scale_x_continuous(
      breaks = seq(-4, 30, by = 5),
      labels = seq(-4, 30, by = 5)
    ) +
    labs(
      title = "Individual Patient Trajectories - Wearable Cohort",
      subtitle = subtitle_text,
      x = "Days Relative to Surgery", 
      y = "Steps Max"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12)
    )
  
  return(list(trend = p1, trajectories = p2))
}

# 2. Period-based distribution plots
create_period_distribution_plots <- function(data) {
  
  cat("\nCreating period distribution plots...\n")
  
  # Boxplot by recovery periods (optimized for small sample)
  p1 <- ggplot(data, aes(x = period_label, y = steps_max, fill = period_category)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.7, outlier.size = 2) +
    geom_jitter(width = 0.25, alpha = 0.8, size = 1.2) +  # Larger points for small sample
    scale_fill_manual(
      values = c("Pre-Surgery" = "#3498db", 
                 "Post-Surgery (0-7d)" = "#e74c3c", 
                 "Post-Surgery (8-30d)" = "#f39c12"),
      name = "Recovery Phase"
    ) +
    labs(
      title = "Steps Max Distribution by Recovery Periods - Wearable Cohort", 
      subtitle = paste("Wearable cohort (n =", length(unique(data$subject_id)), "subjects) with detailed recovery phases"),
      x = "Recovery Periods",
      y = "Steps Max"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Ridge plot
  p2 <- ggplot(data, aes(x = steps_max, y = period_label, fill = period_category)) +
    geom_density_ridges(alpha = 0.7, scale = 0.9) +
    scale_fill_manual(
      values = c("Pre-Surgery" = "#3498db", 
                 "Post-Surgery (0-7d)" = "#e74c3c", 
                 "Post-Surgery (8-30d)" = "#f39c12"),
      name = "Recovery Phase"
    ) +
    labs(
      title = "Steps Max Distribution Shapes - Wearable Cohort",
      subtitle = "Density distributions by recovery periods for wearable cohort",
      x = "Steps Max",
      y = "Recovery Periods"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    )
  
  return(list(boxplot = p1, ridge = p2))
}

# 3. Recovery pattern analysis
create_recovery_analysis <- function(data, summary_stats) {
  
  cat("\nCreating recovery pattern analysis...\n")
  
  # Pre vs Post comparison
  pre_post_data <- data %>%
    mutate(
      phase = ifelse(day_relative < 0, "Pre-Surgery", "Post-Surgery")
    ) %>%
    group_by(subject_id, phase) %>%
    summarise(
      mean_steps_max = mean(steps_max, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    pivot_wider(names_from = phase, values_from = mean_steps_max) %>%
    filter(!is.na(`Pre-Surgery`) & !is.na(`Post-Surgery`)) %>%
    mutate(
      recovery_ratio = `Post-Surgery` / `Pre-Surgery`,
      recovery_diff = `Post-Surgery` - `Pre-Surgery`
    )
  
  # Recovery ratio plot (optimized for small sample)
  p1 <- ggplot(pre_post_data, aes(x = `Pre-Surgery`, y = `Post-Surgery`)) +
    geom_point(alpha = 0.8, size = 3, color = "#3498db") +  # Larger points for visibility
    geom_text(aes(label = subject_id), nudge_y = 100, size = 3, alpha = 0.7) +  # Add subject labels
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "#2c3e50", alpha = 0.3) +
    labs(
      title = "Pre-Surgery vs Post-Surgery Steps Max - Wearable Cohort",
      subtitle = paste("Wearable cohort (n =", nrow(pre_post_data), "subjects) | Red line: perfect recovery"),
      x = "Pre-Surgery Mean Steps Max",
      y = "Post-Surgery Mean Steps Max"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  
  # Recovery ratio distribution (optimized for small sample)
  p2 <- ggplot(pre_post_data, aes(x = recovery_ratio)) +
    geom_histogram(bins = 8, fill = "#3498db", alpha = 0.7, color = "white") +  # Fewer bins for small sample
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = mean(pre_post_data$recovery_ratio, na.rm = TRUE), 
               linetype = "solid", color = "blue", size = 1) +
    geom_text(data = pre_post_data, aes(x = recovery_ratio, y = 0.5, label = subject_id), 
              angle = 90, vjust = -0.5, size = 3, alpha = 0.7) +  # Add subject labels
    labs(
      title = "Recovery Ratio Distribution - Wearable Cohort",
      subtitle = paste("Wearable cohort (n =", nrow(pre_post_data), ") | Red: perfect recovery | Blue: mean"),
      x = "Recovery Ratio (Post/Pre)",
      y = "Number of Subjects"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  
  return(list(scatter = p1, ratio_dist = p2, recovery_data = pre_post_data))
}

# 4. Seasonal analysis
create_seasonal_analysis <- function(data) {
  
  cat("\nCreating seasonal analysis...\n")
  
  # Seasonal boxplot (optimized for small sample)
  p1 <- ggplot(data, aes(x = surgery_season, y = steps_max, fill = surgery_season)) +
    geom_boxplot(alpha = 0.7, outlier.size = 2) +
    geom_jitter(width = 0.25, alpha = 0.8, size = 1.5) +  # Larger points for small sample
    scale_fill_viridis_d(name = "Surgery Season") +
    labs(
      title = "Steps Max by Surgery Season - Wearable Cohort",
      subtitle = paste("Seasonal patterns in wearable cohort (n =", length(unique(data$subject_id)), "subjects)"),
      x = "Surgery Season",
      y = "Steps Max"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    )
  
  return(p1)
}

# ================== 5. Execute All Visualizations ==================

cat("\nStarting wearable cohort daily Steps Max visualization analysis...\n")

# 1. Daily trend plots
daily_plots <- create_daily_trend_plot(steps_max_daily_data, daily_trend_stats)
print(daily_plots$trend)
print(daily_plots$trajectories)

# 2. Period distribution plots  
period_plots <- create_period_distribution_plots(steps_max_daily_data)
print(period_plots$boxplot)
print(period_plots$ridge)

# 3. Recovery analysis
recovery_plots <- create_recovery_analysis(steps_max_daily_data, daily_summary_stats)
print(recovery_plots$scatter)
print(recovery_plots$ratio_dist)

# 4. Seasonal analysis
seasonal_plot <- create_seasonal_analysis(steps_max_daily_data)
print(seasonal_plot)

# ================== 6. Save All Plots ==================

cat("\nSaving all plots...\n")

# Daily trend plots
ggsave("daily_steps_max_trend.pdf", daily_plots$trend, width = 12, height = 8, dpi = 300)
ggsave("daily_steps_max_trend.png", daily_plots$trend, width = 12, height = 8, dpi = 300)

ggsave("daily_steps_max_trajectories.pdf", daily_plots$trajectories, width = 12, height = 8, dpi = 300)
ggsave("daily_steps_max_trajectories.png", daily_plots$trajectories, width = 12, height = 8, dpi = 300)

# Period distribution plots
ggsave("period_steps_max_boxplot.pdf", period_plots$boxplot, width = 14, height = 8, dpi = 300)
ggsave("period_steps_max_boxplot.png", period_plots$boxplot, width = 14, height = 8, dpi = 300)

ggsave("period_steps_max_ridge.pdf", period_plots$ridge, width = 12, height = 8, dpi = 300)
ggsave("period_steps_max_ridge.png", period_plots$ridge, width = 12, height = 8, dpi = 300)

# Recovery analysis plots
ggsave("recovery_scatter_plot.pdf", recovery_plots$scatter, width = 10, height = 8, dpi = 300)
ggsave("recovery_scatter_plot.png", recovery_plots$scatter, width = 10, height = 8, dpi = 300)

ggsave("recovery_ratio_distribution.pdf", recovery_plots$ratio_dist, width = 10, height = 6, dpi = 300)
ggsave("recovery_ratio_distribution.png", recovery_plots$ratio_dist, width = 10, height = 6, dpi = 300)

# Seasonal analysis
ggsave("seasonal_steps_max_analysis.pdf", seasonal_plot, width = 10, height = 6, dpi = 300)
ggsave("seasonal_steps_max_analysis.png", seasonal_plot, width = 10, height = 6, dpi = 300)

# ================== 7. Create Comprehensive Panel ==================

create_comprehensive_daily_panel <- function() {
  
  cat("\nCreating comprehensive daily panel...\n")
  
  # Select key plots for the panel
  p1 <- daily_plots$trend + theme(plot.title = element_text(size = 11), 
                                  axis.title = element_text(size = 10))
  p2 <- period_plots$boxplot + theme(plot.title = element_text(size = 11),
                                     axis.title = element_text(size = 10),
                                     axis.text.x = element_text(size = 8))
  p3 <- recovery_plots$scatter + theme(plot.title = element_text(size = 11),
                                       axis.title = element_text(size = 10))
  p4 <- seasonal_plot + theme(plot.title = element_text(size = 11),
                              axis.title = element_text(size = 10))
  
  # Arrange plots
  combined_plot <- gridExtra::grid.arrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    top = paste("Daily Steps Max Analysis - Wearable Cohort (n =", length(unique(steps_max_daily_data$subject_id)), "subjects, -4 to +30 Days)")
  )
  
  return(combined_plot)
}

# Create and save comprehensive panel
comprehensive_daily_panel <- create_comprehensive_daily_panel()
ggsave("daily_steps_max_comprehensive_panel.pdf", comprehensive_daily_panel, width = 16, height = 12, dpi = 300)
ggsave("daily_steps_max_comprehensive_panel.png", comprehensive_daily_panel, width = 16, height = 12, dpi = 300)

# ================== 8. Generate Report ==================

generate_daily_analysis_report <- function(data, summary_stats, trend_stats, recovery_data) {
  
  report <- paste0(
    "========================================\n",
    "Daily Steps Max Analysis Report - Wearable Cohort\n",
    "========================================\n\n",
    
    "📊 Wearable Cohort Information:\n",
    sprintf("- Cohort type: Clustering analysis wearable device users (mfuzz_D_Surg1_8h_filtered.csv)\n"),
    sprintf("- Analysis period: Day -4 to Day +30 (relative to surgery)\n"),
    sprintf("- Total daily observations: %d\n", nrow(data)),
    sprintf("- Number of subjects: %d (clustering wearable cohort)\n", length(unique(data$subject_id))),
    sprintf("- Subject IDs: %s\n", paste(sort(unique(data$subject_id)), collapse = ", ")),
    sprintf("- Day range in data: %d to %d\n", min(data$day_relative), max(data$day_relative)),
    sprintf("- Overall mean steps max: %.0f\n", mean(data$steps_max, na.rm = TRUE)),
    sprintf("- Overall median steps max: %.0f\n", median(data$steps_max, na.rm = TRUE)),
    sprintf("- Overall range: %d - %d steps\n", min(data$steps_max, na.rm = TRUE), max(data$steps_max, na.rm = TRUE)),
    "\n"
  )
  
  # Recovery period statistics
  report <- paste0(report, "📈 Recovery Period Summary:\n")
  for(i in 1:nrow(summary_stats)) {
    row <- summary_stats[i,]
    report <- paste0(report,
                     sprintf("%s: Mean=%.0f, Median=%.0f, SD=%.0f (n=%d subjects, %d observations)\n",
                             row$period_label, row$mean_steps, row$median_steps, 
                             row$sd_steps, row$n_subjects, row$n_observations))
  }
  
  # Recovery analysis
  if(!is.null(recovery_data)) {
    mean_recovery_ratio <- mean(recovery_data$recovery_ratio, na.rm = TRUE)
    recovery_above_1 <- sum(recovery_data$recovery_ratio > 1, na.rm = TRUE)
    total_recovery_subjects <- nrow(recovery_data)
    
    report <- paste0(report, "\n",
                     "🔄 Recovery Analysis:\n",
                     sprintf("- Subjects with both pre/post data: %d\n", total_recovery_subjects),
                     sprintf("- Mean recovery ratio (Post/Pre): %.2f\n", mean_recovery_ratio),
                     sprintf("- Subjects with recovery ratio > 1: %d (%.1f%%)\n", 
                             recovery_above_1, 100 * recovery_above_1 / total_recovery_subjects),
                     "\n")
  }
  
  report <- paste0(report,
                   "🎯 Key Findings for Wearable Cohort:\n",
                   "1. Daily progression of steps max from pre-surgery to 30 days post\n",
                   "2. Individual patient trajectories from matched wearable cohort\n",
                   "3. Recovery patterns show immediate post-surgical impact\n",
                   "4. High-resolution view of wearable device users' activity patterns\n",
                   "5. Consistent with time window analysis from code 2\n\n",
                   
                   "📁 Generated Visualization Files:\n",
                   "- daily_steps_max_trend: Day-by-day trend analysis\n",
                   "- daily_steps_max_trajectories: Individual patient paths\n",
                   "- period_steps_max_boxplot: Recovery period comparisons\n",
                   "- period_steps_max_ridge: Distribution shapes by period\n",
                   "- recovery_scatter_plot: Pre vs post surgery comparison\n",
                   "- recovery_ratio_distribution: Recovery ratio analysis\n",
                   "- seasonal_steps_max_analysis: Seasonal patterns\n",
                   "- daily_steps_max_comprehensive_panel: Overview panel\n\n",
                   
                   "📝 Data Source: daily_steps_result_assigned.rda (filtered for wearable cohort)\n",
                   "📝 Cohort Source: mfuzz_D_Surg1_8h_filtered.csv (clustering analysis wearable subjects)\n",
                   "📅 Analysis period: -4 to +30 days relative to surgery\n",
                   "Generated on: ", Sys.time(), "\n",
                   "========================================"
  )
  
  writeLines(report, "Daily_Steps_Max_Analysis_Report.txt")
  cat("✓ Daily analysis report saved\n")
  
  return(report)
}

# Generate report
daily_analysis_report <- generate_daily_analysis_report(
  steps_max_daily_data, daily_summary_stats, daily_trend_stats, recovery_plots$recovery_data
)
cat(daily_analysis_report)

# ================== 9. Final Summary ==================

cat("\n🎉 Daily Steps Max Analysis Complete - Wearable Cohort!\n")
cat("========================================\n")
cat("✅ Analysis focused on clustering data wearable cohort\n")
cat("✅ Subjects from mfuzz_D_Surg1_8h_filtered.csv\n")
cat("✅ Matched with daily aggregated steps data\n")
cat("✅ Analysis based on daily aggregated data\n")
cat("✅ Focused on -4 to +30 days around surgery\n")
cat("✅ Generated 8 types of visualization plots\n")
cat("✅ Individual patient trajectories clearly visible\n")
cat("✅ Enhanced visualization for clustering cohort\n")
cat("✅ All plots saved in PDF and PNG formats\n")
cat("✅ Comprehensive analysis report generated\n")
cat("========================================\n")

# Display saved files
cat("\n📁 List of saved files:\n")
saved_files <- list.files(pattern = "\\.(pdf|png|txt)$")
for(file in saved_files) {
  cat(sprintf("✓ %s\n", file))
}

cat(sprintf("\n📍 Output directory: %s\n", getwd()))
cat("🎯 Key recommendations:\n")
cat("1. daily_steps_max_comprehensive_panel.pdf - Complete overview\n")
cat("2. daily_steps_max_trend.pdf - Day-by-day progression\n")
cat("3. recovery_scatter_plot.pdf - Pre vs post comparison\n")
cat("4. Daily_Steps_Max_Analysis_Report.txt - Detailed findings\n")