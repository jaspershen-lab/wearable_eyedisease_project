# Build model formula
formula_str <- "value ~ group * day + age + gender+ bmi + (1|subject_id)"
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
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result_assigned.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")


##########################
# 2. Data Preparation
##########################
# Load heart rate data to get heart rate IDs
load("3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda")
# Get unique IDs from heart rate data
heart_rate_ids <- unique(heart_rate_data@sample_info$subject_id)

# Get unique IDs from steps data
steps_ids <- unique(daily_steps_result$subject_id)

# Find participants who have both heart rate and steps data
common_ids <- intersect(heart_rate_ids, steps_ids)

# Filter baseline info to include only participants with both heart rate and steps data
baseline_info_filtered <- baseline_info %>% 
  filter(ID %in% common_ids)
cat("Number of participants in filtered baseline data:", nrow(baseline_info_filtered), "\n")
cat("Number of participants with heart rate data:", length(heart_rate_ids), "\n")
cat("Number of participants with steps data:", length(steps_ids), "\n")
cat("Number of participants with both heart rate and steps data:", length(common_ids), "\n")

# Create group information based on diabetes status and surgery type
group_info <- baseline_info_filtered %>%
  mutate(
    # Create diabetes status based on diabetes_history
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    # Use surgery_1..0.PI.1.other. as surgery type
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  dplyr::select(ID, dm_status, surgery_type)

# Display distribution of groups
cat("Diabetes status and surgery type distribution:\n")
print(table(group_info$dm_status, group_info$surgery_type, dnn = c("Diabetes Status", "Surgery Type")))

# Create group variable
baseline_info_filtered <- baseline_info_filtered %>%
  left_join(group_info, by = "ID") %>%
  mutate(
    # Create a combined group variable
    group = case_when(
      dm_status == "No Diabetes" & surgery_type == 0 ~ "No_Diabetes_Surgery0",
      dm_status == "Diabetes" & surgery_type == 1 ~ "Diabetes_Surgery1",
      TRUE ~ "Other"
    )
  )

# Print group distribution
cat("\nGroup distribution:\n")
print(table(baseline_info_filtered$group))

# Create simplified info with ID and group for merging
group_status <- data.frame(
  subject_id = baseline_info_filtered$ID,
  group = baseline_info_filtered$group,
  dm_status = baseline_info_filtered$dm_status,
  surgery_type = baseline_info_filtered$surgery_type,
  age = baseline_info_filtered$age,
  gender = baseline_info_filtered$gender
)

# Print surgery type distribution
cat("Number of participants by surgery type:\n")
print(table(baseline_info_filtered$surgery_type, useNA = "ifany"))

# Create analysis results directory
dir.create("3_data_analysis/4_correlation_analysis/steps/steps_group_comparison", recursive = TRUE)
setwd("3_data_analysis/4_correlation_analysis/steps/steps_group_comparison")

##########################
# 3. Merge Steps and Surgery Type Data
##########################
# Merge steps data with group information
group_steps_data <- daily_steps_result %>%
  left_join(group_status, by = "subject_id") %>%
  filter(!is.na(group))  # Remove participants with missing group

# Filter to only include the two groups of interest
group_steps_data <- group_steps_data %>%
  filter(group %in% c("No_Diabetes_Surgery0", "Diabetes_Surgery1"))

# Print information about the final dataset
cat("Final number of participants with complete data:", n_distinct(group_steps_data$subject_id), "\n")
cat("Number of participants in No_Diabetes_Surgery0 group:", sum(group_steps_data$group == "No_Diabetes_Surgery0", na.rm = TRUE), "\n")
cat("Number of participants in Diabetes_Surgery1 group:", sum(group_steps_data$group == "Diabetes_Surgery1", na.rm = TRUE), "\n")

##########################
# 4. Convert Wide Data to Long Format
##########################
cat("\nConverting data to long format...\n")

# Identify demographic columns
demo_cols <- c("subject_id", "surgery_date", "group", "dm_status", "surgery_type", "age", "gender")
actual_demo_cols <- names(group_steps_data)[names(group_steps_data) %in% demo_cols]

# Get all steps columns (excluding demographic columns)
steps_cols <- names(group_steps_data)[!names(group_steps_data) %in% actual_demo_cols]

# Create long format data
long_data <- group_steps_data %>%
  pivot_longer(
    cols = all_of(steps_cols),
    names_to = "variable",
    values_to = "value"
  ) %>%
  # Extract day and statistic type
  mutate(
    day = as.numeric(str_extract(variable, "-?\\d+")),
    stat_type = str_extract(variable, "(total|mean|max|median)")
  ) %>%
  # Remove missing values and filter for days -10 to 30
  filter(!is.na(value), day >= -7, day <= 30)

# Print long data summary
cat("Long data summary:\n")
cat("Total observations:", nrow(long_data), "\n")
cat("Unique days:", n_distinct(long_data$day), "\n")
cat("Statistical metrics:", paste(unique(long_data$stat_type), collapse = ", "), "\n")

##########################
# 5. Compare Steps Metrics by Surgery Type
##########################

# Group data by day, stat_type
comparison_results <- long_data %>%
  group_by(day, stat_type) %>%
  # Perform t-test for each group
  summarize(
    n_diabetes_surgery1 = sum(group == "Diabetes_Surgery1", na.rm = TRUE),
    n_no_diabetes_surgery0 = sum(group == "No_Diabetes_Surgery0", na.rm = TRUE),
    mean_diabetes_surgery1 = mean(value[group == "Diabetes_Surgery1"], na.rm = TRUE),
    mean_no_diabetes_surgery0 = mean(value[group == "No_Diabetes_Surgery0"], na.rm = TRUE),
    sd_diabetes_surgery1 = sd(value[group == "Diabetes_Surgery1"], na.rm = TRUE),
    sd_no_diabetes_surgery0 = sd(value[group == "No_Diabetes_Surgery0"], na.rm = TRUE),
    # Perform t-test only if there are at least 2 observations in each group
    p_value = ifelse(
      sum(group == "Diabetes_Surgery1", na.rm = TRUE) >= 2 & 
        sum(group == "No_Diabetes_Surgery0", na.rm = TRUE) >= 2,
      t.test(value ~ group)$p.value,
      NA
    ),
    # Calculate Cohen's d effect size
    cohens_d = ifelse(
      sum(group == "Diabetes_Surgery1", na.rm = TRUE) >= 2 & 
        sum(group == "No_Diabetes_Surgery0", na.rm = TRUE) >= 2,
      (mean(value[group == "Diabetes_Surgery1"], na.rm = TRUE) - 
         mean(value[group == "No_Diabetes_Surgery0"], na.rm = TRUE)) / 
        sqrt((sd(value[group == "Diabetes_Surgery1"], na.rm = TRUE)^2 + 
                sd(value[group == "No_Diabetes_Surgery0"], na.rm = TRUE)^2) / 2),
      NA
    ),
    # Add significance indicator
    significant = p_value < 0.05,
    .groups = "drop"
  )

# Apply FDR correction for multiple testing
comparison_results <- comparison_results %>%
  group_by(stat_type) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "fdr"),
    significant_adjusted = p_adjusted < 0.05
  ) %>%
  ungroup()

# Create summary table of significant findings
significant_findings <- comparison_results %>%
  filter(significant_adjusted == TRUE) %>%
  arrange(stat_type, day)

cat("Found", nrow(significant_findings), "statistically significant differences (FDR-adjusted)\n")

# Save comparison results
write.csv(comparison_results, "group_steps_comparison_results.csv", row.names = FALSE)
write.csv(significant_findings, "group_steps_significant_findings.csv", row.names = FALSE)

##########################
# 6. Visualization: Plot Steps Total by Day
##########################

# Total Steps by day visualization
cat("Creating time series plot for total steps...\n")

# Filter data for total steps
plot_data <- long_data %>%
  filter(
    stat_type == "total"
  )

# Calculate averages by day and group
avg_data <- plot_data %>%
  group_by(day, group) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - 1.96 * se_value,
    upper_ci = mean_value + 1.96 * se_value
  )

# Create plot
p <- ggplot(avg_data, aes(x = day, y = mean_value, color = group, group = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = group), 
              alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("No_Diabetes_Surgery0" = "#3182bd", "Diabetes_Surgery1" = "#e6550d"),
                     labels = c("No Diabetes Surgery 0", "Diabetes Surgery 1")) +
  scale_fill_manual(values = c("No_Diabetes_Surgery0" = "#3182bd", "Diabetes_Surgery1" = "#e6550d"),
                    labels = c("No Diabetes Surgery 0", "Diabetes Surgery 1")) +
  labs(
    title = paste("Mean daily steps by group"),
    subtitle = "Total daily steps",
    x = "Days relative to surgery",
    y = "Steps (total)",
    color = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Display and save plot
print(p)
file_name <- "daily_steps_total_by_group.png"
ggsave(file_name, p, width = 10, height = 6, dpi = 300)
cat("Plot saved as:", file_name, "\n")

# Do the same for mean steps
cat("Creating time series plot for mean steps...\n")

# Filter data for mean steps
plot_data <- long_data %>%
  filter(
    stat_type == "mean"
  )

# Calculate averages by day and group
avg_data <- plot_data %>%
  group_by(day, group) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - 1.96 * se_value,
    upper_ci = mean_value + 1.96 * se_value
  )

# Create plot
p <- ggplot(avg_data, aes(x = day, y = mean_value, color = group, group = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = group), 
              alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("No_Diabetes_Surgery0" = "#3182bd", "Diabetes_Surgery1" = "#e6550d"),
                     labels = c("No Diabetes Surgery 0", "Diabetes Surgery 1")) +
  scale_fill_manual(values = c("No_Diabetes_Surgery0" = "#3182bd", "Diabetes_Surgery1" = "#e6550d"),
                    labels = c("No Diabetes Surgery 0", "Diabetes Surgery 1")) +
  labs(
    title = paste("Mean steps intensity by group"),
    subtitle = "Average steps per measurement",
    x = "Days relative to surgery",
    y = "Steps (mean)",
    color = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Display and save plot
print(p)
file_name <- "daily_steps_mean_by_group.png"
ggsave(file_name, p, width = 10, height = 6, dpi = 300)
cat("Plot saved as:", file_name, "\n")

##########################
# 7. Visualization: Statistical Significance
##########################
# Create p-value plots for different steps metrics
for(stat_type_filter in c("total", "mean", "max", "median")) {
  cat("Creating p-value plot for", stat_type_filter, "...\n")
  
  # Filter data for comparison between Cataract and No Surgery
  plot_data <- comparison_results %>%
    filter(
      stat_type == stat_type_filter,
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
      title = paste("Statistical significance of surgery type effect on steps", stat_type_filter),
      subtitle = "Comparison: Cataract vs No Surgery",
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
  file_name <- paste0("pvalues_", stat_type_filter, "_by_surgery.png")
  ggsave(file_name, p, width = 10, height = 6, dpi = 300)
  cat("Plot saved as:", file_name, "\n")
}

##########################
# 8. Visualization: Period Comparison
##########################

# Create boxplots comparing periods for different steps metrics
for(stat_type_filter in c("total", "mean", "max", "median")) {
  cat("Creating period comparison plot for", stat_type_filter, "...\n")
  
  # Filter data
  plot_data <- long_data %>%
    filter(
      stat_type == stat_type_filter
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
              aes(x = period, y = value, fill = group)) +
    geom_boxplot(outlier.size = 1, width = 0.7, position = position_dodge(0.8)) +
    scale_fill_manual(
      values = c("No_Diabetes_Surgery0" = "#3182bd", "Diabetes_Surgery1" = "#e6550d"),
      labels = c("No Diabetes Surgery 0", "Diabetes Surgery 1")
    ) +
    labs(
      title = paste("Steps", stat_type_filter, "by surgical period and group"),
      x = "Surgical period",
      y = paste("Steps", stat_type_filter),
      fill = "Group"
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
    aes(group = group),
    label = "p.format",
    method = "t.test",
    label.y.npc = 0.95
  )
  
  # Display and save plot
  print(p)
  file_name <- paste0("period_comparison_", stat_type_filter, "_by_group.png")
  ggsave(file_name, p, width = 10, height = 6, dpi = 300)
  cat("Plot saved as:", file_name, "\n")
}

##########################
# 9. Longitudinal Analysis
##########################
# Analyze different metrics
metrics <- list(
  list(stat = "total"),
  list(stat = "mean"),
  list(stat = "max"),
  list(stat = "median")
)

# Run analysis for each metric
longitudinal_results <- list()
for(metric in metrics) {
  cat("Analyzing", metric$stat, "...\n")
  
  # Filter data
  filtered_data <- long_data %>%
    filter(stat_type == metric$stat) %>%
    mutate(
      time_period = case_when(
        day < 0 ~ "Pre-surgery",
        day >= 0 & day <= 30 ~ "Post-surgery",
        TRUE ~ NA_character_
      )
    )
  
  # Build model formula
  formula_str <- "value ~ surgery_type * day + age + gender + (1|subject_id)"
  
  # Try running mixed effects model
  tryCatch({
    model <- lmer(as.formula(formula_str), data = filtered_data)
    
    # Display model results summary
    cat("Fixed effects:\n")
    print(summary(model)$coefficients)
    
    # Get estimated marginal means for surgery type
    emm <- emmeans(model, ~ surgery_type)
    cat("Marginal means:\n")
    print(emm)
    
    contrasts <- contrast(emm, "pairwise")
    cat("Contrasts:\n")
    print(contrasts)
    
    # Save results
    longitudinal_results[[paste(metric$stat)]] <- list(
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
      model <- lmer(value ~ surgery_type + day + age + (1|subject_id), data = filtered_data)
      
      # Display model results summary
      cat("Simplified model fixed effects:\n")
      print(summary(model)$coefficients)
      
      # Get estimated marginal means for surgery type
      emm <- emmeans(model, ~ surgery_type)
      cat("Marginal means:\n")
      print(emm)
      
      contrasts <- contrast(emm, "pairwise")
      cat("Contrasts:\n")
      print(contrasts)
      
      # Save results
      longitudinal_results[[paste(metric$stat)]] <- list(
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
save(group_steps_data, long_data, comparison_results, longitudinal_results, significant_findings,
     file = "group_steps_analysis.rda")

##########################
# 11. Create Summary Report
##########################
cat("\nCreating summary report...\n")

# Analysis parameters summary
cat("Analysis Summary Report\n")
cat("====================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d"), "\n")
cat("Total participants analyzed:", n_distinct(group_steps_data$subject_id), "\n")

# Count participants by group
group_counts <- table(group_steps_data$group)
for(i in 1:length(group_counts)) {
  cat("Participants in", names(group_counts)[i], "group:", group_counts[i], "\n")
}

cat("Days analyzed:", min(long_data$day), "to", max(long_data$day), "(relative to surgery)\n")
cat("Steps metrics analyzed:", paste(unique(long_data$stat_type), collapse = ", "), "\n\n")

# Significant findings summary
cat("Significant Findings Summary\n")
cat("--------------------------\n\n")
cat("Total significant findings (FDR-adjusted):", nrow(significant_findings), "\n\n")

if(nrow(significant_findings) > 0) {
  # Group significant findings by metric
  sig_by_metric <- significant_findings %>%
    group_by(stat_type) %>%
    summarize(
      n_findings = n(),
      mean_effect_size = mean(abs(cohens_d), na.rm = TRUE),
      days = paste(sort(day), collapse = ", "),
      .groups = "drop"
    )
  
  # Print summary by metric
  for(i in 1:nrow(sig_by_metric)) {
    cat(sig_by_metric$stat_type[i], ": ", 
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
    cat(i, ". ", top_findings$stat_type[i], ", Day ", 
        top_findings$day[i], ", Cohen's d = ", round(top_findings$cohens_d[i], 2), 
        ", p = ", format(top_findings$p_adjusted[i], digits = 3), "\n", sep = "")
  }
}

# Save summary report
sink("group_steps_analysis_summary.txt")
cat("Analysis Summary Report\n")
cat("====================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d"), "\n")
cat("Total participants analyzed:", n_distinct(group_steps_data$subject_id), "\n")

# Count participants by group
group_counts <- table(group_steps_data$group)
for(i in 1:length(group_counts)) {
  cat("Participants in", names(group_counts)[i], "group:", group_counts[i], "\n")
}

cat("Days analyzed:", min(long_data$day), "to", max(long_data$day), "(relative to surgery)\n")
cat("Steps metrics analyzed:", paste(unique(long_data$stat_type), collapse = ", "), "\n\n")

cat("Significant Findings Summary\n")
cat("--------------------------\n\n")
cat("Total significant findings (FDR-adjusted):", nrow(significant_findings), "\n\n")

if(nrow(significant_findings) > 0) {
  # Group significant findings by metric
  sig_by_metric <- significant_findings %>%
    group_by(stat_type) %>%
    summarize(
      n_findings = n(),
      mean_effect_size = mean(abs(cohens_d), na.rm = TRUE),
      days = paste(sort(day), collapse = ", "),
      .groups = "drop"
    )
  
  # Print summary by metric
  for(i in 1:nrow(sig_by_metric)) {
    cat(sig_by_metric$stat_type[i], ": ", 
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
    cat(i, ". ", top_findings$stat_type[i], ", Day ", 
        top_findings$day[i], ", Cohen's d = ", round(top_findings$cohens_d[i], 2), 
        ", p = ", format(top_findings$p_adjusted[i], digits = 3), "\n", sep = "")
  }
}
sink()

