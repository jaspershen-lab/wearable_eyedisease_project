library(tidyverse)
library(lme4)      # For mixed effects models
library(nlme)      # For extended mixed models capabilities
library(emmeans)   # For marginal means analysis
library(car)       # For Anova function
library(ggstatsplot)
library(ggpubr)
rm(list = ls())
setwd(get_project_wd())
source('1_code/100_tools.R')

#########################
# Step 1: Load data
#########################

# Load steps data
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result_assigned.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

#########################
# Step 2: Calculate device coverage
#########################

# Process baseline date information
baseline_info <- baseline_info %>% mutate(surgery_time_1 = as.Date(surgery_time_1))
hr_df <- heart_rate_data@sample_info %>% as.data.frame()

# Calculate daily hourly coverage
coverage <- hr_df %>%
  left_join(baseline_info %>% dplyr::select(ID, surgery_time_1), by = c("subject_id" = "ID")) %>%
  mutate(
    days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
    day_point = floor(days_to_surgery),
    hour = lubridate::hour(measure_time)
  ) %>%
  filter(day_point >= -7, day_point <= 30) %>%
  group_by(subject_id, day_point) %>%
  summarise(hours_covered = n_distinct(hour), .groups = "drop")

# Create complete grid
all_subjects <- unique(hr_df$subject_id)
all_days <- seq(-7, 30)
complete_grid <- expand.grid(subject_id = all_subjects, day_point = all_days, stringsAsFactors = FALSE)

# Merge and determine if coverage meets threshold
daily_coverage <- complete_grid %>%
  left_join(coverage, by = c("subject_id", "day_point")) %>%
  mutate(hours_covered = coalesce(hours_covered, 0),
         meets_threshold = hours_covered >= 8)

#########################
# Step 3: Filter data based on coverage
#########################

# Get subjects and days that meet coverage requirements
filtered_days <- daily_coverage %>% 
  filter(meets_threshold) %>%
  dplyr::select(subject_id, day_point)

# Create wide format filter
wide_filter <- filtered_days %>%
  mutate(include = TRUE) %>%
  pivot_wider(
    id_cols = subject_id,
    names_from = day_point,
    values_from = include,
    names_prefix = "day_",
    values_fill = FALSE
  )

# Filter steps data
demo_cols <- c("subject_id", "surgery_date")
day_cols <- setdiff(names(daily_steps_result), demo_cols)

filtered_steps_data <- daily_steps_result %>%
  inner_join(wide_filter %>% dplyr::select(subject_id), by = "subject_id")

# Apply coverage filtering to each day column
for(col in day_cols) {
  day <- as.numeric(str_extract(col, "-?\\d+"))
  day_col_name <- paste0("day_", day)
  
  if(day_col_name %in% names(wide_filter)) {
    filtered_steps_data <- filtered_steps_data %>%
      left_join(wide_filter %>% dplyr::select(subject_id, !!sym(day_col_name)), by = "subject_id") %>%
      mutate(!!col := ifelse(!!sym(day_col_name), !!sym(col), NA_real_)) %>%
      dplyr::select(-!!sym(day_col_name))
  }
}

#########################
# Step 4: Create diabetes group information
#########################

# Create diabetes status groups
dm_groups <- baseline_info %>%
  mutate(
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dm_status)) %>%
  mutate(dm_group = dm_status) %>%
  dplyr::select(ID, dm_group, age, gender, bmi)

#########################
# Step 5: Define helper functions for visualization
#########################

# Function to extract contrasts with confidence intervals
extract_contrasts <- function(model, group_var) {
  # Get marginal means for each group and period
  emm <- emmeans(model, specs = as.formula(paste0("~ period | ", group_var)))
  
  # Create contrasts (post - pre)
  contrasts <- contrast(emm, method = "revpairwise") 
  
  # Convert to dataframe
  contrast_df <- as.data.frame(contrasts)
  
  # Add group info
  contrast_df$group <- contrast_df[[group_var]]
  
  # Calculate 95% confidence intervals
  contrast_df$t_crit <- qt(0.975, contrast_df$df)
  contrast_df$lower.CL <- contrast_df$estimate - contrast_df$t_crit * contrast_df$SE
  contrast_df$upper.CL <- contrast_df$estimate + contrast_df$t_crit * contrast_df$SE
  
  return(contrast_df)
}

# Function to add significance markers
add_significance <- function(contrast_df) {
  contrast_df$significance <- ""
  contrast_df$significance[contrast_df$p.value < 0.05] <- "*"
  contrast_df$significance[contrast_df$p.value < 0.01] <- "**"
  contrast_df$significance[contrast_df$p.value < 0.001] <- "***"
  return(contrast_df)
}

# Function to create forest plot with significance markers
create_forest_plot_with_sig <- function(contrast_df, group_col, colors_map, 
                                        title, metric_name = "Steps") {
  # Add space for significance markers
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  p <- ggplot(contrast_df, aes(x = estimate, y = group, color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
    # Add significance markers
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    scale_color_manual(values = colors_map) +
    labs(
      title = title,
      subtitle = "Positive value indicates an increase after surgery",
      x = paste0(metric_name, " post - pre (mean & CI 95%)"),
      y = "",
      caption = "Adjusted for age, gender, and BMI\n* p<0.05, ** p<0.01, *** p<0.001"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  return(p)
}

# Function to extract group comparisons with confidence intervals
extract_group_comparison <- function(model, group_var, reference_group) {
  # Get marginal means for each group and period
  emm <- emmeans(model, specs = c(group_var, "period"))
  
  # Create contrasts between groups within each period
  contrasts <- contrast(emm, method = "revpairwise", by = "period") 
  
  # Convert to dataframe
  contrast_df <- as.data.frame(contrasts)
  
  # Filter to only include comparisons against the reference group
  contrast_df <- contrast_df %>%
    filter(grepl(reference_group, contrast))
  
  # Add more readable comparison names
  contrast_df$comparison <- gsub(paste0(" - ", reference_group), "", contrast_df$contrast)
  
  # Calculate 95% confidence intervals
  contrast_df$t_crit <- qt(0.975, contrast_df$df)
  contrast_df$lower.CL <- contrast_df$estimate - contrast_df$t_crit * contrast_df$SE
  contrast_df$upper.CL <- contrast_df$estimate + contrast_df$t_crit * contrast_df$SE
  
  # Add significance markers
  contrast_df$significance <- ""
  contrast_df$significance[contrast_df$p.value < 0.05] <- "*"
  contrast_df$significance[contrast_df$p.value < 0.01] <- "**"
  contrast_df$significance[contrast_df$p.value < 0.001] <- "***"
  
  return(contrast_df)
}

# Function to create comparison forest plot
create_comparison_forest_plot <- function(contrast_df, colors_map, reference_group,
                                          title, metric_name = "Steps") {
  # Ensure period is an ordered factor
  contrast_df$period <- factor(contrast_df$period, levels = c("Pre-surgery", "Post-surgery"))
  
  # Add space for significance markers
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  # Create forest plot
  p <- ggplot(contrast_df, aes(x = estimate, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3, color = colors_map["Diabetes"]) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, 
                   color = colors_map["Diabetes"]) +
    # Add significance markers
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    facet_wrap(~ period, ncol = 1) +
    labs(
      title = title,
      subtitle = paste0("Positive value indicates higher ", metric_name, " compared to ", reference_group),
      x = paste0("Groups - ", reference_group, " ", metric_name, " (mean & CI 95%)"),
      y = "",
      caption = "Adjusted for age, gender, and BMI\n* p<0.05, ** p<0.01, *** p<0.001"
    ) +
    theme_bw() +
    theme(
      legend.position = "none", 
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  return(p)
}

#########################
# Step 6: Create output directory
#########################
dir.create("3_data_analysis/4_correlation_analysis/steps/steps_mixed_effects_models/diabetes_group", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/4_correlation_analysis/steps/steps_mixed_effects_models/diabetes_group")

#########################
# Step 7: Analyze each steps metric by diabetes group
#########################

# Function to analyze and visualize a specific steps metric
analyze_and_visualize_steps <- function(metric_type) {
  # Merge steps data with diabetes group information
  combined_data <- filtered_steps_data %>%
    left_join(dm_groups, by = c("subject_id" = "ID")) %>%
    filter(!is.na(dm_group))
  
  # Convert to long format
  demo_cols <- c("subject_id", "surgery_date", "dm_group", "age", "gender", "bmi")
  step_cols <- setdiff(names(combined_data), demo_cols)
  
  long_data <- combined_data %>%
    pivot_longer(
      cols = all_of(step_cols),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      day = as.numeric(str_extract(variable, "-?\\d+")),
      metric_type_extracted = str_extract(variable, "(total|mean|max|median)"),
      period = ifelse(day < 0, "Pre-surgery", "Post-surgery"),
      gender = as.factor(gender),
      bmi = as.numeric(bmi),
      age = as.numeric(age)
    ) %>%
    filter(!is.na(value), day >= -7, day <= 30, 
           metric_type_extracted == metric_type)
  
  # Fit mixed effects model
  lme_model <- tryCatch({
    lme(value ~ dm_group * period + age + gender + bmi, 
        data = long_data,
        random = ~1|subject_id,
        weights = varIdent(form = ~1|period),
        method = "REML",
        na.action = na.omit)
  }, error = function(e) {
    # Fallback to a simpler model if the original fails
    message("Error in mixed model, trying simpler version: ", e$message)
    lme(value ~ dm_group * period + age + gender + bmi, 
        data = long_data,
        random = ~1|subject_id,
        method = "REML",
        na.action = na.omit)
  })
  
  # Calculate marginal means
  emm_interaction <- emmeans(lme_model, specs = ~ dm_group | period)
  
  # 1. Time trend plot
  avg_data <- long_data %>%
    group_by(day, dm_group) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      lower_ci = mean_value - 1.96 * se_value,
      upper_ci = mean_value + 1.96 * se_value
    )
  
  # Draw trend plot
  p_trend <- ggplot(avg_data, aes(x = day, y = mean_value, color = dm_group, group = dm_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = dm_group), alpha = 0.2, color = NA) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    scale_color_manual(values = diabetes_colors) + 
    scale_fill_manual(values = diabetes_colors) +
    labs(
      title = paste0("Steps (", metric_type, ") by diabetes group"),
      # subtitle = "Mixed-effects model adjusted for age, gender and BMI",
      x = "Days relative to surgery",
      y = paste0("Steps (", metric_type, ")"),
      color = "Group",
      fill = "Group"
    ) + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  # Save the time trend plot
  ggsave(paste0("dm_group_", metric_type, "_time_trend.pdf"), p_trend, width = 10, height = 6, dpi = 300)
  
  # 2. Create box plot for group comparisons
  dm_plot <- ggplot(long_data, aes(x = dm_group, y = value, fill = dm_group)) +
    # Add boxplot
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    # Add jittered points
    geom_jitter(aes(color = dm_group), width = 0.2, alpha = 0.5, size = 2) +
    # Add mean point in red
    stat_summary(fun = mean, geom = "point", shape = 16, size = 5, color = "darkred") +
    # Add horizontal line at mean with label
    stat_summary(fun = mean, geom = "text", 
                 aes(label = sprintf("μ = %.0f", ..y..)),
                 hjust = -0.3, vjust = -0.5) +
    # Use custom colors
    scale_fill_manual(values = diabetes_colors) +
    scale_color_manual(values = diabetes_colors) +
    # Apply bw theme with customizations
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    # Add labels
    labs(
      title = paste0("Steps (", metric_type, ") by Diabetes Group"),
      x = "Diabetes Status",
      y = paste0("Steps (", metric_type, ")"),
      fill = "Group",
      color = "Group"
    ) +
    # Add statistical comparison
    stat_compare_means(method = "wilcox.test",
                       label.y = max(long_data$value, na.rm = TRUE) * 1.1)
  
  # Save the box plot
  ggsave(paste0("dm_group_", metric_type, "_stats_comparison.pdf"), dm_plot, width = 10, height = 6, dpi = 300)
  
  # 3. Create forest plot for pre-post differences
  dm_contrasts <- extract_contrasts(lme_model, "dm_group")
  dm_contrasts <- add_significance(dm_contrasts)
  
  dm_forest_sig <- create_forest_plot_with_sig(
    dm_contrasts, 
    "dm_group", 
    diabetes_colors,
    paste0("Steps (", metric_type, "): Preoperative and postoperative differences by diabetes status"),
    paste0("Steps (", metric_type, ")")
  )
  
  # Save the forest plot
  ggsave(paste0("dm_group_", metric_type, "_forest_plot_period.pdf"), dm_forest_sig, width = 10, height = 6, dpi = 300)
  
  # 4. Create between-group comparison forest plot
  dm_comparison <- extract_group_comparison(lme_model, "dm_group", "No Diabetes")
  
  dm_forest <- create_comparison_forest_plot(
    dm_comparison,
    diabetes_colors,
    "No Diabetes",
    paste0("Diabetes vs Non-Diabetes Steps (", metric_type, ") Differences"),
    paste0("Steps (", metric_type, ")")
  )
  
  # Save the comparison forest plot
  ggsave(paste0("dm_group_", metric_type, "_comparison_forest.pdf"), dm_forest, width = 10, height = 8, dpi = 300)
  
  # Return model and data
  return(list(model = lme_model, data = long_data))
}

#########################
# Step 8: Apply analysis to total and max steps metrics
#########################

# Analyze total steps
total_steps_results <- analyze_and_visualize_steps("total")
print("Completed total steps analysis and visualization")

# Analyze max steps
max_steps_results <- analyze_and_visualize_steps("max")
print("Completed max steps analysis and visualization")

# Print final message
print("Analysis complete. All visualizations saved in the output directory.")

#########################
# Step 9: Create summary table
#########################

# Function to create a concise summary of findings
create_summary_table <- function() {
  models <- list(
    total = total_steps_results$model,
    max = max_steps_results$model
  )
  
  # Function to extract key statistics for a single model
  extract_model_stats <- function(model_name, model) {
    # Extract interaction term p-value
    anova_result <- anova(model)
    interaction_row <- which(rownames(anova_result) == "dm_group:period")
    interaction_p <- anova_result[interaction_row, "p-value"]
    
    # Extract pre-post differences for each group
    dm_contrasts <- extract_contrasts(model, "dm_group")
    
    # Format p-values with significance markers
    format_p <- function(p) {
      sig <- ""
      if (p < 0.05) sig <- "*"
      if (p < 0.01) sig <- "**" 
      if (p < 0.001) sig <- "***"
      sprintf("%.3f%s", p, sig)
    }
    
    # Get contrasts for each group
    diabetes_row <- which(dm_contrasts$dm_group == "Diabetes")
    no_diabetes_row <- which(dm_contrasts$dm_group == "No Diabetes")
    
    # Create a row of summary data
    data.frame(
      Metric = model_name,
      Interaction_P = format_p(interaction_p),
      Diabetes_Pre_Post_Diff = sprintf("%.1f", dm_contrasts$estimate[diabetes_row]),
      Diabetes_P = format_p(dm_contrasts$p.value[diabetes_row]),
      No_Diabetes_Pre_Post_Diff = sprintf("%.1f", dm_contrasts$estimate[no_diabetes_row]),
      No_Diabetes_P = format_p(dm_contrasts$p.value[no_diabetes_row])
    )
  }
  
  # Apply to each model and combine results
  result <- do.call(rbind, mapply(extract_model_stats, 
                                  names(models), 
                                  models, 
                                  SIMPLIFY = FALSE))
  
  # Save as CSV
  write.csv(result, "steps_diabetes_summary.csv", row.names = FALSE)
  
  return(result)
}

# Create and print summary table
summary_table <- create_summary_table()
print(summary_table)