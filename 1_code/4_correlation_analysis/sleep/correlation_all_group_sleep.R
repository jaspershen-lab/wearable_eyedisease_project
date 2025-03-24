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

# Load sleep data
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/daily_sleep_result.rda")
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

# Filter sleep data
demo_cols <- c("subject_id", "surgery_date")
day_cols <- setdiff(names(daily_sleep_result), demo_cols)

filtered_sleep_data <- daily_sleep_result %>%
  inner_join(wide_filter %>% dplyr::select(subject_id), by = "subject_id")

# Apply coverage filtering to each day column
for(col in day_cols) {
  day <- as.numeric(str_extract(col, "-?\\d+"))
  day_col_name <- paste0("day_", day)
  
  if(day_col_name %in% names(wide_filter)) {
    filtered_sleep_data <- filtered_sleep_data %>%
      left_join(wide_filter %>% dplyr::select(subject_id, !!sym(day_col_name)), by = "subject_id") %>%
      mutate(!!col := ifelse(!!sym(day_col_name), !!sym(col), NA_real_)) %>%
      dplyr::select(-!!sym(day_col_name))
  }
}

#########################
# Step 4: Create surgery group information (3 groups)
#########################

# Create three original groups
original_groups <- baseline_info %>%
  mutate(
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  mutate(
    orig_group = case_when(
      dm_status == "No Diabetes" & surgery_type == 0 ~ "No_Diabetes_Surgery0",
      dm_status == "Diabetes" & surgery_type == 0 ~ "Diabetes_Surgery0", 
      dm_status == "Diabetes" & surgery_type == 1 ~ "Diabetes_Surgery1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(orig_group)) %>%
  dplyr::select(ID, orig_group, age, gender, bmi)


# Display labels for better readability
group_labels <- c(
  "No_Diabetes_Surgery0" = "Non-DM Cataract",
  "Diabetes_Surgery0" = "DM Cataract",
  "Diabetes_Surgery1" = "DM PPV"
)

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
create_forest_plot_with_sig <- function(contrast_df, group_col, colors_map, labels_map,
                                        title, metric_name = "Sleep") {
  # Add space for significance markers
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  # Add display names
  contrast_df$display_name <- labels_map[contrast_df$group]
  
  p <- ggplot(contrast_df, aes(x = estimate, y = display_name, color = group)) +
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
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  return(p)
}

# Function to extract paired group comparisons
extract_paired_comparisons <- function(model, group_var) {
  # Get marginal means
  emm <- emmeans(model, specs = c(group_var, "period"))
  
  # Create pairwise contrasts between all groups
  contrasts <- contrast(emm, method = "pairwise", by = "period") 
  
  # Convert to dataframe
  contrast_df <- as.data.frame(contrasts)
  
  # Extract comparison pairs
  contrast_df$group1 <- gsub(" -.*", "", contrast_df$contrast)
  contrast_df$group2 <- gsub(".*- ", "", contrast_df$contrast)
  
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

# Function to create paired comparison forest plot
create_paired_comparison_forest_plot <- function(contrast_df, colors_map, labels_map,
                                                 title, metric_name = "Sleep") {
  # Ensure period is an ordered factor
  contrast_df$period <- factor(contrast_df$period, levels = c("Pre-surgery", "Post-surgery"))
  
  # Add space for significance markers
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  # Create readable contrast labels with display names
  contrast_df$group1_display <- labels_map[contrast_df$group1]
  contrast_df$group2_display <- labels_map[contrast_df$group2]
  contrast_df$contrast_label <- paste(contrast_df$group1_display, "vs", contrast_df$group2_display)
  
  # Create forest plot
  p <- ggplot(contrast_df, aes(x = estimate, y = contrast_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
    # Add significance markers
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    facet_wrap(~ period, ncol = 1) +
    labs(
      title = title,
      subtitle = paste0("Group differences in ", metric_name),
      x = paste0(metric_name, " difference (mean & CI 95%)"),
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
dir.create("3_data_analysis/4_correlation_analysis/sleep/sleep_mixed_effects_models/surgery_group", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/4_correlation_analysis/sleep/sleep_mixed_effects_models/surgery_group")

#########################
# Step 7: Analyze each sleep metric by surgery group
#########################

# Function to analyze and visualize a specific sleep metric by surgery group
analyze_and_visualize_sleep_by_surgery <- function(metric_type) {
  # Merge sleep data with surgery group information
  combined_data <- filtered_sleep_data %>%
    left_join(original_groups, by = c("subject_id" = "ID")) %>%
    filter(!is.na(orig_group))
  
  # Convert to long format
  demo_cols <- c("subject_id", "surgery_date", "orig_group", "age", "gender", "bmi")
  sleep_cols <- setdiff(names(combined_data), demo_cols)
  
  long_data <- combined_data %>%
    pivot_longer(
      cols = all_of(sleep_cols),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      day = as.numeric(str_extract(variable, "-?\\d+")),
      metric_type_extracted = str_extract(variable, "(light_sleep|deep_sleep|dream_sleep|awake|total_sleep|daytime_sleep|sleep_score)"),
      period = ifelse(day < 0, "Pre-surgery", "Post-surgery"),
      gender = as.factor(gender),
      bmi = as.numeric(bmi),
      age = as.numeric(age)
    ) %>%
    filter(!is.na(value), day >= -7, day <= 30, 
           metric_type_extracted == metric_type)
  
  # Create display labels for groups in data
  long_data <- long_data %>%
    mutate(group_label = group_labels[orig_group])
  
  # Set order for group display
  long_data$group_label <- factor(long_data$group_label, 
                                  levels = c("Non-DM Cataract", "DM Cataract", "DM PPV"))
  
  # Fit mixed effects model
  lme_model <- tryCatch({
    lme(value ~ orig_group * period + age + gender + bmi, 
        data = long_data,
        random = ~1|subject_id,
        weights = varIdent(form = ~1|period),
        method = "REML",
        na.action = na.omit)
  }, error = function(e) {
    # Fallback to a simpler model if the original fails
    message("Error in mixed model, trying simpler version: ", e$message)
    lme(value ~ orig_group * period + age + gender + bmi, 
        data = long_data,
        random = ~1|subject_id,
        method = "REML",
        na.action = na.omit)
  })
  
  # Calculate marginal means
  emm_interaction <- emmeans(lme_model, specs = ~ orig_group | period)
  
  # Make friendly metric name for display
  display_name <- case_when(
    metric_type == "total_sleep" ~ "Total Sleep",
    metric_type == "deep_sleep" ~ "Deep Sleep",
    TRUE ~ metric_type
  )
  
  # 1. Time trend plot
  avg_data <- long_data %>%
    group_by(day, group_label) %>%
    summarize(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      lower_ci = mean_value - 1.96 * se_value,
      upper_ci = mean_value + 1.96 * se_value
    )
  
  # Map display labels to colors
  display_colors <- c(
    "Non-DM Cataract" = original_colors["No_Diabetes_Surgery0"],
    "DM Cataract" = original_colors["Diabetes_Surgery0"],
    "DM PPV" = original_colors["Diabetes_Surgery1"]
  )
  
  # Draw trend plot
  p_trend <- ggplot(avg_data, aes(x = day, y = mean_value, color = group_label, group = group_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = group_label), alpha = 0.2, color = NA) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    scale_color_manual(values = display_colors) + 
    scale_fill_manual(values = display_colors) +
    labs(
      title = paste0(display_name, " by surgery group"),
      subtitle = "Mixed-effects model adjusted for age, gender and BMI",
      x = "Days relative to surgery",
      y = paste0(display_name, " (minutes)"),
      color = "Group",
      fill = "Group"
    ) + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  # Save the time trend plot
  ggsave(paste0("surgery_group_", metric_type, "_time_trend.pdf"), p_trend, width = 10, height = 6, dpi = 300)
  
  # 2. Create box plot for group comparisons
  surgery_plot <- ggplot(long_data, aes(x = group_label, y = value, fill = group_label)) +
    # Add boxplot
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    # Add jittered points
    geom_jitter(aes(color = group_label), width = 0.2, alpha = 0.5, size = 2) +
    # Add mean point in red
    stat_summary(fun = mean, geom = "point", shape = 16, size = 5, color = "darkred") +
    # Add horizontal line at mean with label
    stat_summary(fun = mean, geom = "text", 
                 aes(label = sprintf("μ = %.0f", ..y..)),
                 hjust = -0.3, vjust = -0.5) +
    # Use custom colors
    scale_fill_manual(values = display_colors) +
    scale_color_manual(values = display_colors) +
    # Apply bw theme with customizations
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    # Add labels
    labs(
      title = paste0(display_name, " by Surgery Group"),
      x = "Surgery Group",
      y = paste0(display_name, " (minutes)"),
      fill = "Group",
      color = "Group"
    ) +
    # Add ANOVA comparison
    stat_compare_means(method = "anova",
                       label.y = max(long_data$value, na.rm = TRUE) * 1.1)
  
  # Save the box plot
  ggsave(paste0("surgery_group_", metric_type, "_stats_comparison.pdf"), surgery_plot, 
         width = 10, height = 6, dpi = 300)
  
  # 3. Create forest plot for pre-post differences
  group_contrasts <- extract_contrasts(lme_model, "orig_group")
  group_contrasts <- add_significance(group_contrasts)
  
  group_forest_sig <- create_forest_plot_with_sig(
    group_contrasts, 
    "group", 
    original_colors,
    group_labels,
    paste0(display_name, ": Preoperative and postoperative differences by surgery group"),
    paste0(display_name, " (minutes)")
  )
  
  # Save the forest plot
  ggsave(paste0("surgery_group_", metric_type, "_forest_plot_period.pdf"), group_forest_sig, 
         width = 10, height = 6, dpi = 300)
  
  # 4. Create pairwise comparison forest plot
  paired_comparisons <- extract_paired_comparisons(lme_model, "orig_group")
  
  paired_forest <- create_paired_comparison_forest_plot(
    paired_comparisons,
    original_colors,
    group_labels,
    paste0(display_name, ": Pairwise Group Differences"),
    paste0(display_name, " (minutes)")
  )
  
  # Save the paired comparison forest plot
  ggsave(paste0("surgery_group_", metric_type, "_paired_comparisons.pdf"), paired_forest, 
         width = 12, height = 8, dpi = 300)
  
  # Return model and data
  return(list(model = lme_model, data = long_data))
}

#########################
# Step 8: Apply analysis to total sleep and deep sleep metrics
#########################

# Analyze total sleep
total_sleep_results <- analyze_and_visualize_sleep_by_surgery("total_sleep")
print("Completed total sleep analysis by surgery group")

# Analyze deep sleep
deep_sleep_results <- analyze_and_visualize_sleep_by_surgery("deep_sleep")
print("Completed deep sleep analysis by surgery group")

# Print final message
print("Analysis complete. All visualizations saved in the output directory.")

#########################
# Step 9: Create summary table
#########################

# Function to create a concise summary of findings
create_summary_table <- function() {
  metrics <- c("total_sleep", "deep_sleep")
  models <- list(
    total_sleep = total_sleep_results$model,
    deep_sleep = deep_sleep_results$model
  )
  
  # Function to extract key statistics for interaction
  extract_interaction <- function(model_name, model) {
    # Extract interaction term p-value
    anova_result <- anova(model)
    interaction_row <- which(rownames(anova_result) == "orig_group:period")
    interaction_p <- anova_result[interaction_row, "p-value"]
    
    # Format p-value with significance marker
    format_p <- function(p) {
      sig <- ""
      if (p < 0.05) sig <- "*"
      if (p < 0.01) sig <- "**" 
      if (p < 0.001) sig <- "***"
      sprintf("%.3f%s", p, sig)
    }
    
    # Create a row of summary data
    data.frame(
      Metric = case_when(
        model_name == "total_sleep" ~ "Total Sleep",
        model_name == "deep_sleep" ~ "Deep Sleep",
        TRUE ~ model_name
      ),
      Group_Period_Interaction_P = format_p(interaction_p)
    )
  }
  
  # Function to extract pre-post contrasts for each group
  extract_group_contrasts <- function(model_name, model) {
    # Extract pre-post differences for each group
    group_contrasts <- extract_contrasts(model, "orig_group")
    
    # Format p-values with significance markers
    format_p <- function(p) {
      sig <- ""
      if (p < 0.05) sig <- "*"
      if (p < 0.01) sig <- "**" 
      if (p < 0.001) sig <- "***"
      sprintf("%.3f%s", p, sig)
    }
    
    # Get rows for each group
    nd_s0_row <- which(group_contrasts$orig_group == "No_Diabetes_Surgery0")
    dm_s0_row <- which(group_contrasts$orig_group == "Diabetes_Surgery0")
    dm_s1_row <- which(group_contrasts$orig_group == "Diabetes_Surgery1")
    
    # Create data frame with results
    data.frame(
      Metric = case_when(
        model_name == "total_sleep" ~ "Total Sleep",
        model_name == "deep_sleep" ~ "Deep Sleep",
        TRUE ~ model_name
      ),
      NonDM_Cat_Diff = sprintf("%.1f", group_contrasts$estimate[nd_s0_row]),
      NonDM_Cat_P = format_p(group_contrasts$p.value[nd_s0_row]),
      DM_Cat_Diff = sprintf("%.1f", group_contrasts$estimate[dm_s0_row]),
      DM_Cat_P = format_p(group_contrasts$p.value[dm_s0_row]),
      DM_PPV_Diff = sprintf("%.1f", group_contrasts$estimate[dm_s1_row]),
      DM_PPV_P = format_p(group_contrasts$p.value[dm_s1_row])
    )
  }
  
  # Apply to each model and combine results
  interaction_results <- do.call(rbind, mapply(extract_interaction, 
                                               names(models), 
                                               models, 
                                               SIMPLIFY = FALSE))
  
  contrast_results <- do.call(rbind, mapply(extract_group_contrasts,
                                            names(models),
                                            models,
                                            SIMPLIFY = FALSE))
  
  # Combine the results
  result <- interaction_results %>%
    left_join(contrast_results, by = "Metric")
  
  # Save as CSV
  write.csv(result, "sleep_surgery_group_summary.csv", row.names = FALSE)
  
  return(result)
}

# Create and print summary table
summary_table <- create_summary_table()
print(summary_table)