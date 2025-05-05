library(tidyverse)
library(lme4)      # For mixed effects models
library(nlme)      # For extended mixed models capabilities
library(emmeans)   # For marginal means analysis
library(car)       # For Anova function
library(ggstatsplot)
rm(list = ls())
setwd(get_project_wd())
source('1_code/100_tools.R')


#########################
# Step 1: Load data
#########################

load("3_data_analysis/2_data_analysis/bo/bo_time_periods/daily_bo_result.rda")
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

# Filter blood oxygen data
demo_cols <- c("subject_id", "surgery_date")
day_cols <- setdiff(names(daily_bo_result), demo_cols)

filtered_bo_data <- daily_bo_result %>%
  inner_join(wide_filter %>% dplyr::select(subject_id), by = "subject_id")

# Apply coverage filtering to each day column
for(col in day_cols) {
  day <- as.numeric(str_extract(col, "-?\\d+"))
  day_col_name <- paste0("day_", day)
  
  if(day_col_name %in% names(wide_filter)) {
    filtered_bo_data <- filtered_bo_data %>%
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
# Step 5: Define helper functions
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
                                        title, metric_name = "BO（mean）") {
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

#########################
# Step 6: Create output directory
#########################
dir.create("3_data_analysis/4_correlation_analysis/BO/BO_mixed_effects_models/diabetes_group", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/4_correlation_analysis/BO/BO_mixed_effects_models/diabetes_group")

#########################
# Step 7: Analyze diabetes groups - mean metrics
#########################
# Merge data
combined_data <- filtered_bo_data %>%
  left_join(dm_groups, by = c("subject_id" = "ID")) %>%
  filter(!is.na(dm_group))

# Convert to long format
demo_cols <- c("subject_id", "surgery_date", "dm_group", "age", "gender", "bmi")
bo_cols <- setdiff(names(combined_data), demo_cols)

long_data <- combined_data %>%
  pivot_longer(
    cols = all_of(bo_cols),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    day = as.numeric(str_extract(variable, "-?\\d+")),
    stat_type = str_extract(variable, "(mean|min|max|median|sd|cv|iqr|skew|kurt)"),
    period = ifelse(day < 0, "Pre-surgery", "Post-surgery"),
    gender = as.factor(gender),
    bmi = as.numeric(bmi),
    age = as.numeric(age)
  ) %>%
  filter(!is.na(value), day >= -7, day <= 30)

# Analyze mean statistics
mean_data <- long_data %>% filter(stat_type == "mean")

# Fit mixed effects model
lme_model <- lme(value ~ dm_group * period + age + gender + bmi, 
                 data = mean_data,
                 random = ~1|subject_id,
                 weights = varIdent(form = ~1|period),
                 method = "REML",
                 na.action = na.omit)

# Calculate marginal means
emm_interaction <- emmeans(lme_model, specs = ~ dm_group | period)

#########################
# Step 8: Create visualizations
#########################

# 1. Time trend plot
avg_data <- mean_data %>%
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
  geom_line(size = 1) + geom_point(size = 2) +
  scale_color_manual(values = diabetes_colors) + 
  scale_fill_manual(values = diabetes_colors) +
  labs(
    title = "Mean blood oxygen levels by diabetes group",
    x = "Days relative to surgery",
    y = "Blood oxygen (mean)"
  ) + theme_bw() + theme(legend.position = "bottom")

p_trend
ggsave("dm_group_mean_time_trend.pdf", p_trend, width = 10, height = 6, dpi = 300)

# 2. Create box plot for group comparisons
# Creating a custom boxplot with ggplot2 instead of ggbetweenstats

# First, make sure your data is properly structured
# If you're getting viewport errors, try basic ggplot2 functions

# Here's a solution for your data:
dm_plot <- ggplot(mean_data, aes(x = dm_group, y = value, fill = dm_group)) +
  # Add boxplot
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  # Add jittered points
  geom_jitter(aes(color = dm_group), width = 0.2, alpha = 0.5, size = 2) +
  # Add mean point in red
  stat_summary(fun = mean, geom = "point", shape = 16, size = 5, color = "darkred") +
  # Add horizontal line at mean with label
  stat_summary(fun = mean, geom = "text", 
               aes(label = sprintf("μmean = %.2f", ..y..)),
               hjust = -0.3, vjust = -0.5) +
  # Use your custom colors
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
    title = "Blood Oxygen Levels by Diabetes Group",
    x = "Diabetes Status",
    y = "Blood Oxygen (mean)"
  )

library(ggpubr)

dm_plot<-dm_plot +
  stat_compare_means(method = "wilcox.test",
                    label.y = max(mean_data$value, na.rm = TRUE) + 0.5)


ggsave("dm_group_mean_stats_comparison.pdf", dm_plot, width = 10, height = 6, dpi = 300)

# 3. Create forest plot for pre-post differences
dm_contrasts <- extract_contrasts(lme_model, "dm_group")
dm_contrasts <- add_significance(dm_contrasts)

dm_forest_sig <- create_forest_plot_with_sig(
  dm_contrasts, 
  "dm_group", 
  diabetes_colors,
  "BO level: Preoperative and postoperative differences by diabetes status"
)

dm_forest_sig
ggsave("dm_group_forest_plot_period.pdf", dm_forest_sig, width = 10, height = 6, dpi = 300)

# 4. Create between-group comparison forest plot
# Modify your create_comparison_forest_plot function to handle the color mapping correctly:
create_comparison_forest_plot <- function(contrast_df, colors_map, reference_group,
                                          title, metric_name = "Blood Oxygen (mean)") {
  # Ensure period is an ordered factor
  contrast_df$period <- factor(contrast_df$period, levels = c("Pre-surgery", "Post-surgery"))
  
  # Add space for significance markers
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  # Create forest plot
  p <- ggplot(contrast_df, aes(x = estimate, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    # Use direct color instead of mapping
    geom_point(size = 3, color = colors_map["Diabetes"]) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, 
                   color = colors_map["Diabetes"]) +
    # Add significance markers
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    facet_wrap(~ period, ncol = 1) +
    labs(
      title = title,
      subtitle = paste0("Positive value indicates higher BO level compared to ", reference_group),
      x = paste0("Groups - ", reference_group, " ", metric_name, " (mean & CI 95%)"),
      y = "",
      caption = "Mixed-effects model adjusted for age, gender, and BMI\n* p<0.05, ** p<0.01, *** p<0.001"
    ) +
    theme_bw() +
    theme(
      legend.position = "none", # Remove legend since we're using fixed colors
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
  # Using reference_group as the reference level
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

# Extract comparison
dm_comparison <- extract_group_comparison(lme_model, "dm_group", "No Diabetes")

# Create forest plot
dm_forest <- create_comparison_forest_plot(
  dm_comparison,
  diabetes_colors,
  "No Diabetes",
  "Diabetes vs Non-Diabetes Blood Oxygen Level Differences"
)
dm_forest
ggsave("dm_group_comparison_forest.pdf", dm_forest, width = 10, height = 8, dpi = 300)

#########################
# Step 9: Save results
#########################

# Save model and data for potential further analysis
save(lme_model, mean_data, file = "dm_group_analysis_results.RData")
