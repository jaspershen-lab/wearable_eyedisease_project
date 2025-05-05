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

# Load RHR data
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
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

# Filter RHR data
demo_cols <- c("subject_id", "surgery_date")
day_cols <- setdiff(names(daily_rhr_result), demo_cols)

filtered_rhr_data <- daily_rhr_result %>%
  inner_join(wide_filter %>% dplyr::select(subject_id), by = "subject_id")

# Apply coverage filtering to each day column
for(col in day_cols) {
  day <- as.numeric(str_extract(col, "-?\\d+"))
  day_col_name <- paste0("day_", day)
  
  if(day_col_name %in% names(wide_filter)) {
    filtered_rhr_data <- filtered_rhr_data %>%
      left_join(wide_filter %>% dplyr::select(subject_id, !!sym(day_col_name)), by = "subject_id") %>%
      mutate(!!col := ifelse(!!sym(day_col_name), !!sym(col), NA_real_)) %>%
      dplyr::select(-!!sym(day_col_name))
  }
}

#########################
# Step 4: Create surgery type groups
#########################

# Create non-diabetes cataract vs diabetes PPV groups
nd_cat_vs_dm_ppv_groups <- baseline_info %>%
  mutate(
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  mutate(
    nd_cat_vs_dm_ppv = case_when(
      dm_status == "No Diabetes" & surgery_type == 0 ~ "No_Diabetes_Cataract",
      dm_status == "Diabetes" & surgery_type == 1 ~ "Diabetes_PPV",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(nd_cat_vs_dm_ppv)) %>%
  dplyr::select(ID, nd_cat_vs_dm_ppv, age, gender, bmi)

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

# Modified: Forest plot with significance markers using original group names
create_forest_plot_with_sig <- function(contrast_df, group_col, colors_map, 
                                        title, metric_name = "RHR") {
  # Add space for significance markers
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  p <- ggplot(contrast_df, aes(x = estimate, y = !!sym(group_col), color = !!sym(group_col))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
    # Add significance markers
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    scale_color_manual(
      values = colors_map,
      # Add label conversion
      labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", 
                 "Diabetes_PPV" = "Diabetic PPV")
    ) +
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
    ) +
    # Add y-axis label conversion
    scale_y_discrete(labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", 
                                "Diabetes_PPV" = "Diabetic PPV"))
  
  return(p)
}

# Function to extract group comparisons with confidence intervals
# 修改后的函数，将比较方向从 "No_Diabetes_Cataract - Diabetes_PPV" 改为 "Diabetes_PPV - No_Diabetes_Cataract"
extract_group_comparison <- function(model, group_var, reference_group) {
  # 获取每个组和时期的边际均值
  emm <- emmeans(model, specs = c(group_var, "period"))
  
  # 创建组间对比
  contrasts <- contrast(emm, method = "revpairwise", by = "period") 
  
  # 转换为数据框
  contrast_df <- as.data.frame(contrasts)
  
  # 过滤只包含参考组的比较
  contrast_df <- contrast_df %>%
    filter(grepl(reference_group, contrast))
  
  # 添加更易读的比较名称
  contrast_df$comparison <- gsub(paste0(" - ", reference_group), "", contrast_df$contrast)
  
  # 倒转比较方向，从 "reference - comparison" 变为 "comparison - reference"
  # 这样做会使值的符号相反，显示方向与趋势图一致
  contrast_df$estimate <- -contrast_df$estimate
  contrast_df$contrast <- gsub(paste0(reference_group, " - "), "", contrast_df$contrast)
  contrast_df$contrast <- paste0(contrast_df$contrast, " - ", reference_group)
  
  # 计算95%置信区间
  contrast_df$t_crit <- qt(0.975, contrast_df$df)
  # 注意反转后上下置信区间也要互换
  contrast_df$lower.CL <- -contrast_df$estimate + contrast_df$t_crit * contrast_df$SE
  contrast_df$upper.CL <- -contrast_df$estimate - contrast_df$t_crit * contrast_df$SE
  # 再次反转回来，因为我们已经调整了estimate
  contrast_df$lower.CL <- -contrast_df$lower.CL
  contrast_df$upper.CL <- -contrast_df$upper.CL
  
  # 添加显著性标记
  contrast_df$significance <- ""
  contrast_df$significance[contrast_df$p.value < 0.05] <- "*"
  contrast_df$significance[contrast_df$p.value < 0.01] <- "**"
  contrast_df$significance[contrast_df$p.value < 0.001] <- "***"
  
  return(contrast_df)
}

# 修改后的比较森林图函数，标题和轴标签使用一致的差异方向
create_comparison_forest_plot <- function(contrast_df, colors_map, reference_group,
                                          title, metric_name = "RHR") {
  # 确保时期是一个有序因子
  contrast_df$period <- factor(contrast_df$period, levels = c("Pre-surgery", "Post-surgery"))
  
  # 添加空间用于显示显著性标记
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  # 转换参考组名称用于显示在标题中
  reference_group_display <- ifelse(reference_group == "No_Diabetes_Cataract", 
                                    "Non-Diabetic Cataract", reference_group)
  
  # 创建森林图
  p <- ggplot(contrast_df, aes(x = estimate, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3, color = colors_map["Diabetes_PPV"]) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, 
                   color = colors_map["Diabetes_PPV"]) +
    # 添加显著性标记
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    facet_wrap(~ period, ncol = 1) +
    labs(
      title = title,
      subtitle = paste0("Positive value indicates higher ", metric_name, 
                        " compared to ", reference_group_display),
      x = paste0("Diabetic PPV - ", reference_group_display, " ", metric_name, " (mean & CI 95%)"),
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
    ) +
    # 添加y轴标签转换
    scale_y_discrete(labels = c("Diabetes_PPV" = "Diabetic PPV"))
  
  return(p)
}

#########################
# Step 6: Create output directory
#########################
dir.create("3_data_analysis/4_correlation_analysis/RHR/RHR_mixed_effects_models/surgery_type", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/4_correlation_analysis/RHR/RHR_mixed_effects_models/surgery_type")

#########################
# Step 7: Analyze each RHR metric
#########################

# Function to analyze and visualize a specific RHR metric by surgery type
analyze_and_visualize_rhr_by_surgery <- function(stat_type) {
  # Merge RHR data with surgery group information
  combined_data <- filtered_rhr_data %>%
    left_join(nd_cat_vs_dm_ppv_groups, by = c("subject_id" = "ID")) %>%
    filter(!is.na(nd_cat_vs_dm_ppv))
  
  # Convert to long format
  demo_cols <- c("subject_id", "surgery_date", "nd_cat_vs_dm_ppv", "age", "gender", "bmi")
  rhr_cols <- setdiff(names(combined_data), demo_cols)
  
  long_data <- combined_data %>%
    pivot_longer(
      cols = all_of(rhr_cols),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      day = as.numeric(str_extract(variable, "-?\\d+")),
      stat_type_extracted = str_extract(variable, "(mean|min|max|median|sd|cv|iqr|skew|kurt)"),
      rhr_type = str_extract(variable, "rhr_\\d+"),
      period = ifelse(day < 0, "Pre-surgery", "Post-surgery"),
      gender = as.factor(gender),
      bmi = as.numeric(bmi),
      age = as.numeric(age)
    ) %>%
    filter(!is.na(value), day >= -7, day <= 30, 
           stat_type_extracted == stat_type, rhr_type == "rhr_1")
  
  # Fit mixed effects model
  lme_model <- tryCatch({
    lme(value ~ nd_cat_vs_dm_ppv * period + age + gender + bmi, 
        data = long_data,
        random = ~1|subject_id,
        weights = varIdent(form = ~1|period),
        method = "REML",
        na.action = na.omit)
  }, error = function(e) {
    # Fallback to a simpler model if the original fails
    message("Error in mixed model, trying simpler version: ", e$message)
    lme(value ~ nd_cat_vs_dm_ppv * period + age + gender + bmi, 
        data = long_data,
        random = ~1|subject_id,
        method = "REML",
        na.action = na.omit)
  })
  
  # Calculate marginal means
  emm_interaction <- emmeans(lme_model, specs = ~ nd_cat_vs_dm_ppv | period)
  
  # 1. Time trend plot
  avg_data <- long_data %>%
    group_by(day, nd_cat_vs_dm_ppv) %>%
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
  p_trend <- ggplot(avg_data, aes(x = day, y = mean_value, color = nd_cat_vs_dm_ppv, group = nd_cat_vs_dm_ppv)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = nd_cat_vs_dm_ppv), alpha = 0.2, color = NA) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    scale_color_manual(
      values = surgery_colors,
      labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", "Diabetes_PPV" = "Diabetic PPV")
    ) + 
    scale_fill_manual(
      values = surgery_colors,
      labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", "Diabetes_PPV" = "Diabetic PPV")
    ) +
    labs(
      title = paste0("RHR (", stat_type, ") levels by surgery type"),
      x = "Days relative to surgery",
      y = paste0("RHR (", stat_type, ")"),
      color = "Group",
      fill = "Group"
    ) + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  # Save the time trend plot
  ggsave(paste0("surgery_type_", stat_type, "_time_trend.pdf"), p_trend, width = 10, height = 6, dpi = 300)
  
  # 2. Create box plot for group comparisons
  surgery_plot <- ggplot(long_data, aes(x = nd_cat_vs_dm_ppv, y = value, fill = nd_cat_vs_dm_ppv)) +
    # Add boxplot
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    # Add jittered points
    geom_jitter(aes(color = nd_cat_vs_dm_ppv), width = 0.2, alpha = 0.5, size = 2) +
    # Add mean point in red
    stat_summary(fun = mean, geom = "point", shape = 16, size = 5, color = "darkred") +
    # Add horizontal line at mean with label
    stat_summary(fun = mean, geom = "text", 
                 aes(label = sprintf("μ = %.2f", ..y..)),
                 hjust = -0.3, vjust = -0.5) +
    # Use custom colors
    scale_fill_manual(
      values = surgery_colors,
      labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", "Diabetes_PPV" = "Diabetic PPV")
    ) +
    scale_color_manual(
      values = surgery_colors,
      labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", "Diabetes_PPV" = "Diabetic PPV")
    ) +
    # Apply bw theme with customizations
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    # Add labels
    labs(
      title = paste0("RHR (", stat_type, ") Levels by Surgery Type"),
      x = "Surgery Type",
      y = paste0("RHR (", stat_type, ")"),
      fill = "Group",
      color = "Group"
    ) +
    # Add statistical comparison
    stat_compare_means(method = "wilcox.test",
                       label.y = max(long_data$value, na.rm = TRUE) + 0.5) +
    # Add x-axis label conversion
    scale_x_discrete(labels = c("No_Diabetes_Cataract" = "Non-Diabetic Cataract", 
                                "Diabetes_PPV" = "Diabetic PPV"))
  
  # Save the box plot
  ggsave(paste0("surgery_type_", stat_type, "_stats_comparison.pdf"), surgery_plot, width = 10, height = 6, dpi = 300)
  
  # 3. Create forest plot for pre-post differences
  surgery_contrasts <- extract_contrasts(lme_model, "nd_cat_vs_dm_ppv")
  surgery_contrasts <- add_significance(surgery_contrasts)
  
  surgery_forest_sig <- create_forest_plot_with_sig(
    surgery_contrasts, 
    "group", 
    surgery_colors,  # Use original color mapping
    paste0("RHR (", stat_type, "): Preoperative and postoperative differences by surgery type"),
    paste0("RHR (", stat_type, ")")
  )
  
  # Save the forest plot
  ggsave(paste0("surgery_type_", stat_type, "_forest_plot_period.pdf"), surgery_forest_sig, width = 10, height = 6, dpi = 300)
  
  # 4. Create between-group comparison forest plot
  surgery_comparison <- extract_group_comparison(lme_model, "nd_cat_vs_dm_ppv", "No_Diabetes_Cataract")
  
  surgery_forest <- create_comparison_forest_plot(
    surgery_comparison,  # 使用修改后的比较数据
    surgery_colors,      # 直接使用原始颜色映射
    "No_Diabetes_Cataract",  # 使用原始参考组名
    paste0("Diabetic PPV vs Non-Diabetic Cataract RHR (", stat_type, ") Differences"),
    paste0("RHR (", stat_type, ")")
  )
  
  # 保存比较森林图
  ggsave(paste0("surgery_type_", stat_type, "_comparison_forest.pdf"), surgery_forest, width = 10, height = 8, dpi = 300)
  
  # Return model and data
  return(list(model = lme_model, data = long_data))
}

#########################
# Step 8: Apply analysis to all four RHR metrics
#########################

# Analyze mean RHR
mean_results <- analyze_and_visualize_rhr_by_surgery("mean")
print("Completed mean RHR analysis and visualization")

# Analyze min RHR
min_results <- analyze_and_visualize_rhr_by_surgery("min")
print("Completed min RHR analysis and visualization")

# Analyze max RHR
max_results <- analyze_and_visualize_rhr_by_surgery("max")
print("Completed max RHR analysis and visualization")

# Analyze sd (standard deviation) of RHR
sd_results <- analyze_and_visualize_rhr_by_surgery("sd")
