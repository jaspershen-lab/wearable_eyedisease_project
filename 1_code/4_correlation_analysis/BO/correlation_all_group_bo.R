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

# 修改对比提取函数 - 处理标签问题
extract_group_comparison_with_control <- function(model, group_var, control_group, comparison_groups) {
  # 获取边际均值
  emm <- emmeans(model, specs = c(group_var, "period"))
  
  # 创建每个时期内的组间对比
  contrasts <- contrast(emm, method = "trt.vs.ctrl", by = "period", 
                        ref = which(levels(as.factor(emmeans:::get_emm_option("levels")[[group_var]])) == control_group))
  
  # 转换为数据框
  contrast_df <- as.data.frame(contrasts)
  
  # 过滤只包含目标比较组的对比
  contrast_df <- contrast_df %>%
    filter(contrast %in% paste0(comparison_groups, " - ", control_group))
  
  # 添加易读的对比名称
  contrast_df$comparison <- gsub(paste0(" - ", control_group), "", contrast_df$contrast)
  
  # 计算95%置信区间
  contrast_df$t_crit <- qt(0.975, contrast_df$df)
  contrast_df$lower.CL <- contrast_df$estimate - contrast_df$t_crit * contrast_df$SE
  contrast_df$upper.CL <- contrast_df$estimate + contrast_df$t_crit * contrast_df$SE
  
  # 添加显著性标记
  contrast_df$significance <- ""
  contrast_df$significance[contrast_df$p.value < 0.05] <- "*"
  contrast_df$significance[contrast_df$p.value < 0.01] <- "**"
  contrast_df$significance[contrast_df$p.value < 0.001] <- "***"
  
  return(contrast_df)
}

# 简化版的森林图绘制函数，只显示与对照组的比较
create_diabetes_forest_plot <- function(contrast_df, colors_map, control_group,
                                        title, metric_name = "Blood Oxygen (mean)") {
  # 确保时期是有序因子
  contrast_df$period <- factor(contrast_df$period, levels = c("Pre-surgery", "Post-surgery"))
  
  # 调试打印
  print("对比数据框内容:")
  print(contrast_df)
  
  # 为显著性标记添加空间
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  min_lower <- min(contrast_df$lower.CL, na.rm = TRUE)
  
  # 计算x轴范围
  x_range <- max(abs(max_upper), abs(min_lower)) * 1.2
  x_min <- min(-x_range, min_lower * 1.2)
  x_max <- max(x_range, max_upper * 1.2)
  
  # 为显著性标记留出空间
  sig_position <- max_upper + (x_max - max_upper) * 0.7
  
  # 创建森林图
  p <- ggplot(contrast_df, aes(x = estimate, y = comparison, color = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, size = 1) +
    # 添加显著性标记
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 6, color = "black") +
    facet_wrap(~ period, ncol = 1) +
    # 设置颜色映射
    scale_color_manual(values = colors_map) +
    # 设置适当的x轴范围
    scale_x_continuous(limits = c(x_min, x_max)) +
    labs(
      title = title,
      subtitle = paste0("Positive value indicates higher BO level compared to ", control_group),
      x = paste0("Diabetes Groups - ", control_group, " ", metric_name, " (mean & CI 95%)"),
      y = "",
      caption = "Mixed-effects model adjusted for age, gender, and BMI\n* p<0.05, ** p<0.01, *** p<0.001"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "lightgray")
    )
  
  return(p)
}

#########################
# Step 6: Create output directory
#########################
dir.create("3_data_analysis/4_correlation_analysis/BO/BO_mixed_effects_models/all_group", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/4_correlation_analysis/BO/BO_mixed_effects_models/all_group")

#########################
# Step 7: Analyze diabetes groups - mean metrics
#########################
# Merge data
combined_data <- filtered_bo_data %>%
  left_join(original_groups, by = c("subject_id" = "ID")) %>%
  filter(!is.na(orig_group))

# Convert to long format
demo_cols <- c("subject_id", "surgery_date", "orig_group", "age", "gender", "bmi")
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
lme_model <- lme(value ~ orig_group * period + age + gender + bmi, 
                 data = mean_data,
                 random = ~1|subject_id,
                 weights = varIdent(form = ~1|period),
                 method = "REML",
                 na.action = na.omit)

# Calculate marginal means
emm_interaction <- emmeans(lme_model, specs = ~ orig_group | period)


#####Step 8#####
# 定义对照组和比较组
control_group <- "No_Diabetes_Surgery0"
diabetes_groups <- c("Diabetes_Surgery0", "Diabetes_Surgery1")

# 提取特定的组间对比数据
dm_vs_nodm_comparison <- extract_group_comparison_multi(lme_model, "orig_group", control_group)

# 确保只保留两个糖尿病组的比较结果
dm_vs_nodm_comparison <- dm_vs_nodm_comparison %>%
  filter(comparison %in% diabetes_groups)

# 创建糖尿病组vs非糖尿病组的森林图
diabetes_vs_control_forest <- create_diabetes_forest_plot(
  dm_vs_nodm_comparison,
  group_colors,
  control_group,
  "Diabetes Groups vs Non-Diabetes Control Blood Oxygen Level Differences"
)

# 显示和保存图表
print(diabetes_vs_control_forest)
ggsave("diabetes_vs_control_forest.pdf", diabetes_vs_control_forest, width = 10, height = 8, dpi = 300)

#########################
# Step 9: Save results
#########################

# Save model and data for potential further analysis
save(lme_model, mean_data, file = "surgery_group_analysis_results.RData")
