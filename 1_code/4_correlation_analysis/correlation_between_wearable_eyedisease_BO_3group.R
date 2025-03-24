library(tidyverse)
library(lme4)      # For mixed effects models
library(nlme)      # For extended mixed models capabilities
library(emmeans)   # For marginal means analysis
library(car)       # For Anova function
library(ggstatsplot)
source('1_code/100_tools.R')

# Clean environment variables
rm(list = ls())
setwd(get_project_wd())

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
# Step 4: Create group information
#########################

group_info <- baseline_info %>%
  mutate(
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  dplyr::select(ID, dm_status, surgery_type, age, gender, bmi)

# Create diabetes status groups
dm_groups <- group_info %>%
  filter(!is.na(dm_status)) %>%
  mutate(dm_group = dm_status) %>%
  dplyr::select(ID, dm_group, age, gender, bmi)

# Create non-diabetes cataract vs diabetes PPV groups
nd_cat_vs_dm_ppv_groups <- group_info %>%
  mutate(
    nd_cat_vs_dm_ppv = case_when(
      dm_status == "No Diabetes" & surgery_type == 0 ~ "No_Diabetes_Cataract",
      dm_status == "Diabetes" & surgery_type == 1 ~ "Diabetes_PPV",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(nd_cat_vs_dm_ppv)) %>%
  dplyr::select(ID, nd_cat_vs_dm_ppv, age, gender, bmi)

# Create three original groups
original_groups <- group_info %>%
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

# Store all group information
group_data <- list(
  dm_groups = dm_groups,
  nd_cat_vs_dm_ppv_groups = nd_cat_vs_dm_ppv_groups,
  original_groups = original_groups
)

##pre/post period forest plot function
# 对比提取函数，添加置信区间
extract_contrasts <- function(model, group_var) {
  # 为每个组和时期获取估计边际均值
  emm <- emmeans(model, specs = as.formula(paste0("~ period | ", group_var)))
  
  # 创建对比（术后 - 术前）
  contrasts <- contrast(emm, method = "revpairwise") 
  
  # 转换为数据框
  contrast_df <- as.data.frame(contrasts)
  
  # 添加组信息
  contrast_df$group <- contrast_df[[group_var]]
  
  # 手动计算95%置信区间
  # 使用t分布的临界值（通常为1.96，但在这里使用更精确的基于df的值）
  contrast_df$t_crit <- qt(0.975, contrast_df$df)
  contrast_df$lower.CL <- contrast_df$estimate - contrast_df$t_crit * contrast_df$SE
  contrast_df$upper.CL <- contrast_df$estimate + contrast_df$t_crit * contrast_df$SE
  
  return(contrast_df)
}

# 添加p值到对比结果的函数
add_significance <- function(contrast_df) {
  contrast_df$significance <- ""
  contrast_df$significance[contrast_df$p.value < 0.05] <- "*"
  contrast_df$significance[contrast_df$p.value < 0.01] <- "**"
  contrast_df$significance[contrast_df$p.value < 0.001] <- "***"
  return(contrast_df)
}

# 创建带有显著性标记的森林图
create_forest_plot_with_sig <- function(contrast_df, group_col, colors_map, 
                                        title, metric_name = "BO（mean）") {
  # 添加一点空间用于显著性标记
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  p <- ggplot(contrast_df, aes(x = estimate, y = group, color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
    # 添加显著性标记
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    scale_color_manual(values = colors_map) +
    labs(
      title = title,
      subtitle = "ositive value indicates an increase after surgery",
      x = paste0(metric_name, " post - pre (mean & CI 95%)"),
      y = "",
      caption = "Adjusted fo age gender, and BMI\n* p<0.05, ** p<0.01, *** p<0.001"
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
# Step 5: Create output directory
#########################
dir.create("3_data_analysis/4_correlation_analysis/BO/BO_mixed_effects_models", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/4_correlation_analysis/BO/BO_mixed_effects_models")

#########################
# Step 6: Analyze diabetes groups - mean metrics
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

# Create visualizations
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
    subtitle = "Mixed-effects model adjusted for age, gender and BMI",
    x = "Days relative to surgery",
    y = "Blood oxygen (mean)"
  ) + theme_bw() + theme(legend.position = "bottom")
p_trend
ggsave("dm_group_mean_time_trend.pdf", p_trend, width = 10, height = 6, dpi = 300)

# 2. Pre/Post-surgery comparison
emm_data <- as.data.frame(emm_interaction)
names(emm_data)[1] <- "dm_group"


# 创建箱型图
dm_plot <- ggbetweenstats(
  data = mean_data,
  x = dm_group,
  y = value,
  type = "parametric",  # 参数统计
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  p.adjust.method = "none",
  point.args = list(
    alpha = 0.7,
    size = 2,
    position = position_jitter(width = 0.1)
  ),
  centrality.plotting = TRUE,
  centrality.type = "parametric"  # 使用均值
) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Blood oxygen (mean)") +
  scale_color_manual(values = diabetes_colors) +
  scale_fill_manual(values = diabetes_colors)

dm_plot 

# 保存图表
ggsave("dm_group_mean_stats_comparison.pdf", dm_plot, width = 10, height = 6, dpi = 300)


#####forest plot
dm_contrasts <- extract_contrasts(lme_model, "dm_group")
# 添加显著性标记
dm_contrasts <- add_significance(dm_contrasts)
# 创建森林图
dm_forest_sig <- create_forest_plot_with_sig(
  dm_contrasts, 
  "dm_group", 
  diabetes_colors,
  "BO level: Preoperative and postoperative differences by diabetes status"
)

dm_forest_sig

# 保存图表
ggsave( "dm_group_forest_plot_period.pdf", dm_forest_sig, width = 10, height = 6, dpi = 300)

#########################
# Step 7: Analyze non-diabetes cataract vs diabetes PPV groups
#########################

# Merge data
combined_data <- filtered_bo_data %>%
  left_join(nd_cat_vs_dm_ppv_groups, by = c("subject_id" = "ID")) %>%
  filter(!is.na(nd_cat_vs_dm_ppv))

# Convert to long format
demo_cols <- c("subject_id", "surgery_date", "nd_cat_vs_dm_ppv", "age", "gender", "bmi")
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

# Only analyze mean statistics
mean_data <- long_data %>% filter(stat_type == "mean")

# Fit mixed effects model
lme_model <- lme(value ~ nd_cat_vs_dm_ppv * period + age + gender + bmi, 
                 data = mean_data,
                 random = ~1|subject_id,
                 weights = varIdent(form = ~1|period),
                 method = "REML",
                 na.action = na.omit)

# Calculate marginal means
emm_interaction <- emmeans(lme_model, specs = ~ nd_cat_vs_dm_ppv | period)

# Create visualizations
# 1. Time trend plot
avg_data <- mean_data %>%
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
  geom_line(size = 1) + geom_point(size = 2) +
  scale_color_manual(values = surgery_colors) + 
  scale_fill_manual(values = surgery_colors) +
  labs(
    title = "Mean blood oxygen by surgery group",
    subtitle = "No Diabetes Cataract vs Diabetes PPV",
    x = "Days relative to surgery",
    y = "Blood oxygen (mean)"
  ) + theme_bw() + theme(legend.position = "bottom")

p_trend

ggsave("nd_cat_vs_dm_ppv_mean_time_trend.pdf", p_trend, width = 10, height = 6, dpi = 300)

# 2. Pre/Post-surgery comparison
emm_data <- as.data.frame(emm_interaction)
names(emm_data)[1] <- "nd_cat_vs_dm_ppv"

nd_vs_dm_plot <- ggbetweenstats(
  data = mean_data,  # 确保这里使用的是包含nd_cat_vs_dm_ppv列的数据集
  x = nd_cat_vs_dm_ppv,
  y = value,
  type = "parametric",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  p.adjust.method = "none",
  point.args = list(
    alpha = 0.7,
    size = 2,
    position = position_jitter(width = 0.1)
  ),
  centrality.plotting = TRUE,
  centrality.type = "parametric"
) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Blood oxygen (mean)") +
  scale_color_manual(values = surgery_colors) +
  scale_fill_manual(values = surgery_colors)
nd_vs_dm_plot

# 保存图表
ggsave("nd_cat_vs_dm_ppv_mean_stats_comparison.pdf", nd_vs_dm_plot, width = 10, height = 6, dpi = 300)

#####forest_plot
cat_ppv_contrasts <- extract_contrasts(lme_model, "nd_cat_vs_dm_ppv")
cat_ppv_contrasts <- add_significance(cat_ppv_contrasts)
cat_ppv_forest_sig <- create_forest_plot_with_sig(
  cat_ppv_contrasts, 
  "nd_cat_vs_dm_ppv", 
  surgery_colors,
  "BO level: Preoperative and postoperative differences by type of surgery"
)

cat_ppv_forest_sig 
ggsave("cat_ppv_forest_plot.pdf", cat_ppv_forest_sig, width = 10, height = 6, dpi = 300)

#########################
# Step 8: Analyze three original groups
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

# Only analyze mean statistics
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

# Create visualizations
# 1. Time trend plot
avg_data <- mean_data %>%
  group_by(day, orig_group) %>%
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
p_trend <- ggplot(avg_data, aes(x = day, y = mean_value, color = orig_group, group = orig_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = orig_group), alpha = 0.2, color = NA) +
  geom_line(size = 1) + geom_point(size = 2) +
  scale_color_manual(values = original_colors) + 
  scale_fill_manual(values = original_colors) +
  labs(
    title = "Mean blood oxygen by original groups",
    subtitle = "Mixed-effects model adjusted for age, gender and BMI",
    x = "Days relative to surgery",
    y = "Blood oxygen (mean)"
  ) + theme_bw() + theme(legend.position = "bottom")

p_trend 

ggsave("orig_group_mean_time_trend.pdf", p_trend, width = 10, height = 6, dpi = 300)

# 2. Pre/Post-surgery comparison
emm_data <- as.data.frame(emm_interaction)
names(emm_data)[1] <- "orig_group"

orig_plot <- ggbetweenstats(
  data = mean_data,  # 确保这里使用的是包含orig_group列的数据集
  x = orig_group,
  y = value,
  type = "parametric",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  p.adjust.method = "none",
  point.args = list(
    alpha = 0.7,
    size = 2,
    position = position_jitter(width = 0.1)
  ),
  centrality.plotting = TRUE,
  centrality.type = "parametric"
) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Blood oxygen (mean)") +
  scale_color_manual(values = original_colors) +
  scale_fill_manual(values = original_colors)
orig_plot 
# 保存图表
ggsave("orig_group_mean_stats_comparison.pdf", orig_plot, width = 10, height = 6, dpi = 300)









# 通用对比提取函数
extract_group_comparison <- function(model, group_var, reference_level) {
  # 为每个时期计算组间对比
  # 1. 术前
  emm_formula <- as.formula(paste0("~", group_var, " | period"))
  emm_pre <- emmeans(model, specs = emm_formula, at = list(period = "Pre-surgery"))
  
  # 创建对比
  contrasts_pre <- contrast(emm_pre, method = "trt.vs.ctrl", ref = reference_level)
  contrast_df_pre <- as.data.frame(contrasts_pre)
  contrast_df_pre$period <- "Pre-surgery"
  
  # 2. 术后
  emm_post <- emmeans(model, specs = emm_formula, at = list(period = "Post-surgery"))
  contrasts_post <- contrast(emm_post, method = "trt.vs.ctrl", ref = reference_level)
  contrast_df_post <- as.data.frame(contrasts_post)
  contrast_df_post$period <- "Post-surgery"
  
  # 3. 合并两个时期的结果
  contrast_df <- rbind(contrast_df_pre, contrast_df_post)
  
  # 手动计算95%置信区间
  contrast_df$t_crit <- qt(0.975, contrast_df$df)
  contrast_df$lower.CL <- contrast_df$estimate - contrast_df$t_crit * contrast_df$SE
  contrast_df$upper.CL <- contrast_df$estimate + contrast_df$t_crit * contrast_df$SE
  
  # 添加显著性标记
  contrast_df$significance <- ""
  contrast_df$significance[contrast_df$p.value < 0.05] <- "*"
  contrast_df$significance[contrast_df$p.value < 0.01] <- "**"
  contrast_df$significance[contrast_df$p.value < 0.001] <- "***"
  
  # 提取组名（去掉变量名前缀）
  contrast_df$comparison <- gsub(paste0(group_var, " "), "", contrast_df$contrast)
  
  return(contrast_df)
}

# 创建森林图函数
create_comparison_forest_plot <- function(contrast_df, colors_map, reference_group,
                                          title, metric_name = "血氧（平均值）") {
  # 确保period是有序因子
  contrast_df$period <- factor(contrast_df$period, levels = c("Pre-surgery", "Post-surgery"))
  
  # 添加一点空间用于显著性标记
  max_upper <- max(contrast_df$upper.CL, na.rm = TRUE)
  sig_position <- max_upper + abs(max_upper) * 0.1
  
  # 创建森林图
  p <- ggplot(contrast_df, aes(x = estimate, y = comparison, color = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
    # 添加显著性标记
    geom_text(aes(x = sig_position, label = significance), 
              hjust = 0, size = 5) +
    scale_color_manual(values = colors_map) +
    facet_wrap(~ period, ncol = 1) +
    labs(
      title = title,
      subtitle = paste0("正值表示相比", reference_group, "血氧水平高"),
      x = paste0("各组 - ", reference_group, " ", metric_name, " (平均值 & CI 95%)"),
      y = "",
      caption = "混合效应模型已调整年龄、性别和BMI\n* p<0.05, ** p<0.01, *** p<0.001"
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

########################################################################
# 1. 糖尿病分组对比 (dm_group)
########################################################################

# 提取对比
dm_comparison <- extract_group_comparison(lme_model, "dm_group", "No Diabetes")

# 创建森林图
dm_forest <- create_comparison_forest_plot(
  dm_comparison,
  diabetes_colors,
  "No Diabetes",
  "糖尿病组与非糖尿病组血氧水平差异"
)
dm_forest

# 保存图表
ggsave(file.path(plot_dir, "dm_group_comparison_forest.png"), dm_forest, width = 10, height = 8, dpi = 300)

########################################################################
# 2. 手术类型分组对比 (nd_cat_vs_dm_ppv)
########################################################################

# 提取对比
surgery_comparison <- extract_group_comparison(lme_model, "nd_cat_vs_dm_ppv", "No_Diabetes_Cataract")

# 创建森林图
surgery_forest <- create_comparison_forest_plot(
  surgery_comparison,
  surgery_colors,
  "No_Diabetes_Cataract",
  "糖尿病PPV组与非糖尿病白内障组血氧水平差异"
)
surgery_forest
# 保存图表
ggsave(file.path(plot_dir, "surgery_type_comparison_forest.png"), surgery_forest, width = 10, height = 8, dpi = 300)

########################################################################
# 3. 原始三组分组对比 (orig_group)
########################################################################

# 提取对比
orig_comparison <- extract_group_comparison(lme_model, "orig_group", "No_Diabetes_Surgery0")

# 创建森林图
orig_forest <- create_comparison_forest_plot(
  orig_comparison,
  original_colors,
  "No_Diabetes_Surgery0",
  "各组与无糖尿病Surgery0组血氧水平差异"
)
orig_forest
# 保存图表
ggsave(file.path(plot_dir, "orig_group_comparison_forest.png"), orig_forest, width = 10, height = 8, dpi = 300)
