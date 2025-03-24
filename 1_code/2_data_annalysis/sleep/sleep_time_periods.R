library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 用于计算偏度和峰度
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


######## 读取数据
load(
  "3_data_analysis/1_data_preparation/wearable_data/7_sleep/sleep_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############创建输出目录
dir.create("3_data_analysis/2_data_analysis/sleep/sleep_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/sleep/sleep_time_periods")


###############睡眠时间窗口处理
# 处理基线信息
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# 提取睡眠数据的函数
extract_sleep_data <- function(sleep_data) {
  # 获取所有睡眠指标
  sleep_metrics <- sleep_data@expression_data %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
  
  # 添加样本信息
  sleep_metrics <- sleep_metrics %>%
    mutate(
      sample_id = rownames(.),
      timestamp = sleep_data@sample_info$measure_time,
      subject_id = sleep_data@sample_info$subject_id,
      sleep_start = sleep_data@sample_info$sleep_start_time,
      sleep_end = sleep_data@sample_info$sleep_end_time
    )
  
  # 重命名睡眠指标列
  colnames(sleep_metrics)[1:7] <- c(
    "light_sleep",
    "deep_sleep",
    "dream_sleep",
    "awake",
    "total_sleep",
    "daytime_sleep",
    "sleep_score"
  )
  
  # 添加睡眠时长计算
  sleep_metrics %>%
    mutate(
      sleep_duration = as.numeric(difftime(sleep_end, sleep_start, units = "mins"))
    ) %>%
    arrange(subject_id, sleep_start)
}

# 提取睡眠数据
sleep_data_processed <- extract_sleep_data(sleep_data)

# 更新后的睡眠数据时间窗口处理函数 - 添加新的时间窗口
process_sleep_time_periods <- function(sleep_data_processed, baseline_info_processed) {
  # 计算基本时间信息
  base_data <- sleep_data_processed %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      # 使用sleep_start计算相对于手术时间的小时数
      hours_to_surgery = as.numeric(difftime(sleep_start, surgery_time_1, units = "hours"))
    )
  
  # 定义时间窗口统计函数
  calculate_time_period_stats <- function(data, min_hours, max_hours, period_name) {
    data %>%
      filter(hours_to_surgery >= min_hours & hours_to_surgery < max_hours) %>%
      group_by(subject_id) %>%
      summarise(
        light_sleep_total = sum(light_sleep, na.rm = TRUE),
        deep_sleep_total = sum(deep_sleep, na.rm = TRUE),
        dream_sleep_total = sum(dream_sleep, na.rm = TRUE),
        awake_total = sum(awake, na.rm = TRUE),
        total_sleep_total = sum(total_sleep, na.rm = TRUE),
        daytime_sleep_total = sum(daytime_sleep, na.rm = TRUE),
        sleep_score_mean = mean(sleep_score, na.rm = TRUE),
        sleep_duration_mean = mean(sleep_duration, na.rm = TRUE),
        sleep_duration_median = median(sleep_duration, na.rm = TRUE),
        sleep_duration_min = min(sleep_duration, na.rm = TRUE),
        sleep_duration_max = max(sleep_duration, na.rm = TRUE),
        sleep_duration_sd = sd(sleep_duration, na.rm = TRUE),
        sleep_duration_cv = (sd(sleep_duration, na.rm = TRUE) / mean(sleep_duration, na.rm = TRUE)) * 100,
        n_sleep_sessions = n(),
        .groups = "drop"
      ) %>%
      mutate(time_period = period_name)
  }
  
  # 处理术前时间窗口
  pre_3d <- calculate_time_period_stats(base_data, -72, 0, "pre_surgery_3d")
  pre_3d_to_7d <- calculate_time_period_stats(base_data, -168, -72, "pre_surgery_3d_to_7d")
  pre_7d_all <- calculate_time_period_stats(base_data, -168, 0, "pre_surgery_7d_all")
  pre_all <- calculate_time_period_stats(base_data, -720, 0, "pre_surgery_all")
  
  # 处理术后时间窗口 - 添加新的时间窗口
  # 术后1-3天
  post_1to3d <- calculate_time_period_stats(base_data, 24, 72, "post_surgery_1to3d")
  
  # 术后4-6天
  post_4to6d <- calculate_time_period_stats(base_data, 72, 144, "post_surgery_4to6d")
  
  # 术后6天 (0-144小时)
  post_6d <- calculate_time_period_stats(base_data, 0, 144, "post_surgery_6d")
  
  # 术后1周 (0-168小时) - 保留原有窗口
  post_1w <- calculate_time_period_stats(base_data, 0, 168, "post_surgery_1w")
  
  # 术后7天 (原有窗口)
  post_7d <- calculate_time_period_stats(base_data, 0, 168, "post_surgery_7d")
  
  # 术后7-30天
  post_7d_to_30d <- calculate_time_period_stats(base_data, 168, 720, "post_surgery_7d_to_30d")
  
  # 术后23-30天
  post_day23_to_30 <- calculate_time_period_stats(base_data, 552, 720, "post_surgery_day23_to_30")
  
  # 术后27-30天
  post_day27_to_30 <- calculate_time_period_stats(base_data, 648, 720, "post_surgery_day27_to_30")
  
  # 术后>30天
  post_over_30d <- calculate_time_period_stats(base_data, 720, Inf, "post_surgery_over_30d")
  
  # 合并所有时间窗口结果
  bind_rows(
    pre_3d, pre_3d_to_7d, pre_7d_all, pre_all,
    post_1to3d, post_4to6d, post_6d,
    post_1w, post_7d, post_7d_to_30d,
    post_day23_to_30, post_day27_to_30,
    post_over_30d
  )
}

# 处理睡眠数据
sleep_time_periods <- process_sleep_time_periods(sleep_data_processed, baseline_info_processed)

# 为每个统计量创建宽格式数据的函数
create_wide_format <- function(sleep_time_periods) {
  # 要处理的指标列表
  metrics <- c(
    "light_sleep_total", "deep_sleep_total", "dream_sleep_total", "awake_total", 
    "total_sleep_total", "daytime_sleep_total", "sleep_score_mean", 
    "sleep_duration_mean", "sleep_duration_median", "sleep_duration_min", 
    "sleep_duration_max", "sleep_duration_sd", "sleep_duration_cv", "n_sleep_sessions"
  )
  
  # 处理每个指标
  stats_list <- list()
  
  for(metric in metrics) {
    # 对于n_sleep_sessions，直接使用
    if(metric == "n_sleep_sessions") {
      stats_list[[metric]] <- sleep_time_periods %>%
        mutate(col_name = paste0(time_period, "_", metric)) %>%
        dplyr::select(subject_id, col_name, !!sym(metric)) %>%
        pivot_wider(
          names_from = col_name,
          values_from = !!sym(metric)
        )
    } else {
      # 对于其他指标，使用sleep_前缀
      stats_list[[metric]] <- sleep_time_periods %>%
        mutate(sleep_col_name = paste0(time_period, "_sleep_", gsub("sleep_|_sleep", "", metric))) %>%
        dplyr::select(subject_id, sleep_col_name, !!sym(metric)) %>%
        pivot_wider(
          names_from = sleep_col_name,
          values_from = !!sym(metric)
        )
    }
  }
  
  # 合并所有指标的数据
  result <- stats_list[[1]]
  for(i in 2:length(stats_list)) {
    result <- result %>%
      left_join(stats_list[[i]], by = "subject_id")
  }
  
  return(result)
}

# 创建包含手术日期的最终结果
create_sleep_time_period_results <- function(sleep_time_periods, baseline_info_processed) {
  # 创建宽格式数据
  sleep_wide <- create_wide_format(sleep_time_periods)
  
  # 创建最终结果，包含手术日期
  tibble(subject_id = unique(sleep_time_periods$subject_id)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(sleep_wide, by = "subject_id") %>%
    # 重新排序列，将subject_id和surgery_date放在开头
    dplyr::select(subject_id, surgery_date, everything())
}

# 创建宽格式结果，使用新的列名格式
sleep_time_period_results <- create_sleep_time_period_results(sleep_time_periods, baseline_info_processed)

# 验证列名格式
cat("新的睡眠数据列名格式示例:\n")
print(names(sleep_time_period_results)[1:20])

# 检查新添加时间窗口的列数量
cat("\n检查新添加时间窗口的列数:\n")
cat("术后1-3天(post_surgery_1to3d)列数: ", 
    sum(grepl("post_surgery_1to3d", names(sleep_time_period_results))), "\n")
cat("术后4-6天(post_surgery_4to6d)列数: ", 
    sum(grepl("post_surgery_4to6d", names(sleep_time_period_results))), "\n")
cat("术后6天(post_surgery_6d)列数: ", 
    sum(grepl("post_surgery_6d", names(sleep_time_period_results))), "\n")

# 保存结果
save(sleep_time_period_results, file = "sleep_time_period_results.rda", compress = "xz")
save(sleep_time_periods, file = "sleep_time_periods_long.rda", compress = "xz")

# 计算汇总统计数据
calculate_sleep_summary_stats <- function(sleep_time_periods) {
  sleep_time_periods %>%
    group_by(time_period) %>%
    summarise(
      light_sleep_mean = mean(light_sleep_total, na.rm = TRUE) / 60,
      deep_sleep_mean = mean(deep_sleep_total, na.rm = TRUE) / 60,
      dream_sleep_mean = mean(dream_sleep_total, na.rm = TRUE) / 60,
      awake_mean = mean(awake_total, na.rm = TRUE) / 60,
      total_sleep_mean = mean(total_sleep_total, na.rm = TRUE) / 60,
      sleep_score_mean = mean(sleep_score_mean, na.rm = TRUE),
      sleep_duration_mean_hrs = mean(sleep_duration_mean, na.rm = TRUE) / 60,
      n_subjects = n(),
      .groups = "drop"
    )
}

# 计算汇总统计量
sleep_time_period_summary_stats <- calculate_sleep_summary_stats(sleep_time_periods)
save(sleep_time_period_summary_stats, file = "sleep_time_period_summary_stats.rda", compress = "xz")

# 打印最终数据验证
cat("\n最终数据验证:\n")
cat("唯一受试者数量: ", length(unique(sleep_time_periods$subject_id)), "\n")
cat("时间窗口数量: ", length(unique(sleep_time_periods$time_period)), "\n")
cat("宽格式中的列数: ", ncol(sleep_time_period_results), "\n")

###############绘制睡眠模式图
# 更新绘图函数，添加新的时间窗口

# 按时间窗口绘制总睡眠时间
plot_total_sleep_by_period <- function(sleep_time_periods) {
  # 设置时间窗口顺序，包含新增窗口
  plot_data <- sleep_time_periods %>%
    mutate(time_period_factor = factor(
      time_period,
      levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
                 "post_surgery_1to3d", "post_surgery_4to6d", "post_surgery_6d",
                 "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
                 "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
      labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
                 "Post 1-3d", "Post 4-6d", "Post 6d",
                 "Post 1w", "Post 7d", "Post 7d to 30d", 
                 "Post day 23-30", "Post day 27-30", "Post >30d")
    )) %>%
    group_by(time_period_factor) %>%
    summarise(
      mean_total_sleep = mean(total_sleep_total, na.rm = TRUE),
      sd_total_sleep = sd(total_sleep_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 转换为小时以便更好地可视化
  plot_data <- plot_data %>%
    mutate(
      mean_total_sleep_hours = mean_total_sleep / 60,
      sd_total_sleep_hours = sd_total_sleep / 60
    )
  
  # 创建图表
  ggplot(plot_data, aes(x = time_period_factor, y = mean_total_sleep_hours)) +
    geom_bar(stat = "identity", fill = "#7CB9E8", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(ymin = mean_total_sleep_hours - sd_total_sleep_hours/60, 
                      ymax = mean_total_sleep_hours + sd_total_sleep_hours/60),
                  width = 0.2, alpha = 0.7) +
    labs(title = "Total Sleep Time by Time Period",
         x = "Time Period",
         y = "Total Sleep Time (hours)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# 睡眠组成图 - 更新以包含新时间窗口
plot_sleep_composition <- function(sleep_time_periods) {
  # 设置时间窗口顺序，包含新增窗口
  plot_data <- sleep_time_periods %>%
    mutate(time_period_factor = factor(
      time_period,
      levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
                 "post_surgery_1to3d", "post_surgery_4to6d", "post_surgery_6d",
                 "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
                 "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
      labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
                 "Post 1-3d", "Post 4-6d", "Post 6d",
                 "Post 1w", "Post 7d", "Post 7d to 30d", 
                 "Post day 23-30", "Post day 27-30", "Post >30d")
    )) %>%
    group_by(time_period_factor) %>%
    summarise(
      light_sleep = mean(light_sleep_total, na.rm = TRUE) / 60,
      deep_sleep = mean(deep_sleep_total, na.rm = TRUE) / 60,
      dream_sleep = mean(dream_sleep_total, na.rm = TRUE) / 60,
      awake = mean(awake_total, na.rm = TRUE) / 60,
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(light_sleep, deep_sleep, dream_sleep, awake),
      names_to = "sleep_type",
      values_to = "hours"
    ) %>%
    mutate(sleep_type = factor(sleep_type, 
                               levels = c("light_sleep", "deep_sleep", "dream_sleep", "awake"),
                               labels = c("Light Sleep", "Deep Sleep", "REM Sleep", "Awake")
    ))
  
  # 创建图表
  ggplot(plot_data, aes(x = time_period_factor, y = hours, fill = sleep_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    labs(title = "Sleep Composition by Time Period",
         x = "Time Period",
         y = "Duration (hours)",
         fill = "Sleep Phase") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}

# 睡眠评分图 - 更新以包含新时间窗口
plot_sleep_score <- function(sleep_time_periods) {
  # 设置时间窗口顺序，包含新增窗口
  plot_data <- sleep_time_periods %>%
    mutate(time_period_factor = factor(
      time_period,
      levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
                 "post_surgery_1to3d", "post_surgery_4to6d", "post_surgery_6d",
                 "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
                 "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
      labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
                 "Post 1-3d", "Post 4-6d", "Post 6d",
                 "Post 1w", "Post 7d", "Post 7d to 30d", 
                 "Post day 23-30", "Post day 27-30", "Post >30d")
    )) %>%
    group_by(time_period_factor) %>%
    summarise(
      mean_score = mean(sleep_score_mean, na.rm = TRUE),
      sd_score = sd(sleep_score_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 创建图表
  ggplot(plot_data, aes(x = time_period_factor, y = mean_score)) +
    geom_bar(stat = "identity", fill = "#4682B4", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
                  width = 0.2, alpha = 0.7) +
    labs(title = "Average Sleep Score by Time Period",
         x = "Time Period",
         y = "Sleep Score") +
    scale_y_continuous(limits = c(0, 100)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# 创建并保存图表
p_total_sleep <- plot_total_sleep_by_period(sleep_time_periods)
p_composition <- plot_sleep_composition(sleep_time_periods)
p_score <- plot_sleep_score(sleep_time_periods)

# 显示图表
print(p_total_sleep)
print(p_composition)
print(p_score)

# 保存图表
ggsave("total_sleep_by_time_period.pdf", p_total_sleep, width = 10, height = 6, dpi = 300)
ggsave("sleep_composition_by_time_period.pdf", p_composition, width = 10, height = 6, dpi = 300)
ggsave("sleep_score_by_time_period.pdf", p_score, width = 10, height = 6, dpi = 300)

# 打印汇总统计数据
print("按时间窗口的睡眠汇总统计:")
print(sleep_time_period_summary_stats)