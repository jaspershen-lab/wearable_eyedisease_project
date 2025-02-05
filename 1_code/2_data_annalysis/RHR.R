library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)
library(lubridate)
library(zoo)
library(rpart)
library(randomForest)
library(caret)


###read data
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")
load("3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda")
baseline_info<-read.csv("2_data/analysis_data/baseline_info.csv")

######
dir.create("3_data_analysis/2_data_analysis/RHR/", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/RHR/")

#########
# heartrate data
hr_expression <- extract_expression_data(heart_rate_data)
hr_sample_info <- extract_sample_info(heart_rate_data)

# steps data
steps_expression <- extract_expression_data(daily_workout_details_data)
steps_sample_info <- extract_sample_info(daily_workout_details_data)

#check
hr_expression[1:1, 1:6]
steps_expression[1:5, 1:6]

# heartrate data clean 
hr_data <- hr_expression[1,] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")  

# steps data clean
steps_data <- steps_expression["steps",] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

# check clean
print("Processed heart rate data:")
head(hr_data)

print("\nProcessed steps data:")
head(steps_data)


str(hr_data)
str(steps_data)


# steps numerical
steps_data <- steps_data %>%
  mutate(steps = as.numeric(trimws(steps)))

# sample_id extract timestamp and subject_id
hr_data <- hr_data %>%
  mutate(
    timestamp = as.POSIXct(stringr::str_extract(sample_id, "\\d{4}-\\d{2}-\\d{2}\\s\\d{2}:\\d{2}:\\d{2}"),
                           format = "%Y-%m-%d %H:%M:%S"),
    subject_id = stringr::str_extract(sample_id, "^[^_]+")
  ) %>%
  arrange(timestamp)

steps_data <- steps_data %>%
  mutate(
    timestamp = as.POSIXct(stringr::str_extract(sample_id, "\\d{4}-\\d{2}-\\d{2}\\s\\d{2}:\\d{2}:\\d{2}"),
                           format = "%Y-%m-%d %H:%M:%S"),
    subject_id = stringr::str_extract(sample_id, "^[^_]+")
  ) %>%
  arrange(timestamp)

print("Processed heart rate data:")
head(hr_data)

print("\nProcessed steps data:")
head(steps_data)

str(hr_data)
str(steps_data)


# 将心率数据和步数数据合并
# 使用左连接保留所有心率数据，对于没有步数记录的时间点，假设步数为0
combined_data <- hr_data %>%
  left_join(steps_data %>% dplyr::select(timestamp, steps), by = "timestamp") %>%
  mutate(steps = replace_na(steps, 0)) %>%
  arrange(subject_id, timestamp)

print("Combined data:")
head(combined_data)

save(combined_data, 
     file = "hr_steps_combined_data.rda", 
     compress = "xz")

# 找到10分钟窗口，先找到连续steps<=1的序列
print("\nLooking at periods with steps <= 1:")
combined_data %>%
  filter(steps <= 1) %>%
  head(10)

# 计算每个steps<=1的序列的持续时间
print("\nChecking consecutive resting periods:")
combined_data %>%
  group_by(subject_id) %>%
  mutate(
    is_resting = steps <= 1,
    resting_group = cumsum(c(TRUE, diff(is_resting) != 0))  # 标记连续的resting序列
  ) %>%
  filter(is_resting) %>%
  group_by(subject_id, resting_group) %>%
  summarise(
    start_time = min(timestamp),
    end_time = max(timestamp),
    duration_mins = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins")),
    n_measurements = n()
  ) %>%
  filter(duration_mins >= 10) %>%  # 只看持续10分钟以上的窗口
  head(10)

# 计算每个静息窗口的平均心率
rhr_by_window <- combined_data %>%
  group_by(subject_id) %>%
  mutate(
    is_resting = steps <= 1,
    resting_group = cumsum(c(TRUE, diff(is_resting) != 0))
  ) %>%
  filter(is_resting) %>%
  group_by(subject_id, resting_group) %>%
  mutate(
    duration_mins = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins"))
  ) %>%
  filter(duration_mins >= 10) %>%  # 只保留10分钟以上的窗口
  summarise(
    start_time = min(timestamp),
    end_time = max(timestamp),
    duration_mins = dplyr::first(duration_mins),
    n_measurements = n(),
    mean_hr = mean(heart_rate),
    min_hr = min(heart_rate),
    max_hr = max(heart_rate),
    sd_hr = sd(heart_rate),
    .groups = "drop"
  )

print("RHR by resting windows:")
head(rhr_by_window)

save(rhr_by_window, 
     file = "rhr_by_window.rda", 
     compress = "xz")



# Convert surgery dates to POSIXct for proper comparison
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)

# Add surgery timing information to rhr_by_window
rhr_by_window_surgery <- rhr_by_window %>%
  left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
  mutate(
    days_to_surgery = as.numeric(difftime(surgery_time_1, start_time, units = "days"))
  ) %>%
  group_by(subject_id) %>%
  mutate(
    # Check if all data is before surgery
    has_future_surgery = all(days_to_surgery > 0),
    min_days = min(days_to_surgery),
    max_days = max(days_to_surgery)
  ) %>%
  mutate(
    # First categorize normally
    time_period = case_when(
      is.na(surgery_time_1) ~ "no_surgery_date",
      days_to_surgery >= 0 & days_to_surgery < 3 ~ "pre_surgery_3d",
      days_to_surgery >= 3 & days_to_surgery < 7 ~ "pre_surgery_3d_to_7d",
      days_to_surgery >= 7 ~ "pre_surgery_over_7d",
      days_to_surgery < 0 & days_to_surgery >= -7 ~ "post_surgery_1w",
      days_to_surgery < -7 & days_to_surgery >= -30 ~ "post_surgery_1m",
      days_to_surgery < -30 ~ "post_surgery_over_1m",
      TRUE ~ "other"
    ),
    # Then for future surgery cases, recategorize based on relative time
    time_period = case_when(
      has_future_surgery & max_days >= 7 ~ case_when(
        days_to_surgery <= min_days + 3 ~ "pre_surgery_3d",
        days_to_surgery <= min_days + 7 ~ "pre_surgery_3d_to_7d",
        TRUE ~ "pre_surgery_over_7d"
      ),
      TRUE ~ time_period
    )
  ) %>%
  ungroup()

head(rhr_by_window_surgery)


# Calculate RHR summaries for different time periods, ensuring continuous coverage
rhr_summary_by_period <- rhr_by_window_surgery %>%
  group_by(subject_id, time_period) %>%
  summarise(
    n_windows = n(),
    total_rest_minutes = sum(duration_mins),
    period_rhr = mean(mean_hr),
    min_rhr = min(mean_hr),
    max_rhr = max(mean_hr),
    sd_rhr = sd(mean_hr),
    earliest_date = min(start_time),
    latest_date = max(start_time),
    .groups = "drop"
  )


head(rhr_summary_by_period )



# Calculate pre-surgery summaries (combining all pre-surgery periods)
rhr_summary_pre_surgery <- rhr_by_window_surgery %>%
  filter(time_period %in% c("pre_surgery_over_7d", "pre_surgery_3d_to_7d", "pre_surgery_3d")) %>%
  group_by(subject_id) %>%
  summarise(
    time_period = "all_pre_surgery",
    n_windows = n(),
    total_rest_minutes = sum(duration_mins),
    period_rhr = mean(mean_hr),
    min_rhr = min(mean_hr),
    max_rhr = max(mean_hr),
    sd_rhr = sd(mean_hr),
    .groups = "drop"
  )

head(rhr_summary_pre_surgery)



# Combine all summaries
rhr_summary_complete <- bind_rows(
  rhr_summary_by_period,
  rhr_summary_pre_surgery
)

head(rhr_summary_complete)


# 每个受试者的总体RHR统计
rhr_summary <- rhr_by_window %>%
  group_by(subject_id) %>%
  summarise(
    n_windows = n(),
    total_rest_minutes = sum(duration_mins),
    overall_rhr = mean(mean_hr),
    min_rhr = min(mean_hr),
    max_rhr = max(mean_hr),
    sd_rhr = sd(mean_hr)
  )

print("\nOverall RHR summary by subject:")
print(rhr_summary)


# Add the overall RHR from original summary
rhr_summary_complete <- rhr_summary_complete %>%
  left_join(
    rhr_summary %>% 
      dplyr::select(subject_id, overall_rhr),
    by = "subject_id"
  )



# Create a wider format summary for easier comparison
rhr_summary_wide <- rhr_summary_complete %>%
  dplyr::select(subject_id, time_period, period_rhr) %>%
  pivot_wider(
    names_from = time_period,
    values_from = period_rhr,
    names_prefix = "rhr_"
  )




library(tidyverse)
library(lubridate)

# 假设我们已经有了rhr_by_window数据
# 添加时间处理的代码
rhr_daily_pattern <- rhr_by_window %>%
  # 提取时间信息
  mutate(
    hour = hour(start_time) + minute(start_time)/60
  ) %>%
  # 按小时分组计算统计量
  group_by(hour = floor(hour)) %>%
  summarise(
    median_rhr = median(mean_hr),
    n_measurements = n(),
    sd_rhr = sd(mean_hr),
    .groups = "drop"
  )

# 计算总体中位数，用于参考线
overall_median_rhr <- median(rhr_by_window$mean_hr)

# 创建可视化
ggplot(rhr_daily_pattern, aes(x = hour, y = median_rhr)) +
  # 添加主要线条
  geom_line(size = 1, color = "steelblue") +
  # 添加参考线
  geom_hline(yintercept = overall_median_rhr, 
             color = "red", 
             linetype = "dashed") +
  # 添加标签
  annotate("text", 
           x = 23, 
           y = overall_median_rhr, 
           label = "Median RHR", 
           vjust = -0.5,
           color = "red") +
  # 设置坐标轴
  scale_x_continuous(breaks = seq(0, 24, 4),
                     limits = c(0, 24),
                     name = "Time of day") +
  scale_y_continuous(name = "Median RHR",
                     limits = c(55, 80)) +
  # 主题设置
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80"),
    plot.title = element_text(hjust = 0.5)
  ) +
  # 添加标题
  ggtitle("Daily RHR Pattern")

# 保存图片
ggsave("daily_rhr_pattern.pdf", 
       width = 10, 
       height = 6, 
       dpi = 300)

# 输出数据用于检查
print("Daily RHR pattern data:")
print(rhr_daily_pattern)
print("\nOverall median RHR:")
print(overall_median_rhr)




# 找到连续steps<=50的序列
print("\nLooking at periods with steps <= 50:")
combined_data %>%
  filter(steps <= 50) %>%
  head(10)

# 计算每个steps<=50的序列的持续时间
print("\nChecking consecutive resting periods:")
combined_data %>%
  group_by(subject_id) %>%
  mutate(
    is_resting = steps <= 50,
    resting_group = cumsum(c(TRUE, diff(is_resting) != 0))  # 标记连续的resting序列
  ) %>%
  filter(is_resting) %>%
  group_by(subject_id, resting_group) %>%
  summarise(
    start_time = min(timestamp),
    end_time = max(timestamp),
    duration_mins = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins")),
    n_measurements = n()
  ) %>%
  filter(duration_mins >= 10) %>%  # 只看持续10分钟以上的窗口
  head(10)

# 计算每个静息窗口的平均心率
rhr_by_window_50 <- combined_data %>%
  group_by(subject_id) %>%
  mutate(
    is_resting = steps <= 50,
    resting_group = cumsum(c(TRUE, diff(is_resting) != 0))
  ) %>%
  filter(is_resting) %>%
  group_by(subject_id, resting_group) %>%
  mutate(
    duration_mins = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins"))
  ) %>%
  filter(duration_mins >= 10) %>%  # 只保留10分钟以上的窗口
  summarise(
    start_time = min(timestamp),
    end_time = max(timestamp),
    duration_mins = dplyr::first(duration_mins),
    n_measurements = n(),
    mean_hr = mean(heart_rate),
    min_hr = min(heart_rate),
    max_hr = max(heart_rate),
    sd_hr = sd(heart_rate),
    .groups = "drop"
  )

print("RHR by resting windows (steps <= 50):")
head(rhr_by_window_50)

save(rhr_by_window_50, 
     file = "rhr_by_window_50.rda", 
     compress = "xz")

#添加手术时间信息
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)

# Add surgery timing information to rhr_by_window_50
rhr_by_window_50_surgery <- rhr_by_window_50 %>%
  left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
  mutate(
    days_to_surgery = as.numeric(difftime(surgery_time_1, start_time, units = "days"))
  ) %>%
  group_by(subject_id) %>%
  mutate(
    has_future_surgery = all(days_to_surgery > 0),
    min_days = min(days_to_surgery),
    max_days = max(days_to_surgery)
  ) %>%
  mutate(
    time_period = case_when(
      is.na(surgery_time_1) ~ "no_surgery_date",
      days_to_surgery >= 0 & days_to_surgery < 3 ~ "pre_surgery_3d",
      days_to_surgery >= 3 & days_to_surgery < 7 ~ "pre_surgery_3d_to_7d",
      days_to_surgery >= 7 ~ "pre_surgery_over_7d",
      days_to_surgery < 0 & days_to_surgery >= -7 ~ "post_surgery_1w",
      days_to_surgery < -7 & days_to_surgery >= -30 ~ "post_surgery_1m",
      days_to_surgery < -30 ~ "post_surgery_over_1m",
      TRUE ~ "other"
    ),
    time_period = case_when(
      has_future_surgery & max_days >= 7 ~ case_when(
        days_to_surgery <= min_days + 3 ~ "pre_surgery_3d",
        days_to_surgery <= min_days + 7 ~ "pre_surgery_3d_to_7d",
        TRUE ~ "pre_surgery_over_7d"
      ),
      TRUE ~ time_period
    )
  ) %>%
  ungroup()

# 计算不同时期的RHR摘要
rhr_summary_by_period_50 <- rhr_by_window_50_surgery %>%
  group_by(subject_id, time_period) %>%
  summarise(
    n_windows = n(),
    total_rest_minutes = sum(duration_mins),
    period_rhr = mean(mean_hr),
    min_rhr = min(mean_hr),
    max_rhr = max(mean_hr),
    sd_rhr = sd(mean_hr),
    earliest_date = min(start_time),
    latest_date = max(start_time),
    .groups = "drop"
  )

# 计算术前摘要
rhr_summary_pre_surgery_50 <- rhr_by_window_50_surgery %>%
  filter(time_period %in% c("pre_surgery_over_7d", "pre_surgery_3d_to_7d", "pre_surgery_3d")) %>%
  group_by(subject_id) %>%
  summarise(
    time_period = "all_pre_surgery",
    n_windows = n(),
    total_rest_minutes = sum(duration_mins),
    period_rhr = mean(mean_hr),
    min_rhr = min(mean_hr),
    max_rhr = max(mean_hr),
    sd_rhr = sd(mean_hr),
    .groups = "drop"
  )

# 合并所有摘要
rhr_summary_complete_50 <- bind_rows(
  rhr_summary_by_period_50,
  rhr_summary_pre_surgery_50
)

# 总体RHR统计
rhr_summary_50 <- rhr_by_window_50 %>%
  group_by(subject_id) %>%
  summarise(
    n_windows = n(),
    total_rest_minutes = sum(duration_mins),
    overall_rhr = mean(mean_hr),
    min_rhr = min(mean_hr),
    max_rhr = max(mean_hr),
    sd_rhr = sd(mean_hr)
  )

# 添加总体RHR
rhr_summary_complete_50 <- rhr_summary_complete_50 %>%
  left_join(
    rhr_summary_50 %>% 
      dplyr::select(subject_id, overall_rhr),
    by = "subject_id"
  )

# 创建更宽的格式用于比较
rhr_summary_wide_50 <- rhr_summary_complete_50 %>%
  dplyr::select(subject_id, time_period, period_rhr) %>%
  pivot_wider(
    names_from = time_period,
    values_from = period_rhr,
    names_prefix = "rhr_50_"
  )

# 保存结果
save(rhr_by_window_50_surgery, 
     rhr_summary_complete_50,
     rhr_summary_wide_50,
     file = "rhr_analysis_50.rda",
     compress = "xz")
