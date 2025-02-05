library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)
library(lubridate)
library(zoo)


###read data
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")
baseline_info <- read.csv("2_data/analysis_data/sample_info.csv")


######
dir.create("3_data_analysis/2_data_analysis/prediction/1_heart_rate", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/prediction/1_heart_rate")




# 处理视力数据
process_va_data <- function(baseline_info) {
  va_1m <- baseline_info %>%
    mutate(
      va_1m = case_when(
        surger_eye == 0 ~ od_corrected_3m,  # 右眼手术
        surger_eye == 1 ~ os_corrected_3m,  # 左眼手术
        surger_eye == 2 ~ pmax(od_corrected_3m, os_corrected_3m, na.rm = TRUE),  
        TRUE ~ NA_real_  # 其他情况
      )
    ) %>%
    dplyr::select(ID, va_1m, surger_eye)
  
  return(va_3m)
}

va_data <- process_va_data(baseline_info)

# 从massdataset中提取数据
heart_rate_matrix <- extract_expression_data(heart_rate_data)  # 获取心率值矩阵
sample_info <- heart_rate_data@sample_info  # 获取sample_info

# 转换数据为数据框格式
heart_rate_df <- data.frame(
  measure_time = sample_info$measure_time,
  subject_id = sample_info$subject_id,
  heart_rate = as.numeric(heart_rate_matrix[1,])
)



#####
analyze_heart_rate <- function(heart_rate_df, va_data) {
  
  # Convert subject_id to character type if not already
  heart_rate_df$subject_id <- as.character(heart_rate_df$subject_id)
  
  # 计算7天窗口的基线特征
  hr_features <- heart_rate_df %>%
    group_by(subject_id) %>%
    summarise(
      # 基础心率统计
      mean_hr = mean(heart_rate, na.rm = TRUE),
      sd_hr = sd(heart_rate, na.rm = TRUE),
      min_hr = min(heart_rate, na.rm = TRUE),
      max_hr = max(heart_rate, na.rm = TRUE),
      
      # 计算每天的平均心率，然后计算7天的变异系数
      daily_cv = {
        # Calculate daily means within the group
        daily_means <- summarise(
          group_by(
            data.frame(
              heart_rate = heart_rate,
              date = as.Date(measure_time)
            ),
            date
          ),
          daily_mean = mean(heart_rate, na.rm = TRUE)
        )$daily_mean
        
        # Calculate CV if we have data
        if(length(daily_means) > 0) {
          (sd(daily_means, na.rm = TRUE) / mean(daily_means, na.rm = TRUE) * 100)
        } else {
          NA_real_
        }
      },
      
      # 计算心率范围
      hr_range = max_hr - min_hr,
      
      # 计算异常心率的比例 (超过2个标准差)
      abnormal_ratio = {
        z_scores <- scale(heart_rate)
        sum(abs(z_scores) > 2, na.rm = TRUE) / length(z_scores)
      }
    )
  
  # 合并心率特征和视力数据
  final_data <- hr_features %>%
    left_join(va_data %>% 
                mutate(subject_id = as.character(ID)) %>%  # Ensure ID is character type
                dplyr::select(-ID), 
              by = "subject_id")
  
  # 构建预测模型
  model <- lm(va_1m ~ mean_hr + sd_hr + daily_cv + hr_range + abnormal_ratio,
              data = final_data)
  
  # 返回结果列表
  return(list(
    features = final_data,
    model = model,
    preprocessed_data = heart_rate_df
  ))
}

# Add some diagnostic prints
print("Number of unique subject_ids:")
print(length(unique(heart_rate_df$subject_id)))

print("\nFirst few rows of heart_rate_df:")
print(head(heart_rate_df))

# Run the analysis
results <- analyze_heart_rate(heart_rate_df, va_data)

# 检查结果
summary(results$model)

# 查看特征数据
head(results$features)

# 评估模型性能
predicted_va <- predict(results$model, results$features)
rmse <- sqrt(mean((results$features$va_1m - predicted_va)^2, na.rm = TRUE))
r2 <- summary(results$model)$r.squared

cat("\n模型性能：\n")
cat("RMSE:", rmse, "\n")
cat("R-squared:", r2, "\n")

# 1. 可视化预测结果与实际值的对比
library(ggplot2)
library(tidyverse)

# 创建预测值与实际值的比较图
ggplot(data = results$features, aes(x = va_1m, y = predicted_va)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "预测视力值 vs 实际视力值",
       x = "实际视力值",
       y = "预测视力值") +
  theme_minimal() +
  coord_equal()

# 2. 特征重要性可视化
coef_data <- data.frame(
  feature = c("mean_hr", "sd_hr", "daily_cv", "hr_range", "abnormal_ratio"),
  coefficient = coef(results$model)[-1],  # 排除截距
  p_value = summary(results$model)$coefficients[-1, "Pr(>|t|)"]
) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ ".",
    TRUE ~ "ns"
  ))

# 绘制特征重要性条形图
ggplot(coef_data, aes(x = reorder(feature, abs(coefficient)), y = coefficient)) +
  geom_col(aes(fill = p_value < 0.05)) +
  geom_text(aes(label = significance), vjust = -0.5) +
  labs(title = "特征重要性分析",
       x = "特征",
       y = "回归系数",
       fill = "统计显著") +
  theme_minimal() +
  coord_flip()

# 3. 残差分析图
residuals <- resid(results$model)
fitted_values <- fitted(results$model)

par(mfrow = c(2,2))
plot(results$model)
par(mfrow = c(1,1))

# 4. 相关性热图
library(corrplot)
cor_matrix <- cor(results$features %>% 
                    dplyr::select(mean_hr, sd_hr, daily_cv, hr_range, abnormal_ratio, va_1m),
                  use = "complete.obs")
corrplot(cor_matrix, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45)

# 5. 各特征与视力的散点图
features_long <- results$features %>%
  dplyr::select(mean_hr, sd_hr, daily_cv, hr_range, abnormal_ratio, va_1m) %>%
  gather(key = "feature", value = "value", -va_1m)

ggplot(features_long, aes(x = value, y = va_1m)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~feature, scales = "free_x") +
  labs(title = "各特征与视力的关系",
       x = "特征值",
       y = "视力值") +
  theme_minimal()

