# 术后OCTA指标改善预测模型
# 使用术前7天可穿戴设备数据预测术后1个月OCTA指标改善值(T2-T0)

library(tidyverse)
library(caret)
library(gridExtra)
library(pROC)
library(glmnet)
library(ggplot2)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

#-------------------------------
# 1. 加载数据
#-------------------------------

# 加载基础数据
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 创建输出目录
dir.create("3_data_analysis/5_presurgery_analysis/octa_prediction/results", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/octa_prediction/results")

#-------------------------------
# 2. 处理OCTA数据 - 区分血流和厚度
#-------------------------------

# 处理OCTA血流数据
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

# 处理每个患者数据的函数
process_patient_data <- function(patient_data, time_points = c("T0", "T2")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 处理每个时间点
  for(suffix in time_points) {
    # 选择当前时间点的列
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # 选择数据并重命名列
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# 处理所有患者的数据
process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# 分别处理血流和厚度数据
octa_bloodflow_features <- process_all_patients(octa_bloodflow_features)
octa_thickness_features <- process_all_patients(octa_thickness_features)

#-------------------------------
# 3. 计算OCTA指标的改善值(T2-T0)
#-------------------------------

# 计算改善值函数
calculate_improvement <- function(data, data_type = "unknown") {
  # 初始化结果数据框
  result <- data %>% dplyr::select(ID)
  
  # 对每个T0变量找到对应的T2变量并计算差值
  vars_T0 <- names(data)[grep("_T0$", names(data))]
  improvement_count <- 0
  
  for(t0_var in vars_T0) {
    t2_var <- gsub("_T0$", "_T2", t0_var)
    
    if(t2_var %in% names(data)) {
      # 构建改善值变量名
      imp_var <- gsub("_T0$", "_improvement", t0_var)
      
      # 计算改善值(T2-T0)
      result[[imp_var]] <- data[[t2_var]] - data[[t0_var]]
      improvement_count <- improvement_count + 1
    }
  }
  
  cat("计算了", improvement_count, "个", data_type, "改善值变量\n")
  return(result)
}

# 计算血流和厚度改善值
bloodflow_improvement <- calculate_improvement(octa_bloodflow_features, "血流")
thickness_improvement <- calculate_improvement(octa_thickness_features, "厚度")

#-------------------------------
# 4. 准备可穿戴设备预测变量
#-------------------------------

# 筛选术前数据(-7天至-1天)
pre_surgery_data <- combined_data %>%
  filter(as.numeric(day) >= -7 & as.numeric(day) <= -1)

# 按受试者ID分组并计算术前均值
patient_preop_means <- pre_surgery_data %>%
  group_by(subject_id) %>%
  summarise(
    # 选择关键均值指标
    mean_rhr_1 = mean(mean_rhr_1, na.rm = TRUE),     # 平均静息心率
    mean_bo = mean(mean_bo, na.rm = TRUE),           # 平均血氧饱和度
    steps_total = mean(steps_total, na.rm = TRUE),   # 平均总步数
    total_sleep = mean(total_sleep, na.rm = TRUE),   # 平均总睡眠时间
    
    # 保留其他人口统计变量
    dm_2 = dplyr::first(dm_2),                    # 糖尿病状态
    age = dplyr::first(age),                      # 年龄
    gender = dplyr::first(gender),                # 性别
    bmi = dplyr::first(bmi),                      # BMI
    hypertension_2 = dplyr::first(hypertension_2) # 高血压状态
  )

# 中位数填补缺失值
data_imputed <- patient_preop_means
for(col in names(data_imputed)[2:(ncol(data_imputed)-5)]) {  # 跳过subject_id和人口统计变量
  if(sum(is.na(data_imputed[[col]])) > 0) {
    data_imputed[[col]][is.na(data_imputed[[col]])] <- median(data_imputed[[col]], na.rm = TRUE)
  }
}

#-------------------------------
# 5. 合并可穿戴数据与OCTA改善值数据
#-------------------------------

# 合并厚度数据
model_data_thickness <- data_imputed %>%
  left_join(thickness_improvement, by = c("subject_id" = "ID"))

# 合并血流数据
model_data_bloodflow <- data_imputed %>%
  left_join(bloodflow_improvement, by = c("subject_id" = "ID"))

#-------------------------------
# 6. 定义要分析的变量
#-------------------------------

# 血流变量 - PA_前缀为血流灌注区域(Perfusion Area)
significant_bloodflow_pairs <- list(
  list(improvement_var = "PA_Vitreous_0_6_improvement")
)

# 厚度变量
significant_thickness_pairs <- list(
  list(improvement_var = "Thickness_OuterRetina_0_6_improvement"),
  list(improvement_var = "Thickness_Retina_0_6_improvement"),
  list(improvement_var = "Thickness_RNFL_0_6_improvement")
)

# 定义可穿戴设备的预测变量
wearable_vars <- list(
  "mean_rhr_1",    # 静息心率
  "mean_bo",       # 血氧
  "steps_total"    # 步数
)

#-------------------------------
# 7. 构建预测模型函数
#-------------------------------

# 预测模型构建函数 - 适用于血流和厚度改善
build_prediction_models <- function(model_data, target_pairs, model_type = "unknown") {
  # 查看模型数据的变量名
  cat("模型数据变量名：", paste(names(model_data)[1:10], collapse = ", "), "...\n")
  
  # 创建数据框来存储每个目标变量和预测器组合的模型性能指标
  performance_df <- data.frame(
    target_var = character(),
    predictor = character(),
    r_squared = numeric(),
    adj_r_squared = numeric(),
    cv_r_squared = numeric(),
    cv_rmse = numeric(),
    sample_size = numeric(),
    predictor_coef = numeric(),
    predictor_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 对每个目标变量进行建模
  for(pair in target_pairs) {
    improvement_var <- pair$improvement_var
    
    # 检查改善值变量是否存在
    if(!(improvement_var %in% names(model_data))) {
      cat(model_type, "改善值变量不存在:", improvement_var, "\n")
      next
    }
    
    cat("\n===== 构建", model_type, "改善预测模型:", improvement_var, "=====\n")
    
    # 准备基本模型数据 - 移除含有NA的行
    base_vars <- c(improvement_var, "age", "gender", "bmi", "dm_2")
    base_data <- model_data %>%
      select(all_of(base_vars)) %>%
      na.omit()
    
    if(nrow(base_data) < 10) {  # 确保有足够的样本
      cat("样本量不足:", improvement_var, "\n")
      next
    }
    
    cat("样本量:", nrow(base_data), "\n")
    
    # 检查base_data中是否存在所需变量
    cat("base_data中的变量:", paste(names(base_data), collapse=", "), "\n")
    
    # 检查并调整wearable_vars列表，确保使用实际存在的变量名
    wearable_exists <- sapply(wearable_vars, function(var) var %in% names(model_data))
    if(!all(wearable_exists)) {
      cat("警告: 以下可穿戴设备变量在数据中不存在:", 
          paste(wearable_vars[!wearable_exists], collapse=", "), "\n")
      # 尝试修正 - 从模型数据中获取变量名
      actual_vars <- names(model_data)[2:4]  # 假设前几个变量是可穿戴设备变量
      cat("使用实际存在的变量:", paste(actual_vars, collapse=", "), "\n")
      wearable_vars <- as.list(actual_vars)
    }
    
    # 准备基本模型数据 - 现在包含可穿戴设备变量
    base_vars <- c(improvement_var, unlist(wearable_vars), "age", "gender", "bmi", "dm_2")
    base_data <- model_data %>%
      select(all_of(base_vars)) %>%
      na.omit()
    
    # 定义交叉验证函数
    cv_model_performance <- function(model_formula, model_data) {
      # 创建训练控制对象
      train_control <- trainControl(
        method = "cv",
        number = 5,
        savePredictions = "final"
      )
      
      # 使用caret的train函数进行交叉验证
      set.seed(123)
      cv_model <- tryCatch({
        train(
          model_formula,
          data = model_data,
          method = "lm",
          trControl = train_control
        )
      }, error = function(e) {
        cat("交叉验证错误:", e$message, "\n")
        return(NULL)
      })
      
      if(is.null(cv_model)) {
        # 如果交叉验证失败，使用常规线性模型
        lm_model <- lm(model_formula, data = model_data)
        return(list(
          r_squared = summary(lm_model)$r.squared,
          rmse = sqrt(mean(lm_model$residuals^2))
        ))
      }
      
      # 返回交叉验证性能指标
      return(list(
        r_squared = cv_model$results$Rsquared,
        rmse = cv_model$results$RMSE
      ))
    }
    
    # 合并前三个变量为公式字符串
    wearable_formula_str <- paste(unlist(wearable_vars), collapse = " + ")
    
    # 1. 所有预测变量模型
    all_formula <- as.formula(paste(
      improvement_var, "~", 
      wearable_formula_str, " + ",
      "age + gender + bmi + dm_2"
    ))
    
    all_model <- lm(all_formula, data = base_data)
    all_summary <- summary(all_model)
    all_cv <- cv_model_performance(all_formula, base_data)
    
    # 添加到性能数据框
    performance_df <- rbind(performance_df, data.frame(
      target_var = improvement_var,
      predictor = "All Predictors",
      r_squared = all_summary$r.squared,
      adj_r_squared = all_summary$adj.r.squared,
      cv_r_squared = all_cv$r_squared,
      cv_rmse = all_cv$rmse,
      sample_size = nrow(base_data),
      predictor_coef = NA,  # 对于全部预测变量模型，不记录单个系数
      predictor_p_value = NA,
      stringsAsFactors = FALSE
    ))
    
    # 2. 对每个可穿戴设备变量构建单独的模型
    for(wearable_var in wearable_vars) {
      # 创建模型公式
      var_formula <- as.formula(paste(
        improvement_var, "~", 
        wearable_var, "+", 
        "age + gender + bmi + dm_2"
      ))
      
      # 构建模型
      var_model <- lm(var_formula, data = base_data)
      var_summary <- summary(var_model)
      var_cv <- cv_model_performance(var_formula, base_data)
      
      # 获取预测变量的系数和p值
      coef_data <- coef(var_summary)
      var_coef <- coef_data[wearable_var, 1]  # 系数
      var_p_value <- coef_data[wearable_var, 4]  # p值
      
      # 添加到性能数据框
      performance_df <- rbind(performance_df, data.frame(
        target_var = improvement_var,
        predictor = wearable_var,
        r_squared = var_summary$r.squared,
        adj_r_squared = var_summary$adj.r.squared,
        cv_r_squared = var_cv$r_squared,
        cv_rmse = var_cv$rmse,
        sample_size = nrow(base_data),
        predictor_coef = var_coef,
        predictor_p_value = var_p_value,
        stringsAsFactors = FALSE
      ))
      
      # 打印单个预测变量模型信息
      cat("\n--- 单一预测变量:", wearable_var, "---\n")
      cat("R²:", var_summary$r.squared, "\n")
      cat("调整后R²:", var_summary$adj.r.squared, "\n")
      cat("交叉验证R²:", var_cv$r_squared, "\n")
      cat("交叉验证RMSE:", var_cv$rmse, "\n")
      cat(wearable_var, "系数:", var_coef, "p值:", var_p_value, "\n")
    }
    
    # 创建预测器性能对比图
    plot_data <- performance_df %>%
      filter(target_var == improvement_var) %>%
      select(predictor, cv_r_squared)
    
    predictor_names <- c("All Predictors", unlist(wearable_vars))
    plot_data$predictor <- factor(plot_data$predictor, levels = predictor_names)
    
    # 安全地定义颜色映射
    color_values <- c("All Predictors" = "#513a52")
    for(i in seq_along(wearable_vars)) {
      var_name <- unlist(wearable_vars)[i]
      color_values[var_name] <- c("#82a1bf", "#faaa93", "#feefc4")[i]
    }
    
    p <- ggplot(plot_data, aes(x = predictor, y = cv_r_squared, fill = predictor)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = color_values) +
      labs(
        title = paste(model_type, "改善预测:", improvement_var),
        x = "预测变量",
        y = "交叉验证 R²",
        fill = "预测变量"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    
    # 保存图表
    chart_filename <- paste0(model_type, "_prediction_performance_", gsub("_improvement", "", improvement_var), ".pdf")
    ggsave(chart_filename, p, width = 8, height = 6)
    
    # 为最佳模型创建预测散点图
    best_row <- which.max(plot_data$cv_r_squared)
    best_predictor <- as.character(plot_data$predictor[best_row])
    
    if(best_predictor == "All Predictors") {
      best_model <- all_model
    } else {
      best_formula <- as.formula(paste(
        improvement_var, "~", 
        best_predictor, "+", 
        "age + gender + bmi + dm_2"
      ))
      best_model <- lm(best_formula, data = base_data)
    }
    
    # 创建预测散点图
    fitted_values <- predict(best_model)
    scatter_data <- data.frame(
      Actual = base_data[[improvement_var]],
      Predicted = fitted_values
    )
    
    p_scatter <- ggplot(scatter_data, aes(x = Actual, y = Predicted)) +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = paste(model_type, "改善预测:", improvement_var),
        subtitle = paste("最佳预测变量:", best_predictor),
        x = "实际改善值",
        y = "预测改善值"
      ) +
      theme_bw() +
      annotate("text", x = min(scatter_data$Actual), y = max(scatter_data$Predicted),
               label = paste("R² =", sprintf("%.3f", summary(best_model)$r.squared),
                             "\nRMSE =", sprintf("%.3f", sqrt(mean(best_model$residuals^2)))),
               hjust = 0, vjust = 1)
    
    # 保存散点图
    scatter_filename <- paste0(model_type, "_best_prediction_scatter_", gsub("_improvement", "", improvement_var), ".pdf")
    ggsave(scatter_filename, p_scatter, width = 8, height = 6)
  }
  
  return(performance_df)
}

#-------------------------------
# 8. 运行预测模型
#-------------------------------

cat("\n开始构建血流改善预测模型...\n")
bloodflow_performance <- build_prediction_models(model_data_bloodflow, significant_bloodflow_pairs, "bloodflow")

cat("\n开始构建厚度改善预测模型...\n")
thickness_performance <- build_prediction_models(model_data_thickness, significant_thickness_pairs, "thickness")

#-------------------------------
# 9. 保存模型性能结果
#-------------------------------

# 合并所有性能结果
all_performance <- rbind(
  bloodflow_performance,
  thickness_performance
)

# 保存到CSV文件
write.csv(all_performance, "octa_improvement_prediction_performance.csv", row.names = FALSE)

# 创建血流和厚度模型性能比较图
for(model_type in c("bloodflow", "thickness")) {
  # 筛选当前模型类型的数据
  current_data <- all_performance %>%
    filter(grepl(paste0("^", ifelse(model_type == "bloodflow", "PA|VD|SVD", "Thickness")), target_var))
  
  if(nrow(current_data) > 0) {
    # 准备绘图数据
    plot_data <- data.frame()
    
    # 获取所有目标变量
    unique_targets <- unique(current_data$target_var)
    
    # 对每个目标变量提取预测器性能
    for(target in unique_targets) {
      target_data <- current_data %>% 
        filter(target_var == target)
      
      # 确保所有预测器都存在
      for(pred in c("All Predictors", unlist(wearable_vars))) {
        if(!pred %in% target_data$predictor) {
          # 如果预测器不存在，添加NA行
          target_data <- rbind(target_data, data.frame(
            target_var = target,
            predictor = pred,
            cv_r_squared = NA,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      plot_data <- rbind(plot_data, target_data)
    }
    
    # 简化变量名用于图表显示
    plot_data$short_name <- gsub("_improvement$", "", plot_data$target_var)
    plot_data$short_name <- gsub("^(PA|VD|SVD|Thickness)_", "", plot_data$short_name)
    
    # 创建图表
    predictor_names <- c("All Predictors", unlist(wearable_vars))
    plot_data$predictor <- factor(plot_data$predictor, levels = predictor_names)
    
    # 安全地定义颜色映射
    color_values <- c("All Predictors" = "#513a52")
    for(i in seq_along(wearable_vars)) {
      var_name <- unlist(wearable_vars)[i]
      color_values[var_name] <- c("#82a1bf", "#faaa93", "#feefc4")[i]
    }
    
    p <- ggplot(plot_data, aes(x = short_name, y = cv_r_squared, fill = predictor)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      scale_fill_manual(values = color_values) +
      labs(
        title = paste(ifelse(model_type == "bloodflow", "血流", "厚度"), "改善预测模型性能比较"),
        x = "OCTA指标",
        y = "交叉验证 R²",
        fill = "预测变量"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    
    # 保存图表
    comparison_filename <- paste0(model_type, "_models_performance_comparison.pdf")
    ggsave(comparison_filename, p, width = 12, height = 8)
  }
}

cat("\n所有模型性能结果已保存到 'octa_improvement_prediction_performance.csv'\n")
cat("过程完成！\n")