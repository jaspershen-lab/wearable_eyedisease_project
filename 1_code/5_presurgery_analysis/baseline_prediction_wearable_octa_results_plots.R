library(tidyverse)
library(ggpubr)     # 用于统计测试
library(pheatmap)   # 用于热图
library(RColorBrewer)
library(forestplot)
library(scales)     # 用于更好的比例尺
library(r4projects)
setwd(get_project_wd())  
rm(list = ls())
source("1_code/100_tools.R")

# 第1步: 加载第一步分析生成的数据
comparison_data <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/comparison_data.csv")
thickness_correlations <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/thickness_correlations.csv")
bloodflow_correlations <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/bloodflow_correlations.csv")
thickness_regression <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/thickness_regression_results.csv")
bloodflow_regression <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/bloodflow_regression_results.csv")
model_performance <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/model_performance.csv")
long_format_data <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/long_format_data.csv")
top_model_predictions <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/top_model_predictions.csv")
wearable_thickness_data <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/wearable_thickness_data.csv")
wearable_bloodflow_data <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/wearable_bloodflow_data.csv")

# 显著相关性
significant_thickness <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/significant_thickness_correlations.csv")
significant_bloodflow <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/significant_bloodflow_correlations.csv")
top_thickness_vars <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/top_thickness_vars.csv")
top_bloodflow_vars <- read.csv("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/top_bloodflow_vars.csv")

# 创建输出目录
dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/plots", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv/plots")

# 2. 热图函数
create_heatmap <- function(data, title, filename, significant_only = FALSE, 
                           width = 8, height = 7, units = "in") {
  if(significant_only) data <- data %>% filter(P_value < 0.05)
  
  plot_data <- data %>% 
    filter(!is.na(Correlation)) %>%
    mutate(
      significance = case_when(
        P_value < 0.001 ~ "***", P_value < 0.01 ~ "**", 
        P_value < 0.05 ~ "*", TRUE ~ ""
      ),
      text_label = sprintf("%.2f%s", Correlation, significance),
      Wearable_Label = case_when(
        Wearable_Variable == "mean_hr" ~ "Mean HR",
        Wearable_Variable == "mean_rhr" ~ "Resting HR",
        Wearable_Variable == "mean_bo" ~ "Blood Oxygen",
        Wearable_Variable == "total_steps" ~ "Steps",
        Wearable_Variable == "total_sleep" ~ "Sleep",
        TRUE ~ Wearable_Variable
      )
    )
  
  p <- ggplot(plot_data, 
              aes(x = paste(Wearable_Label, Group, sep = "\n"), 
                  y = OCTA_Variable, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = text_label), size = 2.5) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", 
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = title) +
    theme_custom +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  # 保存为PDF而非PNG
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height, 
         units = units, device = pdf)
  
  # 也保存PNG版本
  ggsave(paste0(filename, ".png"), p, width = width, height = height, 
         units = units, dpi = 300)
}

# 创建热图，使用自定义尺寸
create_heatmap(thickness_correlations, "Thickness Correlations", "thickness_heatmap",
               width = 10, height = 8)
create_heatmap(bloodflow_correlations, "Blood Flow Correlations", "bloodflow_heatmap",
               width = 10, height = 8)
create_heatmap(thickness_correlations, "Significant Thickness Correlations", 
               "thickness_sig_heatmap", significant_only = TRUE,
               width = 8, height = 6)
create_heatmap(bloodflow_correlations, "Significant Blood Flow Correlations", 
               "bloodflow_sig_heatmap", significant_only = TRUE,
               width = 8, height = 6)

# 3. 森林图函数
# 简单森林图函数
simple_forest_plot <- function(data, title) {
  # 筛选有显著相关性的数据
  sig_data <- data %>%
    filter(P_value < 0.05) %>%
    mutate(
      # 添加置信区间 (如果不存在)
      z = if("z" %in% names(.)) z else 0.5 * log((1 + Correlation) / (1 - Correlation)),
      se = if("se" %in% names(.)) se else 1 / sqrt(N - 3),
      CI_Lower = if("CI_Lower" %in% names(.)) CI_Lower else tanh(z - 1.96 * se),
      CI_Upper = if("CI_Upper" %in% names(.)) CI_Upper else tanh(z + 1.96 * se),
      
      # 美化标签
      Wearable_Label = case_when(
        Wearable_Variable == "mean_hr" ~ "Mean HR",
        Wearable_Variable == "mean_rhr" ~ "Resting HR",
        Wearable_Variable == "mean_bo" ~ "Blood O₂",
        Wearable_Variable == "total_steps" ~ "Steps",
        Wearable_Variable == "total_sleep" ~ "Sleep",
        TRUE ~ Wearable_Variable
      ),
      # 创建标签
      Pair = paste(OCTA_Variable, "-", Wearable_Label),
      Label = paste(Group, Pair)
    ) %>%
    arrange(Pair, desc(Group))
  
  # 绘制图表
  p <- ggplot(sig_data, 
              aes(y = Label, x = Correlation, 
                  xmin = CI_Lower, xmax = CI_Upper, 
                  color = Group)) +
    geom_point(size = 3) +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = diabetes_colors) +
    labs(title = title, x = "Correlation Coefficient", y = "") +
    theme_custom
  
  # 返回图表对象
  return(p)
}

# 创建并显示厚度相关性森林图
p1 <- simple_forest_plot(thickness_correlations, "Thickness Correlations")
print(p1)
p2 <- simple_forest_plot(bloodflow_correlations, "Blood Flow Correlations")
print(p2)

ggsave("Thickness_Correlations_forset_plot.pdf",plot=p1,width = 7, height = 8)
ggsave("BloodFlow_Correlations_forset_plot.pdf",plot=p2,width = 7, height = 8)

# 4. 散点图函数
create_scatterplots <- function(correlations, var_type, n_plots = 5) {
  # 加载原始数据
  data_file <- paste0("wearable_", var_type, "_data.csv")
  if(file.exists(data_file)) {
    original_data <- read.csv(data_file)
    
    # 前n个显著相关性
    top_correlations <- correlations %>%
      filter(P_value < 0.05) %>%
      arrange(P_value) %>%
      head(n_plots)
    
    for(i in 1:nrow(top_correlations)) {
      pair <- top_correlations[i, ]
      octa_var <- pair$OCTA_Variable
      wearable_var <- pair$Wearable_Variable
      group <- pair$Group
      
      # 准备数据
      plot_data <- original_data %>% 
        filter(ifelse(group == "Diabetes", dm_2 == 1, dm_2 == 0)) %>%
        dplyr::select(subject_id, all_of(c(octa_var, wearable_var))) %>%
        na.omit()
      
      wearable_label <- case_when(
        wearable_var == "mean_hr" ~ "Mean HR",
        wearable_var == "mean_rhr" ~ "Resting HR",
        wearable_var == "mean_bo" ~ "Blood Oxygen",
        wearable_var == "total_steps" ~ "Steps",
        wearable_var == "total_sleep" ~ "Sleep",
        TRUE ~ wearable_var
      )
      
      p <- ggplot(plot_data, aes_string(x = wearable_var, y = octa_var)) +
        geom_point(alpha = 0.7, color = ifelse(group == "Diabetes", "#D6604D", "#4393C3")) +
        geom_smooth(method = "lm", color = ifelse(group == "Diabetes", "#D6604D", "#4393C3")) +
        labs(title = paste(octa_var, "vs", wearable_label),
             subtitle = sprintf("%s: r = %.2f, p = %.3f", group, pair$Correlation, pair$P_value),
             x = wearable_label, y = octa_var) +
        theme_custom
      
      print(p)
      filename <- paste0("scatter_", var_type, "_", tolower(gsub(" ", "_", group)), "_", i)
      ggsave(paste0(filename, ".png"), p, width = 7, height = 5, dpi = 300)
    }
  }
}

# 创建散点图
create_scatterplots(significant_thickness, "thickness")
create_scatterplots(significant_bloodflow, "bloodflow")

# 5. 回归系数森林图函数
create_reg_forest_with_significance <- function(data, title, filename) {
  # 只选择显著的回归系数
  sig_data <- data %>%
    filter(P_value < 0.05,
           Variable %in% c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep"))
  
  if(nrow(sig_data) > 0) {
    plot_data <- sig_data %>%
      mutate(
        # 添加显著性标记
        significance = case_when(
          P_value < 0.001 ~ "***",
          P_value < 0.01 ~ "**",
          P_value < 0.05 ~ "*",
          TRUE ~ ""
        ),
        # 美化变量标签
        Wearable_Label = case_when(
          Variable == "mean_hr" ~ "Mean HR",
          Variable == "mean_rhr" ~ "Resting HR",
          Variable == "mean_bo" ~ "Blood Oxygen",
          Variable == "total_steps" ~ "Steps",
          Variable == "total_sleep" ~ "Sleep",
          TRUE ~ Variable
        ),
        Pair_Label = paste(OCTA_Variable, "-", Wearable_Label),
        Group_Pair = paste(Group, Pair_Label),
        # 为标签添加显著性标记
        Group_Pair_Sig = paste0(Group_Pair, " ", significance)
      ) %>%
      arrange(Pair_Label, desc(Group))
    
    # 创建森林图
    p <- ggplot(plot_data, 
                aes(y = Group_Pair, x = Coefficient, 
                    xmin = Lower_CI, xmax = Upper_CI, 
                    color = Group)) +
      geom_point(size = 3) +
      geom_errorbarh(height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = diabetes_colors) +
      # 添加显著性标记作为额外的文本
      geom_text(aes(label = significance, x = Upper_CI + 
                      0.05 * (max(Upper_CI) - min(Lower_CI))), 
                hjust = 0, size = 4, fontface = "bold") +
      labs(title = title, 
           x = "Regression Coefficient", 
           caption = "sig: * p<0.05, ** p<0.01, *** p<0.001") +
      theme_custom +
      # 确保有足够的空间显示显著性标记
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))
    
    print(p)
    ggsave(paste0(filename, ".pdf"), p, width = 10, 
           height = max(8, length(unique(plot_data$Pair_Label)) * 0.6),
           dpi = 300)
    ggsave(paste0(filename, ".png"), p, width = 10, 
           height = max(8, length(unique(plot_data$Pair_Label)) * 0.6),
           dpi = 300)
    
  }
}

# 创建带显著性标记的回归森林图
create_reg_forest_with_significance(thickness_regression, "Thickness Regression Coefficients", "thickness_reg_forest_sig")
create_reg_forest_with_significance(bloodflow_regression, "Blood Flow Regression Coefficients", "bloodflow_reg_forest_sig")

# 创建回归森林图
create_reg_forest(thickness_regression, "Thickness Regression Coefficients", "thickness_reg_forest")
create_reg_forest(bloodflow_regression, "Blood Flow Regression Coefficients", "bloodflow_reg_forest")

# 6. 模型性能图函数
create_performance_plot <- function(data, metric = "R_squared", title, filename) {
  p <- ggplot(data, aes(x = reorder(OCTA_Variable, !!sym(metric)), 
                        y = !!sym(metric), 
                        fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = diabetes_colors) +
    labs(title = title, 
         y = ifelse(metric == "R_squared", "R²", "RMSE"), 
         x = "OCTA Variables") +
    theme_custom +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  ggsave(paste0(filename, ".png"), p, width = 9, height = 5, dpi = 300)
}

# 分离数据
thickness_model_perf <- model_performance %>% filter(grepl("Thickness", OCTA_Variable))
bloodflow_model_perf <- model_performance %>% filter(!grepl("Thickness", OCTA_Variable))

# 创建模型性能图
create_performance_plot(thickness_model_perf, "R_squared", 
                        "Thickness Model R²", "thickness_r2")
create_performance_plot(thickness_model_perf, "RMSE", 
                        "Thickness Model RMSE", "thickness_rmse")
create_performance_plot(bloodflow_model_perf, "R_squared", 
                        "Blood Flow Model R²", "bloodflow_r2")
create_performance_plot(bloodflow_model_perf, "RMSE", 
                        "Blood Flow Model RMSE", "bloodflow_rmse")




# 7. 预测值-实际值对比图（使用留一法交叉验证）
library(caret)

# 准备LOOCV预测结果的函数
perform_loocv <- function(data, octa_var, wearable_vars, group_filter) {
  # 确保所有变量名存在于数据中
  all_vars <- c(octa_var, wearable_vars)
  missing_vars <- all_vars[!all_vars %in% colnames(data)]
  if(length(missing_vars) > 0) {
    warning(paste("Missing variables in dataset:", paste(missing_vars, collapse=", ")))
    return(NULL)
  }
  
  # 过滤特定组的数据
  group_data <- data %>% 
    filter(!!group_filter) %>%
    select(subject_id, all_of(all_vars)) %>%
    drop_na()  # 删除任何包含NA的行
  
  # 确保有足够的数据点
  if(nrow(group_data) < 3) {
    warning(paste("Too few data points for", octa_var, ". Only", nrow(group_data), "complete cases."))
    return(NULL)  # 返回NULL如果数据点不足
  }
  
  # 初始化结果向量
  n <- nrow(group_data)
  actual_values <- group_data[[octa_var]]
  predicted_values <- numeric(n)
  
  # 对每个样本执行留一法交叉验证
  for(i in 1:n) {
    # 分离测试样本和训练样本
    test_sample <- group_data[i, ]
    train_samples <- group_data[-i, ]
    
    # 构建公式字符串，确保所有预测变量存在
    valid_predictors <- intersect(wearable_vars, colnames(train_samples))
    
    # 如果没有有效的预测变量，则跳过
    if(length(valid_predictors) == 0) {
      warning(paste("No valid predictors for", octa_var))
      return(NULL)
    }
    
    formula_str <- paste(octa_var, "~", paste(valid_predictors, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
    # 训练模型
    tryCatch({
      model <- lm(formula_obj, data = train_samples)
      # 预测
      predicted_values[i] <- predict(model, newdata = test_sample)
    }, error = function(e) {
      warning(paste("Error in model fitting for sample", i, ":", e$message))
      predicted_values[i] <<- NA
    })
  }
  
  # 删除任何NA预测
  valid_indices <- !is.na(predicted_values)
  if(sum(valid_indices) < 3) {
    warning(paste("Too few valid predictions for", octa_var))
    return(NULL)
  }
  
  actual_values <- actual_values[valid_indices]
  predicted_values <- predicted_values[valid_indices]
  subjects <- group_data$subject_id[valid_indices]
  
  # 计算性能指标
  rmse_val <- sqrt(mean((predicted_values - actual_values)^2))
  r_squared <- cor(predicted_values, actual_values)^2
  
  # 返回结果
  return(data.frame(
    actual = actual_values,
    predicted = predicted_values,
    subject_id = subjects,
    r_squared = r_squared,
    rmse = rmse_val
  ))
}

# 运行LOOCV并创建预测图
create_loocv_prediction_plots <- function(data_df, octa_vars, wearable_vars, var_type) {
  # 使用传入的数据框
  original_data <- data_df
  
  # 检查变量是否存在
  print(paste("检查OCTA变量:", paste(head(octa_vars), collapse=", ")))
  print(paste("检查可穿戴设备变量:", paste(head(wearable_vars), collapse=", ")))
  
  missing_octa <- octa_vars[!octa_vars %in% colnames(original_data)]
  missing_wear <- wearable_vars[!wearable_vars %in% colnames(original_data)]
  
  if(length(missing_octa) > 0) {
    warning(paste("以下OCTA变量不在数据集中:", paste(missing_octa, collapse=", ")))
    octa_vars <- octa_vars[octa_vars %in% colnames(original_data)]
  }
  
  if(length(missing_wear) > 0) {
    warning(paste("以下可穿戴设备变量不在数据集中:", paste(missing_wear, collapse=", ")))
    wearable_vars <- wearable_vars[wearable_vars %in% colnames(original_data)]
  }
  
  if(length(octa_vars) == 0 || length(wearable_vars) == 0) {
    stop("没有有效的变量可以进行分析")
  }
  
  # 处理每个组和每个OCTA变量
  groups <- c("Diabetes", "Control")
  
  for(group in groups) {
    # 设置组过滤条件
    group_filter <- if(group == "Diabetes") {
      quote(dm_2 == 1)
    } else {
      quote(dm_2 == 0)
    }
    
    for(octa_var in octa_vars) {
      print(paste("处理", group, "组的", octa_var))
      
      # 执行LOOCV
      loocv_results <- perform_loocv(original_data, octa_var, wearable_vars, group_filter)
      
      # 检查是否有足够的结果
      if(!is.null(loocv_results) && nrow(loocv_results) >= 3) {
        n_points <- nrow(loocv_results)
        r_squared <- loocv_results$r_squared[1]  # 所有行都是相同的
        rmse <- loocv_results$rmse[1]            # 所有行都是相同的
        
        print(paste("  成功完成LOOCV，数据点数:", n_points, 
                    "R²:", round(r_squared, 3), 
                    "RMSE:", round(rmse, 3)))
        
        # 计算适当的坐标轴范围
        range_min <- min(min(loocv_results$predicted, na.rm = TRUE), 
                         min(loocv_results$actual, na.rm = TRUE))
        range_max <- max(max(loocv_results$predicted, na.rm = TRUE), 
                         max(loocv_results$actual, na.rm = TRUE))
        
        # 给范围添加一些余量
        buffer <- (range_max - range_min) * 0.1
        plot_min <- max(0, range_min - buffer)  # 确保不小于0
        plot_max <- range_max + buffer
        
        # 创建散点图
        p <- ggplot(loocv_results, aes(x = predicted, y = actual)) +
          geom_point(alpha = 0.7, color = ifelse(group == "Diabetes", "#D6604D", "#4393C3"), size = 3) +
          geom_smooth(method = "lm", color = "black", linetype = "solid", se = FALSE) +
          geom_segment(aes(x = plot_min, y = plot_min, 
                           xend = plot_max, yend = plot_max), 
                       linetype = "dashed") +
          scale_x_continuous(limits = c(plot_min, plot_max)) +
          scale_y_continuous(limits = c(plot_min, plot_max)) +
          labs(title = paste("LOOCV: Predicted vs Actual", octa_var),
               subtitle = paste(group, "Group -", n_points, "data points"),
               x = "Predicted Values (LOOCV)", y = "Actual Values") +
          annotate("text", x = plot_min + 0.05 * (plot_max - plot_min), 
                   y = plot_max - 0.05 * (plot_max - plot_min),
                   label = sprintf("LOOCV R² = %.3f\nLOOCV RMSE = %.3f", 
                                   r_squared, rmse),
                   hjust = 0, vjust = 1) +
          theme_custom +
          coord_fixed(ratio = 1)
        
        print(p)
        filename <- paste0("loocv_", var_type, "_", tolower(gsub(" ", "_", group)), "_", 
                           tolower(gsub("[^a-zA-Z0-9]", "_", octa_var)), "_pred_actual")
        ggsave(paste0(filename, ".png"), p, width = 7, height = 7, dpi = 300)
        ggsave(paste0(filename, ".pdf"), p, width = 7, height = 7, units = "in", device = pdf)
      } else {
        warning(paste("无法为", octa_var, "在", group, "组中创建LOOCV预测"))
      }
    }
  }
}

# 首先检查数据结构
print("检查可用的数据框和变量:")
print(paste("wearable_thickness_data行数:", nrow(wearable_thickness_data)))
print(paste("wearable_bloodflow_data行数:", nrow(wearable_bloodflow_data)))

# 获取变量名称
thickness_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")
bloodflow_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")

# 检查是否存在significant_thickness和significant_bloodflow数据框
if(exists("significant_thickness") && is.data.frame(significant_thickness)) {
  print("使用significant_thickness中的变量")
  thickness_octa_vars <- unique(significant_thickness$OCTA_Variable)
} else {
  print("significant_thickness不存在，使用替代方法确定OCTA变量")
  # 尝试从列名中查找OCTA变量
  possible_thickness_vars <- colnames(wearable_thickness_data)[
    !colnames(wearable_thickness_data) %in% 
      c("subject_id", "mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep", 
        "dm_2", "age", "gender", "bmi", "hypertension_2")
  ]
  thickness_octa_vars <- possible_thickness_vars[grep("Thickness|thickness", possible_thickness_vars)]
  if(length(thickness_octa_vars) == 0) {
    # 如果找不到包含"Thickness"的变量，使用一般方法
    thickness_octa_vars <- possible_thickness_vars[1:min(5, length(possible_thickness_vars))]
  }
}

if(exists("significant_bloodflow") && is.data.frame(significant_bloodflow)) {
  print("使用significant_bloodflow中的变量")
  bloodflow_octa_vars <- unique(significant_bloodflow$OCTA_Variable)
} else {
  print("significant_bloodflow不存在，使用替代方法确定OCTA变量")
  # 尝试从列名中查找血流变量
  possible_bloodflow_vars <- colnames(wearable_bloodflow_data)[
    !colnames(wearable_bloodflow_data) %in% 
      c("subject_id", "mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep", 
        "dm_2", "age", "gender", "bmi", "hypertension_2")
  ]
  bloodflow_octa_vars <- possible_bloodflow_vars[grep("Flow|flow|Vessel|vessel", possible_bloodflow_vars)]
  if(length(bloodflow_octa_vars) == 0) {
    # 如果找不到包含"Flow"的变量，使用一般方法
    bloodflow_octa_vars <- possible_bloodflow_vars[1:min(5, length(possible_bloodflow_vars))]
  }
}

print(paste("厚度OCTA变量:", paste(head(thickness_octa_vars, 3), collapse=", "), "..."))
print(paste("血流OCTA变量:", paste(head(bloodflow_octa_vars, 3), collapse=", "), "..."))

# 执行LOOCV并创建预测图
create_loocv_prediction_plots(wearable_thickness_data, thickness_octa_vars, thickness_vars, "thickness")
create_loocv_prediction_plots(wearable_bloodflow_data, bloodflow_octa_vars, bloodflow_vars, "bloodflow")
