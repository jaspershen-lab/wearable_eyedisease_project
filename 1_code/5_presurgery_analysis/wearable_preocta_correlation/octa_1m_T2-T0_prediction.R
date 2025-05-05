library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
library(lmerTest)
library(corrplot)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# 输入输出目录设置
input_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data_grouped"
octa_dir <- "3_data_analysis/5_presurgery_analysis/octa_change_prediction/t2_t0_selected_parameters"
output_dir <- "3_data_analysis/5_presurgery_analysis/octa_change_prediction/t2_t0_prediction_models"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义标准时间窗口
time_window_mapping <- list(
  "pre_3d" = "pre_3d",             # 术前3天
  "pre_3d_7d" = "pre_3d_7d",       # 术前3-7天
  "pre_7d_all" = "pre_7d_all",     # 术前7天及以上
  "post_7d" = "post_7d",           # 术后7天内
  "post_7d_30d" = "post_7d_30d",   # 术后7-30天
  "post_day23_30" = "post_day23_30", # 术后23-30天
  "post_day27_30" = "post_day27_30", # 术后27-30天
  "post_1w" = "post_1w",           # 术后1周
  "post_over_30d" = "post_over_30d", # 术后超过30天
  "pre_all" = "pre_all"            # 所有术前数据
)

# 定义时间窗口
time_windows <- names(time_window_mapping)

# 获取所有CSV文件
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# 将文件分类为糖尿病组手术眼和非糖尿病组非手术眼
d_surg1_files <- grep("D_Surg1\\.csv$", all_files, value = TRUE)
nod_surg0_files <- grep("NoD_Surg0\\.csv$", all_files, value = TRUE)

cat("找到的D_Surg1文件数量:", length(d_surg1_files), "\n")
cat("找到的NoD_Surg0文件数量:", length(nod_surg0_files), "\n")

# 读取OCTA变化数据
octa_changes <- read.csv(file.path(octa_dir, "selected_octa_changes.csv"))

# 显示OCTA变化数据中的变量
cat("OCTA变化数据中的变量:\n")
print(names(octa_changes))

# 定义要预测的OCTA参数
octa_params <- c(
  "PA_Choriocapillaris_change",
  "SVD_NerveFiber_change",
  "PA_PED_change",
  "PA_Vitreous_change",
  "Thickness_GCL.IPL_change",
  "Thickness_Retina_change",
  "Thickness_RNFL_change"
)

# 定义可穿戴设备变量
wearable_vars <- c("min_rhr_steps_1", "bo_mean", "steps_total", "total_sleep")

# 控制变量
control_vars <- c("age", "gender", "bmi", "season")

# 识别文件的时间窗口
get_time_window <- function(file_path) {
  file_name <- basename(file_path)
  
  # 检查特定模式，从最具体到最一般
  for (pattern in names(time_window_mapping)) {
    if (grepl(paste0("^", pattern, "_"), file_name)) {
      return(pattern)
    }
  }
  
  # 如果没有匹配到预定义窗口，返回NA
  cat("警告: 无法识别文件", file_name, "的时间窗口\n")
  return(NA)
}

# 定义一个函数处理每个组和每个OCTA参数
process_group_param <- function(files, group_name, octa_param) {
  # 创建数据框来存储每个时间窗口模型的性能指标
  performance_df <- data.frame(
    time_window = character(),
    octa_parameter = character(),
    r_squared = numeric(),
    adj_r_squared = numeric(),
    rmse = numeric(),
    hr_only_r2 = numeric(),
    bo_only_r2 = numeric(),
    steps_only_r2 = numeric(),
    sleep_only_r2 = numeric(),
    all_r2 = numeric(),
    sample_size = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 基于文件名跟踪我们已经处理过的窗口
  processed_windows <- character(0)
  
  # 对每个时间窗口构建模型
  cat(paste0("\n开始为", group_name, "的", octa_param, "构建预测模型...\n"))
  
  for (i in 1:length(files)) {
    file_path <- files[i]
    file_name <- basename(file_path)
    current_window <- get_time_window(file_path)
    
    # 调试信息
    cat("处理文件:", file_name, "\n")
    cat("  识别的时间窗口:", if(is.na(current_window)) "未识别" else current_window, "\n")
    
    if (is.na(current_window) || !current_window %in% time_windows) {
      cat("  跳过: 不在指定的时间窗口列表中\n")
      next
    }
    
    # 检查是否已处理过该窗口
    if (current_window %in% processed_windows) {
      cat("  跳过: 已经处理过", current_window, "窗口\n")
      next
    }
    
    # 记录已处理的窗口
    processed_windows <- c(processed_windows, current_window)
    
    # 读取时间窗口数据
    period_data <- tryCatch({
      read.csv(file_path)
    }, error = function(e) {
      cat(paste0("  错误: 无法读取文件: ", e$message, "\n"))
      return(NULL)
    })
    
    if (is.null(period_data)) next
    
    # 合并OCTA变化数据和时间窗口数据
    merged_data <- octa_changes %>%
      inner_join(period_data, by = c("ID" = "subject_id"))
    
    # 检查合并后的数据大小
    cat(paste0("  合并后的数据行数: ", nrow(merged_data), "\n"))
    
    if (nrow(merged_data) < 10) {
      cat(paste0("  警告: 时间窗口 ", current_window, " 的有效数据少于10行，可能导致模型不可靠\n"))
      next
    }
    
    # 检查OCTA参数是否存在
    if (!octa_param %in% names(merged_data)) {
      cat(paste0("  警告: 数据中不存在OCTA参数 ", octa_param, "\n"))
      next
    }
    
    # 检查预测变量是否存在
    missing_predictors <- wearable_vars[!wearable_vars %in% names(merged_data)]
    if (length(missing_predictors) > 0) {
      cat(paste0("  警告: 缺少以下预测变量: ", paste(missing_predictors, collapse=", "), "\n"))
      # 从wearable_vars中移除缺失的变量
      available_predictors <- setdiff(wearable_vars, missing_predictors)
      if (length(available_predictors) == 0) {
        cat("  错误: 没有可用的预测变量，跳过此时间窗口\n")
        next
      }
      cat(paste0("  使用可用的预测变量: ", paste(available_predictors, collapse=", "), "\n"))
    } else {
      available_predictors <- wearable_vars
    }
    
    # 检查控制变量
    available_controls <- intersect(control_vars, names(merged_data))
    if (length(available_controls) < length(control_vars)) {
      missing_controls <- setdiff(control_vars, available_controls)
      cat(paste0("  警告: 缺少以下控制变量: ", paste(missing_controls, collapse=", "), "\n"))
    }
    
    # 创建完整数据集，移除NA值
    model_vars <- c(octa_param, available_predictors, available_controls)
    model_data <- merged_data[, model_vars]
    model_data <- model_data[complete.cases(model_data), ]
    
    # 检查最终数据集大小
    if (nrow(model_data) < 10) {
      cat(paste0("  警告: 移除NA值后数据不足10行，跳过此时间窗口\n"))
      next
    }
    
    cat(paste0("  最终模型数据行数: ", nrow(model_data), "\n"))
    
    # 转换分类变量为因子
    if ("gender" %in% names(model_data)) model_data$gender <- factor(model_data$gender)
    if ("season" %in% names(model_data)) model_data$season <- factor(model_data$season)
    
    # 设置随机种子确保可重复性
    set.seed(123)
    
    # 构建完整模型
    full_model_formula <- as.formula(paste0(octa_param, " ~ ", paste(c(available_predictors, available_controls), collapse=" + ")))
    cat(paste0("  完整模型公式: ", deparse(full_model_formula), "\n"))
    
    full_model <- tryCatch({
      lm(full_model_formula, data = model_data)
    }, error = function(e) {
      cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
      return(NULL)
    })
    
    if (is.null(full_model)) next
    
    # 计算交叉验证R²的函数
    cv_r2 <- function(model_formula) {
      tryCatch({
        # 使用trainControl进行交叉验证
        train_control <- trainControl(
          method = "cv",
          number = 5,
          savePredictions = "final"
        )
        
        # 使用caret的train函数进行交叉验证
        set.seed(123)
        cv_model <- train(
          model_formula,
          data = model_data,
          method = "lm",
          trControl = train_control
        )
        
        # 返回交叉验证R²
        return(cv_model$results$Rsquared)
      }, error = function(e) {
        cat(paste0("    交叉验证错误: ", e$message, "\n"))
        # 如果交叉验证失败，尝试拟合常规线性模型并返回R²
        regular_model <- lm(model_formula, data = model_data)
        return(summary(regular_model)$r.squared)
      })
    }
    
    # 为每个可用的预测变量单独计算模型
    cat("  计算各预测变量的单独效果...\n")
    predictor_r2 <- list()
    
    # 计算完整模型的交叉验证R²
    all_r2 <- cv_r2(full_model_formula)
    cat(paste0("    完整模型交叉验证R²: ", round(all_r2, 4), "\n"))
    
    # 为每个预测变量计算单独模型
    for (pred in available_predictors) {
      pred_formula <- as.formula(paste0(octa_param, " ~ ", pred, " + ", paste(available_controls, collapse=" + ")))
      pred_r2 <- cv_r2(pred_formula)
      predictor_r2[[pred]] <- pred_r2
      cat(paste0("    ", pred, " 单独预测R²: ", round(pred_r2, 4), "\n"))
    }
    
    # 将性能指标添加到数据框
    perf_row <- data.frame(
      time_window = current_window,
      octa_parameter = octa_param,
      r_squared = summary(full_model)$r.squared,
      adj_r_squared = summary(full_model)$adj.r.squared,
      rmse = sqrt(mean(full_model$residuals^2)),
      sample_size = nrow(model_data),
      all_r2 = all_r2,
      stringsAsFactors = FALSE
    )
    
    # 添加每个预测变量的R²
    for (pred in wearable_vars) {
      if (pred %in% names(predictor_r2)) {
        perf_row[[paste0(gsub("_", "", pred), "_only_r2")]] <- predictor_r2[[pred]]
      } else {
        perf_row[[paste0(gsub("_", "", pred), "_only_r2")]] <- NA
      }
    }
    
    # 补充可能缺失的列
    missing_cols <- setdiff(
      c("hr_only_r2", "bo_only_r2", "steps_only_r2", "sleep_only_r2"),
      names(perf_row)
    )
    for (col in missing_cols) {
      perf_row[[col]] <- NA
    }
    
    # 添加到性能数据框
    performance_df <- rbind(performance_df, perf_row)
    
    cat(paste0("  完成时间窗口 ", current_window, " 的模型。R²: ", round(summary(full_model)$r.squared, 4), 
               ", 交叉验证R²: ", round(all_r2, 4), "\n"))
  }
  
  # 排序性能数据框，按照time_windows中的顺序
  if (nrow(performance_df) > 0) {
    performance_df$time_window <- factor(performance_df$time_window, levels = time_windows)
    performance_df <- performance_df[order(performance_df$time_window), ]
  }
  
  return(performance_df)
}

# 存储所有OCTA参数的模型性能
all_performance_d <- data.frame()
all_performance_nod <- data.frame()

# 对每个OCTA参数处理
for (octa_param in octa_params) {
  # 处理D_Surg1组
  cat(paste0("\n=== 处理 D_Surg1 组的 ", octa_param, " ===\n"))
  perf_d <- process_group_param(d_surg1_files, "D_Surg1", octa_param)
  if (nrow(perf_d) > 0) {
    all_performance_d <- rbind(all_performance_d, perf_d)
    
    # 为当前OCTA参数保存单独的结果
    write.csv(perf_d, file.path(output_dir, paste0("model_performance_D_Surg1_", octa_param, ".csv")), row.names = FALSE)
  }
  
  # 处理NoD_Surg0组
  cat(paste0("\n=== 处理 NoD_Surg0 组的 ", octa_param, " ===\n"))
  perf_nod <- process_group_param(nod_surg0_files, "NoD_Surg0", octa_param)
  if (nrow(perf_nod) > 0) {
    all_performance_nod <- rbind(all_performance_nod, perf_nod)
    
    # 为当前OCTA参数保存单独的结果
    write.csv(perf_nod, file.path(output_dir, paste0("model_performance_NoD_Surg0_", octa_param, ".csv")), row.names = FALSE)
  }
}

# 保存所有结果
write.csv(all_performance_d, file.path(output_dir, "all_octa_model_performance_D_Surg1.csv"), row.names = FALSE)
write.csv(all_performance_nod, file.path(output_dir, "all_octa_model_performance_NoD_Surg0.csv"), row.names = FALSE)


# 为每个OCTA参数创建可视化
for (octa_param in octa_params) {
  # D_Surg1组的可视化
  d_param_data <- all_performance_d[all_performance_d$octa_parameter == octa_param,]
  
  if(nrow(d_param_data) > 0) {
    # 创建绘图数据
    plot_data_d <- data.frame(
      time_window = rep(d_param_data$time_window, 5),
      predictor = c(
        rep("All Predictors", nrow(d_param_data)),
        rep("RHR", nrow(d_param_data)),
        rep("BO", nrow(d_param_data)),
        rep("Steps", nrow(d_param_data)),
        rep("Sleep", nrow(d_param_data))
      ),
      r2 = c(
        d_param_data$all_r2,
        d_param_data$hr_only_r2,
        d_param_data$bo_only_r2,
        d_param_data$steps_only_r2,
        d_param_data$sleep_only_r2
      )
    )
    
    # 绘制图形
    p_d <- ggplot(plot_data_d, aes(x = time_window, y = r2, color = predictor, group = predictor)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = c("All Predictors" = "#513a52", 
                                    "RHR" = "#82a1bf", 
                                    "BO" = "#faaa93", 
                                    "Steps" = "#feefc4",
                                    "Sleep" = "#77dd77")) +
      labs(
        title = paste0("Predictor Performance for ", gsub("_change", "", octa_param), " - D_Surg1 Group"),
        x = "Time window relative to surgery",
        y = "Cross-validated R²",
        color = "Predictors"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95")
      )
    
    # 保存图形
    ggsave(file.path(output_dir, paste0("predictor_performance_D_Surg1_", octa_param, ".png")), 
           p_d, width = 10, height = 6, dpi = 300)
  }
  
  # NoD_Surg0组的可视化
  nod_param_data <- all_performance_nod[all_performance_nod$octa_parameter == octa_param,]
  
  if(nrow(nod_param_data) > 0) {
    # 创建绘图数据
    plot_data_nod <- data.frame(
      time_window = rep(nod_param_data$time_window, 5),
      predictor = c(
        rep("All Predictors", nrow(nod_param_data)),
        rep("RHR", nrow(nod_param_data)),
        rep("BO", nrow(nod_param_data)),
        rep("Steps", nrow(nod_param_data)),
        rep("Sleep", nrow(nod_param_data))
      ),
      r2 = c(
        nod_param_data$all_r2,
        nod_param_data$hr_only_r2,
        nod_param_data$bo_only_r2,
        nod_param_data$steps_only_r2,
        nod_param_data$sleep_only_r2
      )
    )
    
    # 绘制图形
    p_nod <- ggplot(plot_data_nod, aes(x = time_window, y = r2, color = predictor, group = predictor)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = c("All Predictors" = "#513a52", 
                                    "RHR" = "#82a1bf", 
                                    "BO" = "#faaa93", 
                                    "Steps" = "#feefc4",
                                    "Sleep" = "#77dd77")) +
      labs(
        title = paste0("Predictor Performance for ", gsub("_change", "", octa_param), " - NoD_Surg0 Group"),
        x = "Time window relative to surgery",
        y = "Cross-validated R²",
        color = "Predictors"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95")
      )
    
    # 保存图形
    ggsave(file.path(output_dir, paste0("predictor_performance_NoD_Surg0_", octa_param, ".png")), 
           p_nod, width = 10, height = 6, dpi = 300)
  }
}

# 创建热图展示不同时间窗口的预测性能
# D_Surg1组热图
if(nrow(all_performance_d) > 0) {
  heatmap_data_d <- all_performance_d %>%
    mutate(octa_parameter = gsub("_change", "", octa_parameter)) %>%
    dplyr::select(time_window, octa_parameter, all_r2)
  
  # 转换为宽格式
  heatmap_data_d_wide <- heatmap_data_d %>%
    pivot_wider(names_from = time_window, values_from = all_r2)
  
  # 转换回长格式但保持时间窗口顺序
  heatmap_data_d <- heatmap_data_d_wide %>%
    pivot_longer(cols = -octa_parameter, 
                 names_to = "time_window", 
                 values_to = "all_r2") %>%
    mutate(time_window = factor(time_window, levels = time_windows))
  
  # 创建热图
  p_heatmap_d <- ggplot(heatmap_data_d, aes(x = time_window, y = octa_parameter, fill = all_r2)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", mid = "#bce4d8", high = "#02818a", 
                         midpoint = 0.3, na.value = "gray90", limits = c(0, 1),
                         name = "R²") +
    geom_text(aes(label = ifelse(is.na(all_r2), "", sprintf("%.2f", all_r2))), 
              color = "black", size = 3) +
    labs(
      title = "Predictive Performance for OCTA Parameters - D_Surg1 Group",
      x = "Time Window",
      y = "OCTA Parameter",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 保存热图
  ggsave(file.path(output_dir, "heatmap_prediction_performance_D_Surg1.png"), 
         p_heatmap_d, width = 12, height = 8, dpi = 300)
}

# NoD_Surg0组热图
if(nrow(all_performance_nod) > 0) {
  heatmap_data_nod <- all_performance_nod %>%
    mutate(octa_parameter = gsub("_change", "", octa_parameter)) %>%
    dplyr::select(time_window, octa_parameter, all_r2)
  
  # 转换为宽格式
  heatmap_data_nod_wide <- heatmap_data_nod %>%
    pivot_wider(names_from = time_window, values_from = all_r2)
  
  # 转换回长格式但保持时间窗口顺序
  heatmap_data_nod <- heatmap_data_nod_wide %>%
    pivot_longer(cols = -octa_parameter, 
                 names_to = "time_window", 
                 values_to = "all_r2") %>%
    mutate(time_window = factor(time_window, levels = time_windows))
  
  # 创建热图
  p_heatmap_nod <- ggplot(heatmap_data_nod, aes(x = time_window, y = octa_parameter, fill = all_r2)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", mid = "#bce4d8", high = "#02818a", 
                         midpoint = 0.3, na.value = "gray90", limits = c(0, 1),
                         name = "R²") +
    geom_text(aes(label = ifelse(is.na(all_r2), "", sprintf("%.2f", all_r2))), 
              color = "black", size = 3) +
    labs(
      title = "Predictive Performance for OCTA Parameters - NoD_Surg0 Group",
      x = "Time Window",
      y = "OCTA Parameter",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 保存热图
  ggsave(file.path(output_dir, "heatmap_prediction_performance_NoD_Surg0.png"), 
         p_heatmap_nod, width = 12, height = 8, dpi = 300)
}
