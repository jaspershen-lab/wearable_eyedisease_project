# -----------------------------------------------------
# 分析不同步数模式聚类对OCTA改善的影响
# -----------------------------------------------------
library(tidyverse)
library(ggplot2)
library(rstatix)      # 用于统计测试
library(ggpubr)       # 用于创建出版质量图表
library(ggstatsplot)  # 用于带统计信息的图表
library(gtsummary)    # 用于创建汇总表格
library(broom)        # 用于整理模型结果
library(pROC)         # 用于ROC分析
library(forestplot)   # 用于森林图
library(gridExtra)    # 用于多图布局

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# -----------------------------------------------------
# 1. 加载必要的数据
# -----------------------------------------------------
# 加载糖尿病组(Group2)的步数聚类结果
# 可以根据需求选择带或不带SH012的聚类结果
group2_clusters <- readRDS("3_data_analysis/6_clustering_modeling/kshape/time_series/1m/steps/Group2_Diabetes_all_clusters.rds")
# group2_clusters <- readRDS("3_data_analysis/6_clustering_modeling/kshape/time_series/steps_1m_no_SH012/Group2_Diabetes_NoSH012_all_clusters.rds")

# 提取k=2的聚类结果(也可以替换为k=3的结果)
cluster_k2 <- group2_clusters$k2$clusters
cluster_assignments <- data.frame(
  subject_id = names(cluster_k2),
  steps_cluster = as.factor(cluster_k2)
)

# 加载OCTA数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 创建结果目录
# 可以根据使用的数据集选择适当的输出目录
# output_dir <- "3_data_analysis/6_clustering_modeling/kshape/octa_analysis/steps_1m/T2_T0"
output_dir <- "3_data_analysis/6_clustering_modeling/kshape/octa_analysis/1m/steps/T2_T0"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# 2. 处理基线和OCTA数据
# -----------------------------------------------------
# 处理OCTA数据函数
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
octa_bloodflow_processed <- process_all_patients(octa_bloodflow_features)
octa_thickness_processed <- process_all_patients(octa_thickness_features)

# -----------------------------------------------------
# 3. 计算OCTA指标的改善值(T2-T0)
# -----------------------------------------------------
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
bloodflow_improvement <- calculate_improvement(octa_bloodflow_processed, "血流")
thickness_improvement <- calculate_improvement(octa_thickness_processed, "厚度")

# -----------------------------------------------------
# 4. 过滤只包含目标OCTA层的参数
# -----------------------------------------------------
# 过滤血流参数函数
filter_bloodflow_layers <- function(data) {
  # 定义要关注的血流层
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  
  # 构建匹配模式
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*0_21_improvement$")
  
  # 获取匹配的列名
  params_of_interest <- names(data)[grep(pattern, names(data))]
  
  # 如果没有找到匹配的参数，输出警告
  if(length(params_of_interest) == 0) {
    warning("未找到指定的血流层参数！")
    return(list(data = data, params = character(0)))
  } else {
    cat("找到", length(params_of_interest), "个目标血流层的参数:\n")
    cat(paste(params_of_interest, collapse = "\n"), "\n")
  }
  
  # 保留ID和目标参数
  filtered_data <- data %>%
    dplyr::select(ID, all_of(params_of_interest))
  
  return(list(
    data = filtered_data,
    params = params_of_interest
  ))
}

# 过滤厚度参数函数
filter_thickness_layers <- function(data) {
  # 定义要关注的厚度层
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  
  # 构建匹配模式
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*0_21_improvement$")
  
  # 获取匹配的列名
  params_of_interest <- names(data)[grep(pattern, names(data))]
  
  # 如果没有找到匹配的参数，输出警告
  if(length(params_of_interest) == 0) {
    warning("未找到指定的厚度层参数！")
    return(list(data = data, params = character(0)))
  } else {
    cat("找到", length(params_of_interest), "个目标厚度层的参数:\n")
    cat(paste(params_of_interest, collapse = "\n"), "\n")
  }
  
  # 保留ID和目标参数
  filtered_data <- data %>%
    dplyr::select(ID, all_of(params_of_interest))
  
  return(list(
    data = filtered_data,
    params = params_of_interest
  ))
}

# 应用过滤器
bloodflow_filtered <- filter_bloodflow_layers(bloodflow_improvement)
thickness_filtered <- filter_thickness_layers(thickness_improvement)

# -----------------------------------------------------
# 5. 合并步数聚类结果和OCTA改善值
# -----------------------------------------------------
# 合并血流数据
bloodflow_with_clusters <- bloodflow_filtered$data %>%
  inner_join(cluster_assignments, by = c("ID" = "subject_id"))

# 合并厚度数据
thickness_with_clusters <- thickness_filtered$data %>%
  inner_join(cluster_assignments, by = c("ID" = "subject_id"))

# -----------------------------------------------------
# 6. 连续变量分析：不同步数模式聚类组之间的OCTA改善差异
# -----------------------------------------------------
# 创建统计分析函数
analyze_octa_by_cluster <- function(data, param_cols, output_prefix) {
  results_list <- list()
  all_p_values <- c()
  param_names <- c()
  
  # 创建汇总表
  summary_df <- data.frame(
    Parameter = character(),
    P_Value = numeric(),
    P_Adjusted = numeric(),
    Significant = logical(),
    Effect_Size = numeric(),
    Mean_Cluster1 = numeric(),
    Mean_Cluster2 = numeric(),
    Mean_Diff = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 为每个参数创建可视化和执行统计测试
  for (param in param_cols) {
    param_name <- gsub("_0_21_improvement", "", param)
    
    # 创建公式
    formula_str <- paste(param, "~ steps_cluster")
    
    # 执行ANOVA
    anova_result <- data %>%
      anova_test(as.formula(formula_str))
    
    # 保存p值
    p_value <- anova_result$p
    all_p_values <- c(all_p_values, p_value)
    param_names <- c(param_names, param_name)
    
    # 计算每个聚类的均值
    cluster_means <- data %>%
      group_by(steps_cluster) %>%
      summarize(
        mean_value = mean(!!sym(param), na.rm = TRUE),
        .groups = "drop"
      )
    
    # 提取均值
    mean_c1 <- cluster_means$mean_value[cluster_means$steps_cluster == 1]
    mean_c2 <- cluster_means$mean_value[cluster_means$steps_cluster == 2]
    mean_diff <- mean_c2 - mean_c1
    
    # 计算效应量 (Cohen's d)
    effect_size <- cohens_d(as.formula(formula_str), data = data)$effsize
    
    # 添加到汇总表
    summary_df <- rbind(summary_df, data.frame(
      Parameter = param_name,
      P_Value = p_value,
      P_Adjusted = NA,  # 将在之后填充
      Significant = NA, # 将在之后填充
      Effect_Size = effect_size,
      Mean_Cluster1 = mean_c1,
      Mean_Cluster2 = mean_c2,
      Mean_Diff = mean_diff,
      stringsAsFactors = FALSE
    ))
    
    # 创建可视化
    p <- ggstatsplot::ggbetweenstats(
      data = data,
      x = steps_cluster,
      y = !!sym(param),
      type = "parametric",
      pairwise.display = "all",
      p.adjust.method = "fdr",
      effsize.type = "unbiased",
      results.subtitle = TRUE,
      xlab = "步数聚类组",
      ylab = paste(param_name, "改善值 (T2-T0)"),
      title = paste0(param_name, " 在不同步数模式组间的差异"),
      centrality.plotting = TRUE,
      centrality.point.args = list(size = 5, color = "darkred"),
      centrality.type = "parametric",
      ggtheme = ggplot2::theme_bw()
    )
    
    # 保存图表
    ggsave(
      file.path(output_dir, paste0(output_prefix, "_", param_name, ".pdf")),
      p,
      width = 10,
      height = 8
    )
    
    # 保存结果
    results_list[[param]] <- list(
      plot = p,
      anova = anova_result
    )
  }
  
  # 多重比较校正
  summary_df$P_Adjusted <- p.adjust(summary_df$P_Value, method = "fdr")
  summary_df$Significant <- summary_df$P_Adjusted < 0.05
  
  # 按校正后的p值排序
  summary_df <- summary_df %>%
    arrange(P_Adjusted)
  
  # 保存汇总表
  write.csv(
    summary_df,
    file.path(output_dir, paste0(output_prefix, "_results_summary.csv")),
    row.names = FALSE
  )
  
  # 创建top显著参数的多面板图
  if (any(summary_df$Significant)) {
    # 选择显著的参数
    sig_params <- summary_df %>%
      filter(Significant) %>%
      arrange(P_Adjusted) %>%
      pull(Parameter)
    
    sig_params_original <- param_cols[match(sig_params, gsub("_0_21_improvement", "", param_cols))]
    
    # 取前4个或全部显著参数
    top_params <- sig_params_original[1:min(4, length(sig_params_original))]
    
    # 创建图表列表
    plot_list <- list()
    
    for (i in 1:length(top_params)) {
      param <- top_params[i]
      param_name <- gsub("_0_21_improvement", "", param)
      
      p <- ggstatsplot::ggbetweenstats(
        data = data,
        x = steps_cluster,
        y = !!sym(param),
        type = "parametric",
        pairwise.display = "all",
        p.adjust.method = "fdr",
        effsize.type = "unbiased",
        results.subtitle = TRUE,
        xlab = "步数聚类组",
        ylab = "改善值 (T2-T0)",
        title = param_name,
        centrality.plotting = TRUE,
        centrality.point.args = list(size = 5, color = "darkred"),
        centrality.type = "parametric",
        ggtheme = ggplot2::theme_bw()
      )
      
      plot_list[[i]] <- p
    }
    
    # 合并图表
    combined_plot <- ggpubr::ggarrange(
      plotlist = plot_list,
      ncol = 2,
      nrow = ceiling(length(plot_list) / 2),
      common.legend = TRUE,
      legend = "bottom"
    )
    
    # 添加总标题
    combined_plot <- ggpubr::annotate_figure(
      combined_plot,
      top = ggpubr::text_grob(
        paste0("显著的OCTA改善差异 - ", if(grepl("bloodflow", output_prefix)) "血流" else "厚度"),
        face = "bold",
        size = 16
      )
    )
    
    # 保存合并图表
    ggsave(
      file.path(output_dir, paste0(output_prefix, "_significant_params.pdf")),
      combined_plot,
      width = 14,
      height = 10
    )
  }
  
  return(list(
    results = results_list,
    summary = summary_df
  ))
}

# 执行血流参数连续变量分析
bloodflow_analysis <- analyze_octa_by_cluster(
  bloodflow_with_clusters,
  bloodflow_filtered$params,
  "diabetes_bloodflow_by_steps_cluster"
)

# 执行厚度参数连续变量分析
thickness_analysis <- analyze_octa_by_cluster(
  thickness_with_clusters,
  thickness_filtered$params,
  "diabetes_thickness_by_steps_cluster"
)

# -----------------------------------------------------
# 7. 二分类变量分析：基于中位数切点计算OR
# -----------------------------------------------------
# 创建二分类分析函数
analyze_octa_binary <- function(data, param_cols, output_prefix) {
  # 初始化结果
  binary_results <- data.frame(
    Parameter = character(),
    OR = numeric(),
    OR_Lower = numeric(),
    OR_Upper = numeric(),
    P_Value = numeric(),
    P_Adjusted = numeric(),
    Significant = logical(),
    Median_Cutoff = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 循环计算每个参数的OR值
  for (param in param_cols) {
    param_name <- gsub("_0_21_improvement", "", param)
    
    # # 计算中位数
    # median_val <- median(data[[param]], na.rm = TRUE)
    # 
    # # 创建二分类变量（1表示高于中位数的改善，0表示低于或等于中位数）
    # binary_var_name <- paste0(param, "_high")
    # data[[binary_var_name]] <- as.factor(ifelse(data[[param]] > median_val, 1, 0))
    
    # 修改后（基于0分类）:
    # 创建二分类变量（1表示正改善[>0]，0表示无改善或负改善[≤0]）
    binary_var_name <- paste0(param, "_positive")
    data[[binary_var_name]] <- as.factor(ifelse(data[[param]] > 0, 1, 0))
    
    # 创建contingency table
    cont_table <- table(data$steps_cluster, data[[binary_var_name]])
    
    # Fisher's exact test
    fisher_test <- fisher.test(cont_table)
    
    # 提取OR值及其置信区间
    or <- fisher_test$estimate
    or_lower <- fisher_test$conf.int[1]
    or_upper <- fisher_test$conf.int[2]
    p_value <- fisher_test$p.value
    
    # 保存结果
    binary_results <- rbind(binary_results, data.frame(
      Parameter = param_name,
      OR = or,
      OR_Lower = or_lower,
      OR_Upper = or_upper,
      P_Value = p_value,
      P_Adjusted = NA,  # 将在后面填充
      Significant = NA, # 将在后面填充
      # Median_Cutoff = median_val,
      stringsAsFactors = FALSE
    ))
    
    # 创建条形图
    cross_tab <- data %>%
      group_by(steps_cluster, !!sym(binary_var_name)) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        percentage = count / sum(count) * 100,
        improvement = ifelse(!!sym(binary_var_name) == 1, "高改善", "低改善")
      )
    
    p <- ggplot(cross_tab, aes(x = steps_cluster, y = percentage, fill = improvement)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = sprintf("%.1f%%", percentage)), 
                position = position_dodge(width = 0.9), vjust = -0.5) +
      labs(
        title = paste0(param_name, " - 不同步数聚类组中高/低改善比例"),
        subtitle = paste0("OR = ", round(or, 2), " (95% CI: ", 
                          round(or_lower, 2), "-", round(or_upper, 2), 
                          "), p = ", format.pval(p_value, digits = 3)),
        x = "步数聚类组",
        y = "百分比 (%)",
        fill = "改善水平"
      ) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1")
    
    # 保存图表
    ggsave(
      file.path(output_dir, paste0(output_prefix, "_binary_", param_name, ".pdf")),
      p,
      width = 8,
      height = 6
    )
  }
  
  # 多重比较校正
  binary_results$P_Adjusted <- p.adjust(binary_results$P_Value, method = "fdr")
  binary_results$Significant <- binary_results$P_Adjusted < 0.05
  
  # 按校正后的p值排序
  binary_results <- binary_results %>%
    arrange(P_Adjusted)
  
  # 保存汇总表
  write.csv(
    binary_results,
    file.path(output_dir, paste0(output_prefix, "_binary_results.csv")),
    row.names = FALSE
  )
  
  # 创建森林图
  if (nrow(binary_results) > 0) {
    # 准备森林图数据
    forest_data <- binary_results %>%
      mutate(
        OR_CI = paste0(format(round(OR, 2), nsmall = 2), " (", 
                       format(round(OR_Lower, 2), nsmall = 2), "-", 
                       format(round(OR_Upper, 2), nsmall = 2), ")"),
        P_Value_Display = ifelse(P_Adjusted < 0.001, "< 0.001", 
                                 format(round(P_Adjusted, 3), nsmall = 3))
      ) %>%
      arrange(P_Adjusted)
    
    # 取前10个参数（如果超过10个）
    if(nrow(forest_data) > 10) {
      forest_data <- forest_data[1:10,]
    }
    
    # 准备森林图
    forest_matrix <- cbind(
      c("参数", forest_data$Parameter),
      c("OR (95% CI)", forest_data$OR_CI),
      c("p值", forest_data$P_Value_Display)
    )
    
    # 森林图
    pdf(file.path(output_dir, paste0(output_prefix, "_forest_plot.pdf")), 
        width = 10, height = 8)
    
    forestplot(
      labeltext = forest_matrix,
      mean = c(NA, forest_data$OR),
      lower = c(NA, forest_data$OR_Lower),
      upper = c(NA, forest_data$OR_Upper),
      is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
      graph.pos = 2,
      hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:3, col = "#444444")),
      graphwidth = unit(50, "mm"),
      colgap = unit(5, "mm"),
      col = fpColors(box = "royalblue", lines = "darkblue", zero = "gray50"),
      boxsize = 0.25,
      line.margin = unit(5, "mm"),
      txt_gp = fpTxtGp(
        ticks = gpar(cex = 0.8),
        xlab = gpar(cex = 0.8),
        label = gpar(cex = 0.9)
      ),
      lwd.ci = 2,
      lwd.zero = 1,
      grid = TRUE,
      zero = 1,
      clip = c(0.1, 5),
      xticks = c(0.1, 0.5, 1, 2, 5),
      xlab = "Odds Ratio (95% CI)",
      title = paste0("步数聚类对OCTA改善的影响 - ", 
                     if(grepl("bloodflow", output_prefix)) "血流参数" else "厚度参数")
    )
    
    dev.off()
  }
  
  return(binary_results)
}

# 执行血流参数二分类分析
bloodflow_binary_analysis <- analyze_octa_binary(
  bloodflow_with_clusters,
  bloodflow_filtered$params,
  "diabetes_bloodflow_by_steps_cluster"
)

# 执行厚度参数二分类分析
thickness_binary_analysis <- analyze_octa_binary(
  thickness_with_clusters,
  thickness_filtered$params,
  "diabetes_thickness_by_steps_cluster"
)

# -----------------------------------------------------
# 8. 生成综合报告
# -----------------------------------------------------
# 创建HTML报告函数
create_html_report <- function() {
  # 血流连续变量分析摘要
  bloodflow_summary <- bloodflow_analysis$summary
  bloodflow_sig <- bloodflow_summary %>% filter(Significant)
  
  # 厚度连续变量分析摘要
  thickness_summary <- thickness_analysis$summary
  thickness_sig <- thickness_summary %>% filter(Significant)
  
  # 血流二分类分析摘要
  bloodflow_binary_sig <- bloodflow_binary_analysis %>% filter(Significant)
  
  # 厚度二分类分析摘要
  thickness_binary_sig <- thickness_binary_analysis %>% filter(Significant)
  
  # 创建HTML报告
  sink(file.path(output_dir, "steps_cluster_octa_analysis_report.html"))
  
  cat('<!DOCTYPE html>
<html>
<head>
  <title>步数模式聚类与OCTA改善关系分析报告</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #2980b9; margin-top: 30px; }
    h3 { color: #3498db; }
    table { border-collapse: collapse; width: 100%; margin: 20px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #f2f2f2; }
    tr:nth-child(even) { background-color: #f9f9f9; }
    .significant { color: #e74c3c; font-weight: bold; }
    .non-significant { color: #7f8c8d; }
    .summary { background-color: #eef7fa; padding: 15px; border-radius: 5px; border-left: 5px solid #3498db; }
  </style>
</head>
<body>
  <h1>步数模式聚类与OCTA改善关系分析报告</h1>
  
  <div class="summary">
    <p>本报告分析了糖尿病患者(Group2)的不同步数恢复模式聚类组之间OCTA参数改善值的差异。分析包括连续变量分析和基于中位数的二分类变量分析。</p>
  </div>
  
  <h2>1. 连续变量分析结果</h2>')
  
  # 血流参数结果
  cat('<h3>1.1 血流参数</h3>')
  
  if(nrow(bloodflow_sig) > 0) {
    cat('<p>显著的血流参数差异:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>P值</th><th>校正后P值</th><th>效应量</th><th>聚类1均值</th><th>聚类2均值</th><th>差异</th></tr>')
    
    for(i in 1:nrow(bloodflow_sig)) {
      cat('<tr>')
      cat('<td>', bloodflow_sig$Parameter[i], '</td>')
      cat('<td>', format(bloodflow_sig$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(bloodflow_sig$P_Adjusted[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_sig$Effect_Size[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_sig$Mean_Cluster1[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_sig$Mean_Cluster2[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_sig$Mean_Diff[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的血流参数差异。</p>')
  }
  
  # 厚度参数结果
  cat('<h3>1.2 厚度参数</h3>')
  
  if(nrow(thickness_sig) > 0) {
    cat('<p>显著的厚度参数差异:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>P值</th><th>校正后P值</th><th>效应量</th><th>聚类1均值</th><th>聚类2均值</th><th>差异</th></tr>')
    
    for(i in 1:nrow(thickness_sig)) {
      cat('<tr>')
      cat('<td>', thickness_sig$Parameter[i], '</td>')
      cat('<td>', format(thickness_sig$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(thickness_sig$P_Adjusted[i], digits=3), '</td>')
      cat('<td>', format(thickness_sig$Effect_Size[i], digits=3), '</td>')
      cat('<td>', format(thickness_sig$Mean_Cluster1[i], digits=3), '</td>')
      cat('<td>', format(thickness_sig$Mean_Cluster2[i], digits=3), '</td>')
      cat('<td>', format(thickness_sig$Mean_Diff[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的厚度参数差异。</p>')
  }
  
  # 二分类变量分析结果
  cat('<h2>2. 二分类变量分析结果 (基于中位数)</h2>')
  
  # 血流参数二分类结果
  cat('<h3>2.1 血流参数</h3>')
  
  if(nrow(bloodflow_binary_sig) > 0) {
    cat('<p>显著的血流参数OR值:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(bloodflow_binary_sig)) {
      cat('<tr>')
      cat('<td>', bloodflow_binary_sig$Parameter[i], '</td>')
      cat('<td>', format(bloodflow_binary_sig$OR[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_sig$OR_Lower[i], digits=3), '-', 
          format(bloodflow_binary_sig$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_sig$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(bloodflow_binary_sig$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的血流参数OR值差异。</p>')
  }
  
  # 厚度参数二分类结果
  cat('<h3>2.2 厚度参数</h3>')
  
  if(nrow(thickness_binary_sig) > 0) {
    cat('<p>显著的厚度参数OR值:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(thickness_binary_sig)) {
      cat('<tr>')
      cat('<td>', thickness_binary_sig$Parameter[i], '</td>')
      cat('<td>', format(thickness_binary_sig$OR[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_sig$OR_Lower[i], digits=3), '-', 
          format(thickness_binary_sig$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_sig$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(thickness_binary_sig$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的厚度参数OR值差异。</p>')
  }
  
  # 结论
  cat('<h2>3. 结论</h2>')
  
  total_sig <- nrow(bloodflow_sig) + nrow(thickness_sig) + 
    nrow(bloodflow_binary_sig) + nrow(thickness_binary_sig)
  
  if(total_sig > 0) {
    cat('<p>本分析发现了步数模式聚类组之间存在显著的OCTA参数改善差异，具体如下：</p>')
    cat('<ul>')
    
    if(nrow(bloodflow_sig) > 0) {
      cat('<li>血流参数连续变量分析：发现', nrow(bloodflow_sig), '个显著差异参数</li>')
    }
    
    if(nrow(thickness_sig) > 0) {
      cat('<li>厚度参数连续变量分析：发现', nrow(thickness_sig), '个显著差异参数</li>')
    }
    
    if(nrow(bloodflow_binary_sig) > 0) {
      cat('<li>血流参数二分类分析：发现', nrow(bloodflow_binary_sig), '个显著OR值</li>')
    }
    
    if(nrow(thickness_binary_sig) > 0) {
      cat('<li>厚度参数二分类分析：发现', nrow(thickness_binary_sig), '个显著OR值</li>')
    }
    
    cat('</ul>')
    
    cat('<p>这些结果表明不同的步数恢复模式可能与术后OCTA参数的改善情况相关。建议进一步研究这些关联的临床意义。</p>')
  } else {
    cat('<p>本分析未发现步数模式聚类组之间存在显著的OCTA参数改善差异。这可能表明：</p>')
    cat('<ul>')
    cat('<li>步数恢复模式与OCTA参数改善可能没有直接关系</li>')
    cat('<li>可能需要更多样本量或更长的随访时间来检测可能的差异</li>')
    cat('<li>可能需要考虑其他潜在的混杂因素</li>')
    cat('</ul>')
  }
  
  cat('<p>建议结合患者的临床特征和其他检查结果进行综合分析。</p>')
  
  cat('</body>
</html>')
  
  sink()
  
  cat("HTML报告已生成：", file.path(output_dir, "steps_cluster_octa_analysis_report.html"), "\n")
}

# 生成HTML报告
create_html_report()

# 打印分析摘要
cat("\n=====================================================\n")
cat("步数模式聚类与OCTA改善关系分析完成\n")
cat("=====================================================\n\n")

cat("分析结果已保存到目录：", output_dir, "\n\n")

cat("血流参数连续变量分析显著结果数量：", sum(bloodflow_analysis$summary$Significant), "\n")
cat("厚度参数连续变量分析显著结果数量：", sum(thickness_analysis$summary$Significant), "\n")
cat("血流参数二分类分析显著结果数量：", sum(bloodflow_binary_analysis$Significant), "\n")
cat("厚度参数二分类分析显著结果数量：", sum(thickness_binary_analysis$Significant), "\n\n")

