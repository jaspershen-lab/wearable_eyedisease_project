# -----------------------------------------------------
# 分析Group1和Group2步数模式聚类对视力改善的影响
# 仅包含1周和1个月视力数据
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
# 加载Group1和Group2的步数聚类结果
group1_clusters <- readRDS("3_data_analysis/6_clustering_modeling/kshape/time_series/1m/steps/Group1_Surgery0_all_clusters.rds")
group2_clusters <- readRDS("3_data_analysis/6_clustering_modeling/kshape/time_series/1m/steps/Group2_Diabetes_all_clusters.rds")

# 提取k=2的聚类结果(也可以替换为k=3的结果)
# Group1
cluster1_k2 <- group1_clusters$k2$clusters
cluster1_assignments <- data.frame(
  subject_id = names(cluster1_k2),
  steps_cluster = as.factor(cluster1_k2)
)

# Group2
cluster2_k2 <- group2_clusters$k2$clusters
cluster2_assignments <- data.frame(
  subject_id = names(cluster2_k2),
  steps_cluster = as.factor(cluster2_k2)
)

# 加载基线信息和视力数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 创建结果目录
output_dir_group1 <- "3_data_analysis/6_clustering_modeling/kshape/vision_analysis/1m/steps/Group1"
output_dir_group2 <- "3_data_analysis/6_clustering_modeling/kshape/vision_analysis/1m/steps/Group2"
dir.create(output_dir_group1, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_group2, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# 2. 准备视力数据，计算改善值
# -----------------------------------------------------
# 创建视力数据集，包括术前和术后(1周和1个月)的数据
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,  # 右眼手术
      surgery_eye_1 == 1 ~ os_corrected_bas,  # 左眼手术
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # 双眼(平均)
      TRUE ~ NA_real_
    ),
    post_vision_1w = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1w,   # 右眼术后1周
      surgery_eye_1 == 1 ~ os_corrected_1w,   # 左眼术后1周
      surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # 双眼术后1周
      TRUE ~ NA_real_
    ),
    post_vision_1m = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1m,   # 右眼术后1个月
      surgery_eye_1 == 1 ~ os_corrected_1m,   # 左眼术后1个月
      surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,   # 双眼术后1个月
      TRUE ~ NA_real_
    ),
    vision_improvement_1w = post_vision_1w - pre_vision,
    vision_improvement_1m = post_vision_1m - pre_vision,
    vision_improved_1m = ifelse(vision_improvement_1m > 0, 1, 0),  # 二分类改善指标
    vision_improved_factor_1m = factor(vision_improved_1m, 
                                       levels = c(0, 1), 
                                       labels = c("无改善", "有改善"))  # 分类变量版本
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision_1w, post_vision_1m,
                vision_improvement_1w, vision_improvement_1m, 
                vision_improved_1m, vision_improved_factor_1m,
                age, gender)

# -----------------------------------------------------
# 3. 合并步数聚类结果和视力改善数据
# -----------------------------------------------------
# 分别合并Group1和Group2的数据
vision_with_clusters_group1 <- vision_data %>%
  inner_join(cluster1_assignments, by = c("ID" = "subject_id"))

vision_with_clusters_group2 <- vision_data %>%
  inner_join(cluster2_assignments, by = c("ID" = "subject_id"))

# -----------------------------------------------------
# 4. 连续变量分析：不同步数模式聚类组之间的视力改善差异
# -----------------------------------------------------
# 统计分析函数
analyze_vision_by_cluster <- function(data, vision_params, output_prefix, output_dir, group_name) {
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
  for (param in vision_params) {
    param_name <- gsub("vision_improvement_", "", param)
    
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
      ylab = paste(param_name, "改善值"),
      title = paste0(group_name, ": ", param_name, " 在不同步数模式组间的差异"),
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
    
    sig_params_original <- vision_params[match(sig_params, gsub("vision_improvement_", "", vision_params))]
    
    # 取前4个或全部显著参数
    top_params <- sig_params_original[1:min(4, length(sig_params_original))]
    
    # 创建图表列表
    plot_list <- list()
    
    for (i in 1:length(top_params)) {
      param <- top_params[i]
      param_name <- gsub("vision_improvement_", "", param)
      
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
        ylab = "视力改善值",
        title = paste(group_name, "- 时间点:", param_name),
        centrality.plotting = TRUE,
        centrality.point.args = list(size = 5, color = "darkred"),
        centrality.type = "parametric",
        ggtheme = ggplot2::theme_bw()
      )
      
      plot_list[[i]] <- p
    }
    
    # 合并图表
    if(length(plot_list) > 0) {
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
          paste(group_name, "- 显著的视力改善差异"),
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
  }
  
  return(list(
    results = results_list,
    summary = summary_df
  ))
}

# 定义视力改善参数（仅1周和1个月）
vision_params <- c("vision_improvement_1w", "vision_improvement_1m")

# 分别执行各组的视力参数连续变量分析
vision_analysis_group1 <- analyze_vision_by_cluster(
  vision_with_clusters_group1,
  vision_params,
  "vision_by_steps_cluster",
  output_dir_group1,
  "Group1"
)

vision_analysis_group2 <- analyze_vision_by_cluster(
  vision_with_clusters_group2,
  vision_params,
  "vision_by_steps_cluster",
  output_dir_group2,
  "Group2"
)

# -----------------------------------------------------
# 5. 二分类变量分析：基于改善与否计算OR
# -----------------------------------------------------
# 二分类分析函数
analyze_vision_binary <- function(data, vision_params, output_prefix, output_dir, group_name) {
  # 初始化结果
  binary_results <- data.frame(
    Parameter = character(),
    OR = numeric(),
    OR_Lower = numeric(),
    OR_Upper = numeric(),
    P_Value = numeric(),
    P_Adjusted = numeric(),
    Significant = logical(),
    Definition = character(),
    stringsAsFactors = FALSE
  )
  
  # 循环计算每个参数的OR值
  for (param in vision_params) {
    param_name <- gsub("vision_improvement_", "", param)
    
    # 尝试多种二分类方法
    binary_methods <- list(
      positive = list(
        name = "改善(>0)",
        func = function(x) ifelse(x > 0, 1, 0)
      ),
      median = list(
        name = "高于中位数",
        func = function(x) {
          med <- median(x, na.rm = TRUE)
          ifelse(x > med, 1, 0)
        }
      ),
      quartile = list(
        name = "高于75%分位数",
        func = function(x) {
          q3 <- quantile(x, 0.75, na.rm = TRUE)
          ifelse(x > q3, 1, 0)
        }
      )
    )
    
    # 尝试每种二分类方法
    for (method_name in names(binary_methods)) {
      method <- binary_methods[[method_name]]
      # 创建二分类变量
      binary_var_name <- paste0(param, "_", method_name)
      data[[binary_var_name]] <- as.factor(method$func(data[[param]]))
      
      # 检查分类是否有效（每个类别至少有一个值）
      if (length(unique(data[[binary_var_name]])) < 2) {
        cat(paste0("警告: ", param_name, " 使用 ", method$name, " 方法分类无效，所有值都在同一类别\n"))
        next
      }
      
      # 创建列联表
      cont_table <- table(data$steps_cluster, data[[binary_var_name]])
      
      # 检查列联表是否有足够维度
      if (nrow(cont_table) < 2 || ncol(cont_table) < 2) {
        cat(paste0("警告: ", param_name, " 使用 ", method$name, " 方法的列联表维度不足\n"))
        next
      }
      
      # 检查每个单元格是否有足够样本
      min_cell <- min(cont_table)
      if (min_cell < 1) {
        cat(paste0("警告: ", param_name, " 使用 ", method$name, " 方法的列联表有空单元格\n"))
        next
      }
      
      # 执行Fisher's exact test
      tryCatch({
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
          Definition = method$name,
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
            title = paste0(group_name, ": ", param_name, " - 不同步数聚类组中视力改善比例"),
            subtitle = paste0(method$name, " OR = ", round(or, 2), " (95% CI: ", 
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
          file.path(output_dir, paste0(output_prefix, "_binary_", param_name, "_", method_name, ".pdf")),
          p,
          width = 8,
          height = 6
        )
      }, error = function(e) {
        cat(paste0("错误: ", param_name, " 使用 ", method$name, " 方法执行Fisher检验时出错: ", e$message, "\n"))
      })
    }
  }
  
  # 如果没有有效结果，返回空数据框
  if (nrow(binary_results) == 0) {
    cat(paste0("警告：", group_name, "没有找到任何有效的二分类分析结果\n"))
    return(binary_results)
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
        Parameter_Full = paste0(Parameter, " (", Definition, ")"),
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
      c("参数", forest_data$Parameter_Full),
      c("OR (95% CI)", forest_data$OR_CI),
      c("p值", forest_data$P_Value_Display)
    )
    
    # 森林图
    pdf(file.path(output_dir, paste0(output_prefix, "_forest_plot.pdf")), 
        width = 12, height = 10)
    
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
      title = paste(group_name, "- 步数聚类对视力改善的影响")
    )
    
    dev.off()
  }
  
  return(binary_results)
}

# 分别执行各组的视力参数二分类分析
vision_binary_analysis_group1 <- analyze_vision_binary(
  vision_with_clusters_group1,
  vision_params,
  "vision_by_steps_cluster",
  output_dir_group1,
  "Group1"
)

vision_binary_analysis_group2 <- analyze_vision_binary(
  vision_with_clusters_group2,
  vision_params,
  "vision_by_steps_cluster",
  output_dir_group2,
  "Group2"
)

# -----------------------------------------------------
# 6. 添加多时间点比较可视化
# -----------------------------------------------------
# 创建多时间点可视化函数
create_combined_timepoint_plot <- function(data, group_name, output_dir) {
  # 确保聚类是因子类型
  data$steps_cluster <- as.factor(data$steps_cluster)
  
  # 转换为长格式以便于绘图
  long_data <- data %>%
    pivot_longer(
      cols = c(vision_improvement_1w, vision_improvement_1m),
      names_to = "timepoint",
      values_to = "improvement"
    ) %>%
    mutate(
      timepoint = factor(
        timepoint,
        levels = c("vision_improvement_1w", "vision_improvement_1m"),
        labels = c("1周", "1个月")
      )
    )
  
  # 创建分面图
  p <- ggplot(long_data, aes(x = steps_cluster, y = improvement, fill = timepoint)) +
    geom_boxplot(width = 0.7, position = position_dodge(0.8)) +
    geom_point(aes(color = timepoint), alpha = 0.4, size = 3,
               position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
    stat_summary(aes(group = timepoint), fun = mean, geom = "point", 
                 shape = 23, size = 4, fill = "white",
                 position = position_dodge(0.8)) +
    labs(
      title = paste(group_name, "- 不同步数模式聚类组的视力改善趋势"),
      x = "步数聚类组",
      y = "视力改善",
      fill = "时间点",
      color = "时间点"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1")
  
  # 保存图表
  ggsave(
    file.path(output_dir, "vision_combined_timepoints.pdf"),
    p,
    width = 10,
    height = 8
  )
  
  return(p)
}

# 为每组创建多时间点比较图
combined_timepoint_plot_group1 <- create_combined_timepoint_plot(
  vision_with_clusters_group1, 
  "Group1", 
  output_dir_group1
)

combined_timepoint_plot_group2 <- create_combined_timepoint_plot(
  vision_with_clusters_group2, 
  "Group2", 
  output_dir_group2
)

# -----------------------------------------------------
# 7. 分析步数聚类成员度与视力改善的相关性
# -----------------------------------------------------
# 分析成员度相关性函数
analyze_membership_correlation <- function(data, vision_params, output_dir, group_name) {
  # 确保membership为数值型
  data$membership <- as.numeric(as.character(data$steps_cluster))
  
  results <- data.frame(
    Parameter = character(),
    Rho = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (param in vision_params) {
    # Spearman相关性测试
    corr_test <- cor.test(data$membership, data[[param]], 
                          method = "spearman", exact = FALSE)
    
    # 添加结果
    results <- results %>%
      rbind(data.frame(
        Parameter = gsub("vision_improvement_", "", param),
        Rho = corr_test$estimate,
        P_value = corr_test$p.value,
        stringsAsFactors = FALSE
      ))
  }
  
  # 按p值排序
  results <- results %>%
    arrange(P_value)
  
  # 保存结果
  write.csv(
    results,
    file.path(output_dir, "membership_correlation_vision.csv"),
    row.names = FALSE
  )
  
  return(results)
}

# 分析各组的成员度相关性
membership_corr_group1 <- analyze_membership_correlation(
  vision_with_clusters_group1,
  vision_params,
  output_dir_group1,
  "Group1"
)

membership_corr_group2 <- analyze_membership_correlation(
  vision_with_clusters_group2,
  vision_params,
  output_dir_group2,
  "Group2"
)

# -----------------------------------------------------
# 8. 生成各组分析报告
# -----------------------------------------------------
# 创建HTML报告函数
create_html_report <- function(vision_analysis, vision_binary_analysis, membership_corr, output_dir, group_name) {
  # 整理分析结果
  vision_summary <- vision_analysis$summary
  vision_sig <- vision_summary %>% filter(Significant)
  
  vision_binary_sig <- vision_binary_analysis %>% filter(Significant)
  
  sig_corr <- membership_corr %>% filter(P_value < 0.05)
  
  # 创建HTML报告
  sink(file.path(output_dir, "steps_cluster_vision_analysis_report.html"))
  
  cat('<!DOCTYPE html>
<html>
<head>
  <title>', group_name, ' - 步数模式聚类与视力改善关系分析报告</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #2980b9; margin-top: 30px; }
    h3 { color: #3498db; }
    h4 { color: #16a085; margin-top: 20px; }
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
  <h1>', group_name, ' - 步数模式聚类与视力改善关系分析报告</h1>
  
  <div class="summary">
    <p>本报告分析了', group_name, '患者的不同步数恢复模式聚类组之间视力改善值的差异。分析包括连续变量分析和基于改善与否的二分类变量分析。</p>
  </div>
  
  <h2>1. 连续变量分析结果</h2>')
  
  # 视力参数结果
  if(nrow(vision_sig) > 0) {
    cat('<p>显著的视力改善差异:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>P值</th><th>校正后P值</th><th>效应量</th><th>聚类1均值</th><th>聚类2均值</th><th>差异</th></tr>')
    
    for(i in 1:nrow(vision_sig)) {
      cat('<tr>')
      cat('<td>', vision_sig$Parameter[i], '</td>')
      cat('<td>', format(vision_sig$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(vision_sig$P_Adjusted[i], digits=3), '</td>')
      cat('<td>', format(vision_sig$Effect_Size[i], digits=3), '</td>')
      cat('<td>', format(vision_sig$Mean_Cluster1[i], digits=3), '</td>')
      cat('<td>', format(vision_sig$Mean_Cluster2[i], digits=3), '</td>')
      cat('<td>', format(vision_sig$Mean_Diff[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的视力改善差异。</p>')
  }
  
  # 二分类变量分析结果
  cat('<h2>2. 二分类变量分析结果</h2>')
  
  if(nrow(vision_binary_sig) > 0) {
    cat('<p>显著的视力改善OR值:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>分类方法</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(vision_binary_sig)) {
      cat('<tr>')
      cat('<td>', vision_binary_sig$Parameter[i], '</td>')
      cat('<td>', vision_binary_sig$Definition[i], '</td>')
      cat('<td>', format(vision_binary_sig$OR[i], digits=3), '</td>')
      cat('<td>', format(vision_binary_sig$OR_Lower[i], digits=3), '-', 
          format(vision_binary_sig$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(vision_binary_sig$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(vision_binary_sig$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的视力改善OR值差异。</p>')
  }
  
  # 成员度相关性结果
  cat('<h2>3. 聚类成员度与视力改善的相关性</h2>')
  cat('<table>')
  cat('<tr><th>参数</th><th>Spearman&#39;s Rho</th><th>P值</th></tr>')
  
  for(i in 1:nrow(membership_corr)) {
    cat('<tr>')
    cat('<td>', membership_corr$Parameter[i], '</td>')
    cat('<td>', format(membership_corr$Rho[i], digits=3), '</td>')
    if(membership_corr$P_value[i] < 0.05) {
      cat('<td class="significant">', format(membership_corr$P_value[i], digits=3), '</td>')
    } else {
      cat('<td>', format(membership_corr$P_value[i], digits=3), '</td>')
    }
    cat('</tr>')
  }
  
  cat('</table>')
  
  # 结论
  cat('<h2>4. 总结与结论</h2>')
  
  total_sig <- nrow(vision_sig) + nrow(vision_binary_sig) + nrow(sig_corr)
  
  cat('<p>本分析探讨了', group_name, '患者中不同步数恢复模式与视力改善之间的关系。总结如下：</p>')
  
  if(total_sig > 0) {
    cat('<p>步数聚类组间存在显著的视力改善差异，包括：</p>')
    cat('<ul>')
    if(nrow(vision_sig) > 0) {
      cat('<li>连续变量分析：', nrow(vision_sig), '个显著结果</li>')
      for(i in 1:nrow(vision_sig)) {
        direction <- ifelse(vision_sig$Mean_Diff[i] > 0, "高于", "低于")
        cat('<li>在', vision_sig$Parameter[i], '时间点，聚类2组视力改善', direction, '聚类1组 (差异=', 
            format(vision_sig$Mean_Diff[i], digits=3), ', p=', format(vision_sig$P_Adjusted[i], digits=3), ')</li>')
      }
    }
    if(nrow(vision_binary_sig) > 0) {
      cat('<li>二分类分析：', nrow(vision_binary_sig), '个显著结果</li>')
      for(i in 1:min(3, nrow(vision_binary_sig))) {
        direction <- ifelse(vision_binary_sig$OR[i] > 1, "更可能", "不太可能")
        cat('<li>', vision_binary_sig$Parameter[i], ' (', vision_binary_sig$Definition[i], 
            '): 聚类2组', direction, '有视力改善 (OR=', 
            format(vision_binary_sig$OR[i], digits=2), ', p=', format(vision_binary_sig$P_Adjusted[i], digits=3), ')</li>')
      }
      if(nrow(vision_binary_sig) > 3) {
        cat('<li>另外还有', nrow(vision_binary_sig) - 3, '个显著的二分类结果</li>')
      }
    }
    if(nrow(sig_corr) > 0) {
      cat('<li>成员度相关性：', nrow(sig_corr), '个显著相关性</li>')
      for(i in 1:nrow(sig_corr)) {
        direction <- ifelse(sig_corr$Rho[i] > 0, "正相关", "负相关")
        cat('<li>步数聚类成员度与', sig_corr$Parameter[i], '视力改善呈', direction, 
            ' (Rho=', format(sig_corr$Rho[i], digits=2), 
            ', p=', format(sig_corr$P_value[i], digits=3), ')</li>')
      }
    }
    cat('</ul>')
    
    cat('<p>这些结果表明', group_name, '患者的步数恢复模式与视力改善存在显著关联。建议进一步研究这些关联的临床意义。</p>')
  } else {
    cat('<p>在', group_name, '患者中，未发现步数聚类组间存在显著的视力改善差异。这可能表明：</p>')
    cat('<ul>')
    cat('<li>步数恢复模式与视力改善可能没有直接关系</li>')
    cat('<li>可能需要更多样本量或更长的随访时间来检测可能的差异</li>')
    cat('<li>可能需要考虑其他潜在的混杂因素</li>')
    cat('</ul>')
  }
  
  cat('<h2>5. 建议</h2>')
  cat('<p>基于本分析结果，我们提出以下建议：</p>')
  cat('<ol>')
  if(total_sig > 0) {
    cat('<li>临床实践中可考虑将步数恢复模式作为评估', group_name, '患者视力改善潜力的参考指标。</li>')
    cat('<li>建议进一步研究步数模式与视力改善的生理机制。</li>')
    cat('<li>考虑将步数模式整合到预测模型中，以更准确预测患者术后视力改善情况。</li>')
  } else {
    cat('<li>对于', group_name, '患者，步数模式可能不是预测视力改善的有效指标。</li>')
    cat('<li>建议探索其他可能影响视力改善的因素。</li>')
  }
  cat('<li>设计更大规模的前瞻性研究，进一步验证步数模式与视力改善的关联性。</li>')
  cat('</ol>')
  
  cat('</body>
</html>')
  
  sink()
  
  cat("HTML报告已生成：", file.path(output_dir, "steps_cluster_vision_analysis_report.html"), "\n")
}

# 为每组生成HTML报告
create_html_report(
  vision_analysis_group1, 
  vision_binary_analysis_group1, 
  membership_corr_group1, 
  output_dir_group1, 
  "Group1"
)

create_html_report(
  vision_analysis_group2, 
  vision_binary_analysis_group2, 
  membership_corr_group2, 
  output_dir_group2, 
  "Group2"
)

# 打印分析摘要
cat("\n=====================================================\n")
cat("步数模式聚类与视力改善关系分析完成\n")
cat("=====================================================\n\n")

cat("Group1分析结果已保存到目录：", output_dir_group1, "\n")
cat("Group2分析结果已保存到目录：", output_dir_group2, "\n\n")

cat("Group1视力连续变量分析显著结果数量：", sum(vision_analysis_group1$summary$Significant), "\n")
cat("Group2视力连续变量分析显著结果数量：", sum(vision_analysis_group2$summary$Significant), "\n\n")

cat("Group1视力二分类分析显著结果数量：", sum(vision_binary_analysis_group1$Significant), "\n")
cat("Group2视力二分类分析显著结果数量：", sum(vision_binary_analysis_group2$Significant), "\n\n")

cat("Group1视力改善与聚类成员度相关性显著结果数量：", sum(membership_corr_group1$P_value < 0.05), "\n")
cat("Group2视力改善与聚类成员度相关性显著结果数量：", sum(membership_corr_group2$P_value < 0.05), "\n\n")

# 保存处理后的数据集，以便进一步分析
saveRDS(vision_with_clusters_group1, file.path(output_dir_group1, "vision_with_steps_clusters.rds"))
saveRDS(vision_with_clusters_group2, file.path(output_dir_group2, "vision_with_steps_clusters.rds"))

