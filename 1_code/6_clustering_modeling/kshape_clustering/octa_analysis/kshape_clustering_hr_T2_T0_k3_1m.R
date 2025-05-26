# -----------------------------------------------------
# 分析不同心率模式聚类(K=3)对OCTA改善的影响
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
library(emmeans)      # 用于多组比较
library(multcomp)     # 用于多重比较

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# -----------------------------------------------------
# 1. 加载必要的数据 - 修改为使用k=3的聚类
# -----------------------------------------------------
# 加载糖尿病组(Group2)的心率聚类结果
group2_clusters <- readRDS("3_data_analysis/6_clustering_modeling/kshape/time_series/1m/Group2_Diabetes_all_clusters.rds")

# 提取k=3的聚类结果
cluster_k3 <- group2_clusters$k3$clusters
cluster_assignments <- data.frame(
  subject_id = names(cluster_k3),
  hr_cluster = as.factor(cluster_k3)
)

# 重新编码集群，确保1是参考组
# 这一步很重要，因为聚类算法可能不会按照1,2,3的顺序标记集群
# cluster_assignments$hr_cluster <- factor(cluster_assignments$hr_cluster, levels = c("1", "2", "3"))

# 打印集群分布，查看每个集群的样本数
print(table(cluster_assignments$hr_cluster))

# 加载OCTA数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

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
# 5. 合并心率聚类结果和OCTA改善值
# -----------------------------------------------------
# 合并血流数据
bloodflow_with_clusters <- bloodflow_filtered$data %>%
  inner_join(cluster_assignments, by = c("ID" = "subject_id"))

# 合并厚度数据
thickness_with_clusters <- thickness_filtered$data %>%
  inner_join(cluster_assignments, by = c("ID" = "subject_id"))

# 创建结果目录
output_dir <- "3_data_analysis/6_clustering_modeling/octa_hr_cluster_analysis_k3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# 6. 连续变量分析：修改为K=3聚类组的多组比较
# -----------------------------------------------------
# 修改统计分析函数以处理K=3聚类组
analyze_octa_by_cluster_k3 <- function(data, param_cols, output_prefix) {
  results_list <- list()
  all_p_values <- c()
  param_names <- c()
  
  # 确保聚类变量为因子类型
  data$hr_cluster <- as.factor(data$hr_cluster)
  
  # 创建汇总表
  summary_df <- data.frame(
    Parameter = character(),
    ANOVA_P_Value = numeric(),
    ANOVA_P_Adjusted = numeric(),
    Significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # 创建用于存储两两比较的数据框
  pairwise_df <- data.frame(
    Parameter = character(),
    Comparison = character(),
    Mean_Diff = numeric(),
    P_Value = numeric(),
    P_Adjusted = numeric(),
    Significant = logical(),
    Effect_Size = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 为每个参数创建可视化和执行统计测试
  for (param in param_cols) {
    param_name <- gsub("_0_21_improvement", "", param)
    
    # 创建公式
    formula_str <- paste(param, "~ hr_cluster")
    
    # 执行ANOVA
    anova_result <- data %>%
      anova_test(as.formula(formula_str))
    
    # 保存p值
    p_value <- anova_result$p
    all_p_values <- c(all_p_values, p_value)
    param_names <- c(param_names, param_name)
    
    # 计算每个聚类的均值
    cluster_means <- data %>%
      group_by(hr_cluster) %>%
      summarize(
        mean_value = mean(!!sym(param), na.rm = TRUE),
        sd_value = sd(!!sym(param), na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    # 添加到汇总表
    summary_df <- rbind(summary_df, data.frame(
      Parameter = param_name,
      ANOVA_P_Value = p_value,
      ANOVA_P_Adjusted = NA,  # 将在之后填充
      Significant = NA, # 将在之后填充
      stringsAsFactors = FALSE
    ))
    
    # 执行事后检验 (Tukey HSD)
    posthoc <- data %>%
      tukey_hsd(as.formula(formula_str))
    
    # 计算对比组1的效应量
    # 所有涉及组1的比较
    group1_comparisons <- posthoc %>%
      filter(group1 == "1" | group2 == "1")
    
    # 收集组2 vs 组1和组3 vs 组1的比较结果
    if (nrow(group1_comparisons) > 0) {
      for (i in 1:nrow(group1_comparisons)) {
        comp_row <- group1_comparisons[i, ]
        
        # 确定比较组 - 添加NA检查
        comparison <- paste(comp_row$group2, "vs", comp_row$group1)
        if (!is.na(comp_row$group1) && comp_row$group1 == "1") {
          comparison <- paste(comp_row$group2, "vs 1")
        } else if (!is.na(comp_row$group2) && comp_row$group2 == "1") {
          comparison <- paste(comp_row$group1, "vs 1")
        }
        
        # 计算效应量 (Cohen's d)
        # 为简单起见，我们可以直接使用均值差除以合并标准差作为效应量估计
        # 获取相关组的均值和标准差
        if (!is.na(comp_row$group1) && comp_row$group1 == "1") {
          group1_stats <- cluster_means[cluster_means$hr_cluster == "1", ]
          group2_stats <- cluster_means[cluster_means$hr_cluster == comp_row$group2, ]
        } else if (!is.na(comp_row$group2) && comp_row$group2 == "1") {
          group1_stats <- cluster_means[cluster_means$hr_cluster == "1", ]
          group2_stats <- cluster_means[cluster_means$hr_cluster == comp_row$group1, ]
        } else {
          # 处理无法确定比较的情况 - 跳过计算效应量
          next
        }
        
        # 计算合并标准差
        pooled_sd <- sqrt(((group1_stats$n - 1) * group1_stats$sd_value^2 + 
                             (group2_stats$n - 1) * group2_stats$sd_value^2) / 
                            (group1_stats$n + group2_stats$n - 2))
        
        # 计算效应量
        effect_size <- abs(group1_stats$mean_value - group2_stats$mean_value) / pooled_sd
        
        # 添加到两两比较数据框
        pairwise_df <- rbind(pairwise_df, data.frame(
          Parameter = param_name,
          Comparison = comparison,
          Mean_Diff = comp_row$estimate,
          P_Value = comp_row$p.adj,
          P_Adjusted = NA,  # 将在之后填充
          Significant = NA, # 将在之后填充
          Effect_Size = effect_size,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # 创建可视化 - 使用ggbetweenstats
    p <- ggstatsplot::ggbetweenstats(
      data = data,
      x = hr_cluster,
      y = !!sym(param),
      type = "parametric",
      pairwise.display = "all",
      p.adjust.method = "fdr",
      effsize.type = "unbiased",
      results.subtitle = TRUE,
      xlab = "心率聚类组",
      ylab = paste(param_name, "改善值 (T2-T0)"),
      title = paste0(param_name, " 在不同心率模式组间的差异"),
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
      anova = anova_result,
      posthoc = posthoc,
      cluster_means = cluster_means
    )
  }
  
  # 多重比较校正 - ANOVA结果
  summary_df$ANOVA_P_Adjusted <- p.adjust(summary_df$ANOVA_P_Value, method = "fdr")
  summary_df$Significant <- summary_df$ANOVA_P_Adjusted < 0.05
  
  # 多重比较校正 - 两两比较结果
  if (nrow(pairwise_df) > 0) {
    pairwise_df$P_Adjusted <- p.adjust(pairwise_df$P_Value, method = "fdr")
    pairwise_df$Significant <- pairwise_df$P_Adjusted < 0.05
  }
  
  # 按校正后的p值排序
  summary_df <- summary_df %>%
    arrange(ANOVA_P_Adjusted)
  
  pairwise_df <- pairwise_df %>%
    arrange(P_Adjusted)
  
  # 保存汇总表
  write.csv(
    summary_df,
    file.path(output_dir, paste0(output_prefix, "_anova_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    pairwise_df,
    file.path(output_dir, paste0(output_prefix, "_pairwise_summary.csv")),
    row.names = FALSE
  )
  
  # 创建top显著参数的多面板图
  if (any(summary_df$Significant)) {
    # 选择显著的参数
    sig_params <- summary_df %>%
      filter(Significant) %>%
      arrange(ANOVA_P_Adjusted) %>%
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
        x = hr_cluster,
        y = !!sym(param),
        type = "parametric",
        pairwise.display = "all",
        p.adjust.method = "fdr",
        effsize.type = "unbiased",
        results.subtitle = TRUE,
        xlab = "心率聚类组",
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
    summary = summary_df,
    pairwise = pairwise_df
  ))
}

# 执行血流参数连续变量分析
bloodflow_analysis <- analyze_octa_by_cluster_k3(
  bloodflow_with_clusters,
  bloodflow_filtered$params,
  "diabetes_bloodflow_by_hr_cluster_k3"
)

# 执行厚度参数连续变量分析
thickness_analysis <- analyze_octa_by_cluster_k3(
  thickness_with_clusters,
  thickness_filtered$params,
  "diabetes_thickness_by_hr_cluster_k3"
)

# -----------------------------------------------------
# 7. 二分类变量分析：修改为计算组2 vs 组1和组3 vs 组1的OR
# -----------------------------------------------------
# 修改二分类分析函数以处理K=3聚类组
# First, make sure to add this package to your library section at the top
library(nnet)  # For multinomial logistic regression

# Modified function to use multinomial logistic regression

analyze_octa_binary_multinomial_k3 <- function(data, param_cols, output_prefix) {
  # 确保聚类变量为因子类型
  data$hr_cluster <- as.factor(data$hr_cluster)
  
  # 初始化结果
  fisher_results <- data.frame(
    Parameter = character(),
    Comparison = character(),
    P_Value = numeric(),
    P_Adjusted = numeric(),
    Significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # 循环计算每个参数的Fisher精确检验结果
  for (param in param_cols) {
    param_name <- gsub("_0_21_improvement", "", param)
    
    # 创建二分类变量（1表示改善>0，0表示改善≤0）
    binary_var_name <- paste0(param, "_improved")
    data[[binary_var_name]] <- ifelse(data[[param]] > 0, 1, 0)
    
    # 计算改善比例
    improve_ratio <- mean(data[[binary_var_name]], na.rm = TRUE)
    cat(param_name, "改善比例:", round(improve_ratio * 100, 1), "%\n")
    
    # 使用更简单的方法 - 对每个聚类组(2和3)单独与参考组(1)进行比较
    ref_group <- "1"  # 参考组
    
    for (compare_group in setdiff(levels(data$hr_cluster), ref_group)) {
      # 创建仅包含参考组和比较组的子集
      subset_data <- data %>% 
        filter(hr_cluster %in% c(ref_group, compare_group))
      
      # 创建列联表
      contingency_table <- table(subset_data$hr_cluster, subset_data[[binary_var_name]])
      
      # 执行Fisher精确检验
      tryCatch({
        fisher_test <- fisher.test(contingency_table)
        p_value <- fisher_test$p.value
        
        # 添加到结果
        fisher_results <- rbind(fisher_results, data.frame(
          Parameter = param_name,
          Comparison = paste(compare_group, "vs", ref_group),
          P_Value = p_value,
          P_Adjusted = NA,  # 稍后填充
          Significant = NA, # 稍后填充
          stringsAsFactors = FALSE
        ))
        
        # 创建每个组别的改善/不改善分布可视化
        plot_data <- subset_data %>%
          group_by(hr_cluster) %>%
          summarise(
            improved = mean(!!sym(binary_var_name) == 1) * 100,
            not_improved = mean(!!sym(binary_var_name) == 0) * 100,
            .groups = "drop"
          ) %>%
          pivot_longer(
            cols = c(improved, not_improved),
            names_to = "improvement_status",
            values_to = "percentage"
          ) %>%
          mutate(
            improvement_status = factor(improvement_status, 
                                        levels = c("not_improved", "improved"),
                                        labels = c("无改善/恶化", "有改善"))
          )
        
        # 创建堆叠条形图
        p <- ggplot(plot_data, aes(x = hr_cluster, y = percentage, fill = improvement_status)) +
          geom_bar(stat = "identity", position = "stack") +
          geom_text(aes(label = sprintf("%.1f%%", percentage)),
                    position = position_stack(vjust = 0.5)) +
          labs(
            title = paste0(param_name, " - ", compare_group, " vs ", ref_group, " 改善比例"),
            subtitle = paste0("分类标准: >0为改善\n",
                              "Fisher精确检验 p-value: ", format.pval(p_value, digits = 3)),
            x = "心率聚类组",
            y = "比例 (%)",
            fill = "改善状态"
          ) +
          theme_bw() +
          scale_fill_brewer(palette = "Set2")
        
        # 保存图表
        ggsave(
          file.path(output_dir, paste0(output_prefix, "_fisher_", param_name, "_", 
                                       compare_group, "vs", ref_group, ".pdf")),
          p,
          width = 8,
          height = 6
        )
        
      }, error = function(e) {
        cat("执行Fisher精确检验时出错:", param_name, "\n")
        cat("错误信息:", e$message, "\n")
      })
    }
  }
  
  # 多重比较校正
  if (nrow(fisher_results) > 0) {
    params <- unique(fisher_results$Parameter)
    
    for (param in params) {
      param_rows <- fisher_results$Parameter == param
      fisher_results$P_Adjusted[param_rows] <- p.adjust(fisher_results$P_Value[param_rows], method = "fdr")
    }
    
    fisher_results$Significant <- fisher_results$P_Adjusted < 0.05
    
    # 按校正后的p值排序
    fisher_results <- fisher_results %>%
      arrange(Parameter, P_Adjusted)
    
    # 保存汇总表
    write.csv(
      fisher_results,
      file.path(output_dir, paste0(output_prefix, "_fisher_results.csv")),
      row.names = FALSE
    )
  } else {
    cat("没有成功执行任何Fisher精确检验\n")
  }
  
  return(fisher_results)
}

# 执行血流参数多项式回归分析
bloodflow_multinomial_analysis <- analyze_octa_binary_multinomial_k3(
  bloodflow_with_clusters,
  bloodflow_filtered$params,
  "diabetes_bloodflow_by_hr_cluster_k3"
)

# 执行厚度参数多项式回归分析
thickness_multinomial_analysis <- analyze_octa_binary_multinomial_k3(
  thickness_with_clusters,
  thickness_filtered$params,
  "diabetes_thickness_by_hr_cluster_k3"
)

# -----------------------------------------------------
# 8. 生成综合报告
# -----------------------------------------------------
# 创建HTML报告函数
create_html_report <- function() {
  # 血流连续变量分析摘要
  bloodflow_summary <- bloodflow_analysis$summary
  bloodflow_sig <- bloodflow_summary %>% filter(Significant)
  
  # 血流两两比较结果
  bloodflow_pairwise <- bloodflow_analysis$pairwise
  bloodflow_pairwise_sig <- bloodflow_pairwise %>% filter(Significant)
  
  # 厚度连续变量分析摘要
  thickness_summary <- thickness_analysis$summary
  thickness_sig <- thickness_summary %>% filter(Significant)
  
  # 厚度两两比较结果
  thickness_pairwise <- thickness_analysis$pairwise
  thickness_pairwise_sig <- thickness_pairwise %>% filter(Significant)
  
  # 血流二分类分析摘要
  bloodflow_binary_2vs1 <- bloodflow_multinomial_analysis %>% 
    filter(Comparison == "2 vs 1" & Significant)
  bloodflow_binary_3vs1 <- bloodflow_multinomial_analysis %>% 
    filter(Comparison == "3 vs 1" & Significant)
  
  # 厚度二分类分析摘要
  thickness_binary_2vs1 <- thickness_multinomial_analysis %>% 
    filter(Comparison == "2 vs 1" & Significant)
  thickness_binary_3vs1 <- thickness_multinomial_analysis %>% 
    filter(Comparison == "3 vs 1" & Significant)
  
  # 创建HTML报告
  sink(file.path(output_dir, "heart_rate_cluster_k3_octa_analysis_report.html"))
  
  cat('<!DOCTYPE html>
<html>
<head>
  <title>心率模式聚类(K=3)与OCTA改善关系分析报告</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
    h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    h2 { color: #2980b9; margin-top: 30px; }
    h3 { color: #3498db; }
    h4 { color: #2c3e50; margin-top: 20px; }
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
  <h1>心率模式聚类(K=3)与OCTA改善关系分析报告</h1>
  
  <div class="summary">
    <p>本报告分析了糖尿病患者(Group2)的不同心率恢复模式聚类组(K=3)之间OCTA参数改善值的差异。分析包括连续变量分析和基于中位数的二分类变量分析，以组1为参照组，分析组2和组3对比组1的差异。</p>
  </div>
  
  <h2>1. 连续变量分析结果</h2>')
  
  # 血流参数ANOVA结果
  cat('<h3>1.1 血流参数ANOVA结果</h3>')
  
  if(nrow(bloodflow_sig) > 0) {
    cat('<p>显著的血流参数ANOVA结果:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>ANOVA P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(bloodflow_sig)) {
      cat('<tr>')
      cat('<td>', bloodflow_sig$Parameter[i], '</td>')
      cat('<td>', format(bloodflow_sig$ANOVA_P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(bloodflow_sig$ANOVA_P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
    
    # 显示有显著ANOVA结果的参数的两两比较
    if(nrow(bloodflow_pairwise_sig) > 0) {
      cat('<h4>显著的血流参数两两比较结果:</h4>')
      cat('<table>')
      cat('<tr><th>参数</th><th>比较组</th><th>均值差</th><th>P值</th><th>校正后P值</th><th>效应量</th></tr>')
      
      for(i in 1:nrow(bloodflow_pairwise_sig)) {
        cat('<tr>')
        cat('<td>', bloodflow_pairwise_sig$Parameter[i], '</td>')
        cat('<td>', bloodflow_pairwise_sig$Comparison[i], '</td>')
        cat('<td>', format(bloodflow_pairwise_sig$Mean_Diff[i], digits=3), '</td>')
        cat('<td>', format(bloodflow_pairwise_sig$P_Value[i], digits=3), '</td>')
        cat('<td class="significant">', format(bloodflow_pairwise_sig$P_Adjusted[i], digits=3), '</td>')
        cat('<td>', format(bloodflow_pairwise_sig$Effect_Size[i], digits=3), '</td>')
        cat('</tr>')
      }
      
      cat('</table>')
    } else {
      cat('<p class="non-significant">尽管ANOVA显示有总体差异，但多重比较校正后没有发现显著的两两比较结果。</p>')
    }
  } else {
    cat('<p class="non-significant">没有发现显著的血流参数ANOVA结果。</p>')
  }
  
  # 厚度参数ANOVA结果
  cat('<h3>1.2 厚度参数ANOVA结果</h3>')
  
  if(nrow(thickness_sig) > 0) {
    cat('<p>显著的厚度参数ANOVA结果:</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>ANOVA P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(thickness_sig)) {
      cat('<tr>')
      cat('<td>', thickness_sig$Parameter[i], '</td>')
      cat('<td>', format(thickness_sig$ANOVA_P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(thickness_sig$ANOVA_P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
    
    # 显示有显著ANOVA结果的参数的两两比较
    if(nrow(thickness_pairwise_sig) > 0) {
      cat('<h4>显著的厚度参数两两比较结果:</h4>')
      cat('<table>')
      cat('<tr><th>参数</th><th>比较组</th><th>均值差</th><th>P值</th><th>校正后P值</th><th>效应量</th></tr>')
      
      for(i in 1:nrow(thickness_pairwise_sig)) {
        cat('<tr>')
        cat('<td>', thickness_pairwise_sig$Parameter[i], '</td>')
        cat('<td>', thickness_pairwise_sig$Comparison[i], '</td>')
        cat('<td>', format(thickness_pairwise_sig$Mean_Diff[i], digits=3), '</td>')
        cat('<td>', format(thickness_pairwise_sig$P_Value[i], digits=3), '</td>')
        cat('<td class="significant">', format(thickness_pairwise_sig$P_Adjusted[i], digits=3), '</td>')
        cat('<td>', format(thickness_pairwise_sig$Effect_Size[i], digits=3), '</td>')
        cat('</tr>')
      }
      
      cat('</table>')
    } else {
      cat('<p class="non-significant">尽管ANOVA显示有总体差异，但多重比较校正后没有发现显著的两两比较结果。</p>')
    }
  } else {
    cat('<p class="non-significant">没有发现显著的厚度参数ANOVA结果。</p>')
  }
  
  # 二分类变量分析结果
  cat('<h2>2. 二分类变量分析结果 (基于中位数)</h2>')
  
  # 血流参数二分类结果 - 组2 vs 组1
  cat('<h3>2.1 血流参数 - 组2 vs 组1</h3>')
  
  if(nrow(bloodflow_binary_2vs1) > 0) {
    cat('<p>显著的血流参数OR值 (组2 vs 组1):</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>中位数切点</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(bloodflow_binary_2vs1)) {
      cat('<tr>')
      cat('<td>', bloodflow_binary_2vs1$Parameter[i], '</td>')
      cat('<td>', format(bloodflow_binary_2vs1$Median_Cutoff[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_2vs1$OR[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_2vs1$OR_Lower[i], digits=3), '-', 
          format(bloodflow_binary_2vs1$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_2vs1$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(bloodflow_binary_2vs1$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的血流参数OR值差异 (组2 vs 组1)。</p>')
  }
  
  # 血流参数二分类结果 - 组3 vs 组1
  cat('<h3>2.2 血流参数 - 组3 vs 组1</h3>')
  
  if(nrow(bloodflow_binary_3vs1) > 0) {
    cat('<p>显著的血流参数OR值 (组3 vs 组1):</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>中位数切点</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(bloodflow_binary_3vs1)) {
      cat('<tr>')
      cat('<td>', bloodflow_binary_3vs1$Parameter[i], '</td>')
      cat('<td>', format(bloodflow_binary_3vs1$Median_Cutoff[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_3vs1$OR[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_3vs1$OR_Lower[i], digits=3), '-', 
          format(bloodflow_binary_3vs1$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(bloodflow_binary_3vs1$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(bloodflow_binary_3vs1$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的血流参数OR值差异 (组3 vs 组1)。</p>')
  }
  
  # 厚度参数二分类结果 - 组2 vs 组1
  cat('<h3>2.3 厚度参数 - 组2 vs 组1</h3>')
  
  if(nrow(thickness_binary_2vs1) > 0) {
    cat('<p>显著的厚度参数OR值 (组2 vs 组1):</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>中位数切点</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(thickness_binary_2vs1)) {
      cat('<tr>')
      cat('<td>', thickness_binary_2vs1$Parameter[i], '</td>')
      cat('<td>', format(thickness_binary_2vs1$Median_Cutoff[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_2vs1$OR[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_2vs1$OR_Lower[i], digits=3), '-', 
          format(thickness_binary_2vs1$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_2vs1$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(thickness_binary_2vs1$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的厚度参数OR值差异 (组2 vs 组1)。</p>')
  }
  
  # 厚度参数二分类结果 - 组3 vs 组1
  cat('<h3>2.4 厚度参数 - 组3 vs 组1</h3>')
  
  if(nrow(thickness_binary_3vs1) > 0) {
    cat('<p>显著的厚度参数OR值 (组3 vs 组1):</p>')
    cat('<table>')
    cat('<tr><th>参数</th><th>中位数切点</th><th>OR值</th><th>95% CI</th><th>P值</th><th>校正后P值</th></tr>')
    
    for(i in 1:nrow(thickness_binary_3vs1)) {
      cat('<tr>')
      cat('<td>', thickness_binary_3vs1$Parameter[i], '</td>')
      cat('<td>', format(thickness_binary_3vs1$Median_Cutoff[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_3vs1$OR[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_3vs1$OR_Lower[i], digits=3), '-', 
          format(thickness_binary_3vs1$OR_Upper[i], digits=3), '</td>')
      cat('<td>', format(thickness_binary_3vs1$P_Value[i], digits=3), '</td>')
      cat('<td class="significant">', format(thickness_binary_3vs1$P_Adjusted[i], digits=3), '</td>')
      cat('</tr>')
    }
    
    cat('</table>')
  } else {
    cat('<p class="non-significant">没有发现显著的厚度参数OR值差异 (组3 vs 组1)。</p>')
  }
  
  # 结论
  cat('<h2>3. 结论</h2>')
  
  total_sig <- nrow(bloodflow_sig) + nrow(thickness_sig) + 
    nrow(bloodflow_binary_2vs1) + nrow(bloodflow_binary_3vs1) +
    nrow(thickness_binary_2vs1) + nrow(thickness_binary_3vs1)
  
  if(total_sig > 0) {
    cat('<p>本分析发现了心率模式聚类组之间存在显著的OCTA参数改善差异，具体如下：</p>')
    cat('<ul>')
    
    if(nrow(bloodflow_sig) > 0) {
      cat('<li>血流参数ANOVA分析：发现', nrow(bloodflow_sig), '个显著差异参数</li>')
      if(nrow(bloodflow_pairwise_sig) > 0) {
        cat('<li>血流参数两两比较：发现', nrow(bloodflow_pairwise_sig), '个显著差异的两两比较</li>')
      }
    }
    
    if(nrow(thickness_sig) > 0) {
      cat('<li>厚度参数ANOVA分析：发现', nrow(thickness_sig), '个显著差异参数</li>')
      if(nrow(thickness_pairwise_sig) > 0) {
        cat('<li>厚度参数两两比较：发现', nrow(thickness_pairwise_sig), '个显著差异的两两比较</li>')
      }
    }
    
    if(nrow(bloodflow_binary_2vs1) > 0) {
      cat('<li>血流参数二分类分析 (组2 vs 组1)：发现', nrow(bloodflow_binary_2vs1), '个显著OR值</li>')
    }
    
    if(nrow(bloodflow_binary_3vs1) > 0) {
      cat('<li>血流参数二分类分析 (组3 vs 组1)：发现', nrow(bloodflow_binary_3vs1), '个显著OR值</li>')
    }
    
    if(nrow(thickness_binary_2vs1) > 0) {
      cat('<li>厚度参数二分类分析 (组2 vs 组1)：发现', nrow(thickness_binary_2vs1), '个显著OR值</li>')
    }
    
    if(nrow(thickness_binary_3vs1) > 0) {
      cat('<li>厚度参数二分类分析 (组3 vs 组1)：发现', nrow(thickness_binary_3vs1), '个显著OR值</li>')
    }
    
    cat('</ul>')
    
    cat('<p>这些结果表明不同的心率恢复模式可能与术后OCTA参数的改善情况相关。建议进一步研究这些关联的临床意义，并探索不同心率模式组的临床特征差异。</p>')
  } else {
    cat('<p>本分析未发现心率模式聚类组之间存在显著的OCTA参数改善差异。这可能表明：</p>')
    cat('<ul>')
    cat('<li>心率恢复模式与OCTA参数改善可能没有直接关系</li>')
    cat('<li>可能需要更多样本量或更长的随访时间来检测可能的差异</li>')
    cat('<li>可能需要考虑其他潜在的混杂因素</li>')
    cat('<li>聚类算法可能需要进一步优化或尝试其他聚类方法</li>')
    cat('</ul>')
  }
  
  cat('<p>建议结合患者的临床特征和其他检查结果进行综合分析，并考虑探索每个聚类组的特征。</p>')
  
  cat('</body>
</html>')
  
  sink()
  
  cat("HTML报告已生成：", file.path(output_dir, "heart_rate_cluster_k3_octa_analysis_report.html"), "\n")
}

# 生成HTML报告
create_html_report()

# 打印分析摘要
cat("\n=====================================================\n")
cat("心率模式聚类(K=3)与OCTA改善关系分析完成\n")
cat("=====================================================\n\n")

cat("分析结果已保存到目录：", output_dir, "\n\n")

cat("血流参数ANOVA显著结果数量：", sum(bloodflow_analysis$summary$Significant), "\n")
cat("厚度参数ANOVA显著结果数量：", sum(thickness_analysis$summary$Significant), "\n")

cat("血流参数两两比较显著结果数量：", sum(bloodflow_analysis$pairwise$Significant), "\n")
cat("厚度参数两两比较显著结果数量：", sum(thickness_analysis$pairwise$Significant), "\n")

cat("血流参数二分类分析 (组2 vs 组1) 显著结果数量：", 
    sum(bloodflow_multinomial_analysis$Significant & bloodflow_multinomial_analysis$Comparison == "2 vs 1"), "\n")
cat("血流参数二分类分析 (组3 vs 组1) 显著结果数量：", 
    sum(bloodflow_multinomial_analysis$Significant & bloodflow_multinomial_analysis$Comparison == "3 vs 1"), "\n")

cat("厚度参数二分类分析 (组2 vs 组1) 显著结果数量：", 
    sum(thickness_multinomial_analysis$Significant & thickness_multinomial_analysis$Comparison == "2 vs 1"), "\n")
cat("厚度参数二分类分析 (组3 vs 组1) 显著结果数量：", 
    sum(thickness_multinomial_analysis$Significant & thickness_multinomial_analysis$Comparison == "3 vs 1"), "\n\n")

cat("查看详细报告请打开生成的HTML文件。\n")
