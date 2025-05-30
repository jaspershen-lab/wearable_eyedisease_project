# 修正分析：使用聚类membership值而非原始指标
# Corrected Analysis: Using Cluster Membership Values

library(tidyverse)
library(ggplot2)

# ================== 1. 加载聚类membership数据 ==================
cat("===== 修正分析：使用聚类membership值 =====\n")

# 加载可穿戴设备聚类结果（包含membership值）
wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/ppv_cluster_results_time_windows.csv", check.names = FALSE)

# 检查数据结构
cat("可穿戴设备聚类数据结构:\n")
print(names(wearable_clusters))
print(head(wearable_clusters))

# 获取OCTA改善数据（重用之前的处理）
# 假设 octa_improvements 已经从之前的分析中获得
# 如果没有，需要重新运行OCTA数据处理部分

# ================== 2. 合并membership和OCTA数据 ==================
# 正确的分析：membership vs OCTA改善
membership_octa_analysis <- wearable_clusters %>%
  left_join(octa_improvements, by = c("subject_id" = "ID"))

cat("\n合并后的分析数据:\n")
cat("总患者数:", nrow(membership_octa_analysis), "\n")
cat("有效membership值:", sum(!is.na(membership_octa_analysis$max_membership)), "\n")

# 检查数据列名
cat("\n数据列名（前20个）:\n")
print(names(membership_octa_analysis)[1:min(20, length(names(membership_octa_analysis)))])

# ================== 3. Membership与OCTA改善的相关分析 ==================
perform_membership_correlation <- function(data, octa_params) {
  results <- data.frame()
  
  membership_var <- "max_membership"
  
  cat(sprintf("\n分析 %s 与OCTA改善参数的相关性...\n", membership_var))
  
  for(octa_param in octa_params) {
    if(!octa_param %in% names(data)) {
      cat(sprintf("跳过参数: %s (不存在)\n", octa_param))
      next
    }
    
    # 创建完整案例数据
    complete_data <- data[!is.na(data[[membership_var]]) & !is.na(data[[octa_param]]), ]
    
    if(nrow(complete_data) >= 3) {
      # Pearson相关
      pearson_test <- try(cor.test(complete_data[[membership_var]], complete_data[[octa_param]], 
                                   method = "pearson"), silent = TRUE)
      
      # Spearman相关
      spearman_test <- try(cor.test(complete_data[[membership_var]], complete_data[[octa_param]], 
                                    method = "spearman"), silent = TRUE)
      
      if(class(pearson_test) != "try-error" && class(spearman_test) != "try-error") {
        # 参数分类
        param_type <- case_when(
          grepl("SVP|ICP|DCP|Choroid", octa_param) ~ "BloodFlow",
          grepl("GCL|INL|Retina", octa_param) ~ "Thickness",
          TRUE ~ "Other"
        )
        
        region <- case_when(
          grepl("0_21", octa_param) ~ "Macular",
          grepl("0_6", octa_param) ~ "Widefield",
          TRUE ~ "Other"
        )
        
        # 效应量分类
        effect_size <- case_when(
          abs(pearson_test$estimate) >= 0.5 ~ "Large",
          abs(pearson_test$estimate) >= 0.3 ~ "Medium",
          abs(pearson_test$estimate) >= 0.1 ~ "Small",
          TRUE ~ "Negligible"
        )
        
        results <- rbind(results, data.frame(
          OCTA_Parameter = octa_param,
          Parameter_Type = param_type,
          Region = region,
          N = nrow(complete_data),
          Pearson_r = as.numeric(pearson_test$estimate),
          Pearson_p = pearson_test$p.value,
          Pearson_CI_lower = pearson_test$conf.int[1],
          Pearson_CI_upper = pearson_test$conf.int[2],
          Spearman_rho = as.numeric(spearman_test$estimate),
          Spearman_p = spearman_test$p.value,
          Effect_Size = effect_size,
          Significant_p05 = pearson_test$p.value < 0.05,
          Trend_p10 = pearson_test$p.value < 0.10,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      cat(sprintf("跳过参数: %s (样本量不足: %d)\n", octa_param, nrow(complete_data)))
    }
  }
  
  # FDR校正
  if(nrow(results) > 0) {
    results$Pearson_p_FDR <- p.adjust(results$Pearson_p, method = "fdr")
    results$Spearman_p_FDR <- p.adjust(results$Spearman_p, method = "fdr")
    results$Significant_FDR <- results$Pearson_p_FDR < 0.05
    
    # 按相关性强度排序
    results <- results %>% arrange(desc(abs(Pearson_r)))
  }
  
  return(results)
}

# 执行membership相关分析
membership_correlations <- perform_membership_correlation(
  membership_octa_analysis, 
  octa_improvement_params
)

# ================== 4. 显示membership分析结果 ==================
cat("\n===== 可穿戴设备Membership与OCTA改善相关分析 =====\n")

if(nrow(membership_correlations) > 0) {
  # 显著结果
  significant_results <- membership_correlations %>%
    filter(Significant_p05 == TRUE) %>%
    arrange(desc(abs(Pearson_r)))
  
  if(nrow(significant_results) > 0) {
    cat("🎯 显著相关结果 (p < 0.05):\n")
    print(significant_results %>%
            dplyr::select(OCTA_Parameter, Parameter_Type, Region, 
                          Pearson_r, Pearson_p, Effect_Size, N))
  } else {
    cat("❌ 未发现显著相关 (p < 0.05)\n")
  }
  
  # 趋势性结果
  trend_results <- membership_correlations %>%
    filter(Trend_p10 == TRUE & abs(Pearson_r) >= 0.4) %>%
    arrange(Pearson_p)
  
  if(nrow(trend_results) > 0) {
    cat("\n📈 趋势性显著结果 (p < 0.10, |r| ≥ 0.4):\n")
    print(trend_results %>%
            dplyr::select(OCTA_Parameter, Parameter_Type, Region, 
                          Pearson_r, Pearson_p, Effect_Size, N))
  }
  
  # 最强相关性（无论是否显著）
  cat("\n🔍 最强相关性前10个:\n")
  top_correlations <- membership_correlations %>%
    head(10) %>%
    mutate(
      Clinical_Relevance = case_when(
        Significant_p05 ~ "Statistically Significant",
        Trend_p10 & Effect_Size %in% c("Large", "Medium") ~ "Clinically Relevant Trend",
        Effect_Size == "Large" ~ "Large Effect Size Only",
        TRUE ~ "Limited Clinical Relevance"
      )
    )
  
  print(top_correlations %>%
          dplyr::select(OCTA_Parameter, Pearson_r, Pearson_p, 
                        Effect_Size, Clinical_Relevance, N))
  
  # 统计总结
  total_params <- nrow(membership_correlations)
  significant_count <- sum(membership_correlations$Significant_p05)
  trend_count <- sum(membership_correlations$Trend_p10)
  large_effect_count <- sum(membership_correlations$Effect_Size == "Large")
  
  cat("\n===== 统计总结 =====\n")
  cat("总分析参数:", total_params, "个\n")
  cat("显著相关 (p<0.05):", significant_count, "个 (", 
      round(significant_count/total_params*100, 1), "%)\n")
  cat("趋势显著 (p<0.10):", trend_count, "个 (", 
      round(trend_count/total_params*100, 1), "%)\n")
  cat("大效应量 (|r|≥0.5):", large_effect_count, "个 (", 
      round(large_effect_count/total_params*100, 1), "%)\n")
  
  # 最强相关的参数
  if(nrow(membership_correlations) > 0) {
    strongest_param <- membership_correlations$OCTA_Parameter[1]
    strongest_r <- round(membership_correlations$Pearson_r[1], 3)
    strongest_p <- round(membership_correlations$Pearson_p[1], 4)
    
    cat("\n最强相关参数:", strongest_param, "\n")
    cat("相关系数: r =", strongest_r, ", p =", strongest_p, "\n")
  }
  
} else {
  cat("❌ 未找到任何有效的相关性结果\n")
}

# ================== 5. 对比两种分析方法 ==================
cat("\n===== 分析方法对比 =====\n")
cat("方法1 - 时间窗口特异性分析（你刚才的结果）:\n")
cat("✓ 使用原始可穿戴设备指标值\n")
cat("✓ 发现8个显著相关 (p < 0.05)\n")
cat("✓ 最强相关: r = -0.704\n")
cat("✓ 提供了时间特异性信息\n\n")

cat("方法2 - Membership预测分析（当前分析）:\n")
if(nrow(membership_correlations) > 0) {
  membership_significant <- sum(membership_correlations$Significant_p05)
  membership_strongest <- max(abs(membership_correlations$Pearson_r))
  
  cat("• 使用聚类membership值\n")
  cat("• 发现", membership_significant, "个显著相关 (p < 0.05)\n")
  cat("• 最强相关: |r| =", round(membership_strongest, 3), "\n")
  cat("• 提供综合预测评分\n\n")
} else {
  cat("• 使用聚类membership值\n")
  cat("• 未发现显著相关\n")
  cat("• 可能需要优化聚类方法或增加样本量\n\n")
}

# ================== 6. 可视化membership相关性 ==================
if(nrow(membership_correlations) > 0 && any(membership_correlations$Trend_p10)) {
  # 创建散点图展示最强的几个相关性
  create_membership_scatterplots <- function(data, corr_results, top_n = 4) {
    plot_list <- list()
    
    top_params <- corr_results %>%
      filter(Trend_p10 == TRUE) %>%
      head(top_n)
    
    for(i in 1:nrow(top_params)) {
      param <- top_params$OCTA_Parameter[i]
      r_value <- round(top_params$Pearson_r[i], 3)
      p_value <- round(top_params$Pearson_p[i], 4)
      n_size <- top_params$N[i]
      
      plot_data <- data[!is.na(data$max_membership) & !is.na(data[[param]]), ]
      
      if(nrow(plot_data) >= 3) {
        param_clean <- gsub("_improvement", "", param)
        param_clean <- gsub("_", " ", param_clean)
        
        p <- ggplot(plot_data, aes(x = max_membership, y = .data[[param]])) +
          geom_point(size = 3, alpha = 0.7, color = "steelblue") +
          geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
          labs(
            title = param_clean,
            subtitle = paste0("r = ", r_value, ", p = ", p_value, ", n = ", n_size),
            x = "可穿戴设备 Membership",
            y = "OCTA 改善值"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 10)
          )
        
        plot_list[[i]] <- p
      }
    }
    
    if(length(plot_list) > 0) {
      combined_plot <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2))
      return(combined_plot)
    }
    
    return(NULL)
  }
  
  # 生成散点图
  membership_plots <- create_membership_scatterplots(membership_octa_analysis, membership_correlations)
  
  if(!is.null(membership_plots)) {
    cat("\n已生成membership相关性散点图\n")
  }
}

# ================== 7. 建议和结论 ==================
cat("\n===== 建议和结论 =====\n")

if(exists("membership_correlations") && nrow(membership_correlations) > 0) {
  membership_significant <- sum(membership_correlations$Significant_p05, na.rm = TRUE)
  
  if(membership_significant > 0) {
    cat("🎯 Membership分析成功!\n")
    cat("可穿戴设备聚类membership确实能预测OCTA改善\n")
  } else {
    cat("💡 Membership分析显示较弱的预测能力\n")
    cat("建议考虑以下策略:\n")
    cat("1. 结合时间窗口特异性和membership分析\n")
    cat("2. 使用复合评分方法\n")
    cat("3. 探索非线性关系\n")
    cat("4. 增加样本量验证结果\n")
  }
} else {
  cat("⚠️ 需要检查数据加载和处理流程\n")
}

cat("\n两种方法各有优势:\n")
cat("• 时间窗口分析: 提供机制洞察和时间特异性\n")
cat("• Membership分析: 提供简化的预测评分\n")
cat("• 建议在论文中同时报告两种方法的结果\n")

# 保存结果
write.csv(membership_correlations, "membership_octa_correlations_corrected.csv", row.names = FALSE)

cat("\n结果已保存到: membership_octa_correlations_corrected.csv\n")