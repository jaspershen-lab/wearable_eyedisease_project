# OCTA聚类与可穿戴设备聚类相关性分析
# 使用改进后的OCTA聚类结果

library(tidyverse)
library(ggplot2)
library(corrplot)
library(r4projects)
library(RColorBrewer)
library(pheatmap)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# 加载可穿戴设备聚类结果
wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/ppv_cluster_results_time_windows.csv")

# 加载改进后的OCTA聚类结果
octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/traditional_clustering/octa_cluster/ppv_octa_clean_only_results.csv")

# 创建相关性分析输出目录
dir.create("3_data_analysis/6_clustering_modeling/correlation_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/correlation_analysis")

# -------------------- 1. 加载聚类结果 --------------------


cat("可穿戴设备聚类数据:\n")
print(head(wearable_clusters))
cat("\nOCTA聚类数据:\n")
print(head(octa_clusters))

# 检查数据基本信息
cat("\n数据基本信息:\n")
cat("可穿戴设备聚类患者数:", nrow(wearable_clusters), "\n")
cat("OCTA聚类患者数:", nrow(octa_clusters), "\n")

cat("\n可穿戴设备聚类分布:\n")
print(table(wearable_clusters$max_cluster))
cat("\nOCTA聚类分布:\n")
print(table(octa_clusters$cluster))

# -------------------- 2. 数据匹配和清理 --------------------
cat("\n========== 数据匹配处理 ==========\n")

# 检查列名
cat("可穿戴设备聚类数据列名:", paste(colnames(wearable_clusters), collapse = ", "), "\n")
cat("OCTA聚类数据列名:", paste(colnames(octa_clusters), collapse = ", "), "\n")

# 先检查数据结构
cat("可穿戴设备数据结构:\n")
str(wearable_clusters)
cat("OCTA数据结构:\n")
str(octa_clusters)

# 合并两个聚类结果 - 分步进行以便调试
cat("开始数据合并...\n")

# 第一步：直接合并
merged_temp <- wearable_clusters %>%
  inner_join(octa_clusters, by = "subject_id")

cat("合并后的列名:", paste(colnames(merged_temp), collapse = ", "), "\n")

# 第二步：重命名（使用字符串而不是变量名）
# 使用更直接的方法重命名
names(merged_temp)[names(merged_temp) == "max_cluster"] <- "wearable_cluster"
names(merged_temp)[names(merged_temp) == "max_membership"] <- "wearable_membership"
names(merged_temp)[names(merged_temp) == "cluster"] <- "octa_cluster"
names(merged_temp)[names(merged_temp) == "method"] <- "octa_method"

# 选择需要的列
merged_data <- merged_temp[, c("subject_id", "wearable_cluster", "wearable_membership", 
                               "octa_cluster", "octa_method", "optimal_k", "silhouette_score")]

cat("匹配后的患者数:", nrow(merged_data), "\n")
cat("匹配的患者ID:", paste(merged_data$subject_id, collapse = ", "), "\n")


# 检查是否有患者在某个聚类中缺失
wearable_only <- setdiff(wearable_clusters$subject_id, merged_data$subject_id)
octa_only <- setdiff(octa_clusters$subject_id, merged_data$subject_id)

if(length(wearable_only) > 0) {
  cat("只在可穿戴设备聚类中的患者:", paste(wearable_only, collapse = ", "), "\n")
}
if(length(octa_only) > 0) {
  cat("只在OCTA聚类中的患者:", paste(octa_only, collapse = ", "), "\n")
}

# -------------------- 3. 描述性统计 --------------------
cat("\n========== 描述性统计 ==========\n")

# 创建交叉表
cross_table <- table(merged_data$wearable_cluster, merged_data$octa_cluster)
colnames(cross_table) <- paste0("OCTA_", colnames(cross_table))
rownames(cross_table) <- paste0("Wearable_", rownames(cross_table))

cat("聚类交叉表:\n")
print(cross_table)

# 计算百分比交叉表
cross_table_pct <- prop.table(cross_table) * 100
cat("\n聚类交叉表（百分比）:\n")
print(round(cross_table_pct, 1))

# 按行计算百分比（每个可穿戴设备聚类的OCTA聚类分布）
cross_table_row_pct <- prop.table(cross_table, margin = 1) * 100
cat("\n按行百分比（可穿戴设备聚类中OCTA聚类的分布）:\n")
print(round(cross_table_row_pct, 1))

# 按列计算百分比（每个OCTA聚类的可穿戴设备聚类分布）
cross_table_col_pct <- prop.table(cross_table, margin = 2) * 100
cat("\n按列百分比（OCTA聚类中可穿戴设备聚类的分布）:\n")
print(round(cross_table_col_pct, 1))

# -------------------- 4. 统计显著性检验 --------------------
cat("\n========== 统计显著性检验 ==========\n")

# 4.1 卡方检验
if(all(cross_table >= 1)) {  # 卡方检验要求期望频数≥1
  chi_test <- chisq.test(cross_table)
  cat("卡方检验结果:\n")
  cat("Chi-square =", round(chi_test$statistic, 4), "\n")
  cat("df =", chi_test$parameter, "\n")
  cat("p-value =", format(chi_test$p.value, scientific = TRUE), "\n")
  
  if(chi_test$p.value < 0.05) {
    cat("结论: 两种聚类方法显著相关 (p < 0.05)\n")
  } else {
    cat("结论: 两种聚类方法无显著相关性 (p ≥ 0.05)\n")
  }
} else {
  cat("卡方检验不适用（期望频数过小）\n")
}

# 4.2 Fisher精确检验（适用于小样本）
fisher_test <- fisher.test(cross_table, simulate.p.value = TRUE, B = 10000)
cat("\nFisher精确检验结果:\n")
cat("p-value =", format(fisher_test$p.value, scientific = TRUE), "\n")

if(fisher_test$p.value < 0.05) {
  cat("结论: 两种聚类方法显著相关 (p < 0.05)\n")
} else {
  cat("结论: 两种聚类方法无显著相关性 (p ≥ 0.05)\n")
}

# 4.3 Cramér's V系数（效应量）
cramers_v <- sqrt(chi_test$statistic / (sum(cross_table) * (min(dim(cross_table)) - 1)))
cat("\nCramér's V系数:", round(cramers_v, 3), "\n")
cat("效应量解释: ")
if(cramers_v < 0.1) {
  cat("微弱关联\n")
} else if(cramers_v < 0.3) {
  cat("弱关联\n")
} else if(cramers_v < 0.5) {
  cat("中等关联\n")
} else {
  cat("强关联\n")
}

# -------------------- 5. 一致性分析 --------------------
cat("\n========== 聚类一致性分析 ==========\n")

# 计算调整兰德指数（Adjusted Rand Index）
if(requireNamespace("mclust", quietly = TRUE)) {
  library(mclust)
  ari <- adjustedRandIndex(merged_data$wearable_cluster, merged_data$octa_cluster)
  cat("调整兰德指数(ARI):", round(ari, 3), "\n")
  cat("ARI解释: ")
  if(ari < 0) {
    cat("聚类结果比随机分组还差\n")
  } else if(ari < 0.2) {
    cat("聚类一致性很低\n")
  } else if(ari < 0.4) {
    cat("聚类一致性较低\n")
  } else if(ari < 0.6) {
    cat("聚类一致性中等\n")
  } else if(ari < 0.8) {
    cat("聚类一致性较高\n")
  } else {
    cat("聚类一致性很高\n")
  }
} else {
  cat("mclust包不可用，跳过ARI计算\n")
}

# 计算简单一致性
total_patients <- nrow(merged_data)
# 定义"一致"的标准（可以根据需要调整）
consistent_pairs <- 0

# 由于聚类数可能不同，我们检查是否主要聚类模式一致
# 这里采用简化的方法：检查主要的聚类模式
for(i in 1:total_patients) {
  for(j in (i+1):total_patients) {
    if(j <= total_patients) {
      # 检查两个患者在两种聚类中是否都被分在一起或都被分开
      same_wearable <- merged_data$wearable_cluster[i] == merged_data$wearable_cluster[j]
      same_octa <- merged_data$octa_cluster[i] == merged_data$octa_cluster[j]
      
      if(same_wearable == same_octa) {
        consistent_pairs <- consistent_pairs + 1
      }
    }
  }
}

total_pairs <- total_patients * (total_patients - 1) / 2
consistency_rate <- consistent_pairs / total_pairs
cat("简单一致性比例:", round(consistency_rate, 3), "\n")

# -------------------- 6. 可视化分析 --------------------
cat("\n========== 创建可视化图表 ==========\n")

# 6.1 交叉表热图
create_heatmap <- function(cross_table, title) {
  # 转换为数据框用于ggplot
  heatmap_data <- as.data.frame(as.table(cross_table))
  names(heatmap_data) <- c("Wearable_Cluster", "OCTA_Cluster", "Count")
  
  p <- ggplot(heatmap_data, aes(x = OCTA_Cluster, y = Wearable_Cluster, fill = Count)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Count), color = "black", size = 4, fontface = "bold") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "患者数") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    ) +
    labs(
      title = title,
      x = "OCTA聚类",
      y = "可穿戴设备聚类"
    )
  
  return(p)
}

# 创建热图
heatmap_plot <- create_heatmap(cross_table, "可穿戴设备聚类 vs OCTA聚类交叉表")
print(heatmap_plot)
ggsave("clustering_correlation_heatmap.pdf", heatmap_plot, width = 8, height = 6)
ggsave("clustering_correlation_heatmap.png", heatmap_plot, width = 8, height = 6, dpi = 300)

# 6.2 散点图（显示个体患者的聚类分配）
scatter_plot <- ggplot(merged_data, aes(x = factor(wearable_cluster), y = factor(octa_cluster))) +
  geom_jitter(aes(color = factor(wearable_cluster)), alpha = 0.7, size = 3, width = 0.2, height = 0.2) +
  geom_text(aes(label = subject_id), vjust = -0.8, size = 3) +
  scale_color_brewer(palette = "Set1", name = "可穿戴设备\n聚类") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(
    title = "患者聚类分配散点图",
    x = "可穿戴设备聚类",
    y = "OCTA聚类"
  )

print(scatter_plot)
ggsave("clustering_scatter_plot.pdf", scatter_plot, width = 10, height = 8)
ggsave("clustering_scatter_plot.png", scatter_plot, width = 10, height = 8, dpi = 300)

# 6.3 桑基图数据准备（如果需要的话）
sankey_data <- merged_data %>%
  count(wearable_cluster, octa_cluster) %>%
  rename(
    source = wearable_cluster,
    target = octa_cluster,
    value = n
  ) %>%
  mutate(
    source = paste0("Wearable_", source),
    target = paste0("OCTA_", target)
  )

cat("桑基图数据:\n")
print(sankey_data)

# -------------------- 7. 详细的患者分组分析 --------------------
cat("\n========== 详细患者分组分析 ==========\n")

# 按聚类组合分析患者
detailed_analysis <- merged_data %>%
  mutate(
    cluster_combination = paste0("W", wearable_cluster, "_O", octa_cluster)
  ) %>%
  group_by(cluster_combination, wearable_cluster, octa_cluster) %>%
  summarise(
    patient_count = n(),
    patients = paste(subject_id, collapse = ", "),
    avg_wearable_membership = mean(wearable_membership, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(patient_count))

cat("详细聚类组合分析:\n")
print(detailed_analysis)

# 保存详细分析结果
write.csv(detailed_analysis, "detailed_clustering_analysis.csv", row.names = FALSE)

# -------------------- 8. 总结报告 --------------------
cat("\n========== 聚类相关性分析总结报告 ==========\n")

# 保存总结结果
summary_results <- list(
  basic_info = list(
    total_patients = nrow(merged_data),
    wearable_clusters = length(unique(merged_data$wearable_cluster)),
    octa_clusters = length(unique(merged_data$octa_cluster))
  ),
  cross_table = cross_table,
  statistical_tests = list(
    chi_square_p = if(exists("chi_test")) chi_test$p.value else NA,
    fisher_p = fisher_test$p.value,
    cramers_v = cramers_v,
    ari = if(exists("ari")) ari else NA,
    consistency_rate = consistency_rate
  ),
  interpretation = list(
    significant_correlation = fisher_test$p.value < 0.05,
    effect_size = if(cramers_v < 0.1) "微弱" else if(cramers_v < 0.3) "弱" else if(cramers_v < 0.5) "中等" else "强",
    consistency_level = if(consistency_rate < 0.4) "低" else if(consistency_rate < 0.6) "中等" else "高"
  )
)

# 保存完整结果
write.csv(merged_data, "merged_clustering_results.csv", row.names = FALSE)

# 打印最终结论
cat("\n🔍 最终结论:\n")
cat("1. 总体相关性: ", ifelse(fisher_test$p.value < 0.05, "显著相关", "无显著相关性"), 
    " (Fisher精确检验 p = ", format(fisher_test$p.value, digits = 3), ")\n")
cat("2. 效应量: ", summary_results$interpretation$effect_size, "关联 (Cramér's V = ", 
    round(cramers_v, 3), ")\n")
cat("3. 聚类一致性: ", summary_results$interpretation$consistency_level, " (", 
    round(consistency_rate * 100, 1), "%)\n")

if(exists("ari")) {
  cat("4. 调整兰德指数: ", round(ari, 3), "\n")
}

cat("\n📊 主要发现:\n")
# 找出最大的聚类组合
max_combination <- detailed_analysis[1, ]
cat("- 最大的聚类组合: 可穿戴设备聚类", max_combination$wearable_cluster, 
    " + OCTA聚类", max_combination$octa_cluster, " (", max_combination$patient_count, "人)\n")
cat("- 涉及患者: ", max_combination$patients, "\n")

cat("\n📈 数据文件已保存:\n")
cat("- merged_clustering_results.csv: 合并的聚类结果\n")
cat("- detailed_clustering_analysis.csv: 详细分组分析\n")
cat("- clustering_correlation_heatmap.png: 相关性热图\n")
cat("- clustering_scatter_plot.png: 患者分布散点图\n")

cat("\n========== 相关性分析完成 ==========\n")