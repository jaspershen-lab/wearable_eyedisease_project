library(tidyverse)
library(lubridate)
library(Biobase)
library(Mfuzz)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

####read data
ppv_data<- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/combined_D_Surg1.csv", check.names = FALSE)
cat_data<-read.csv("3_data_analysis/6_clustering_modeling/data_prepare/combined_NoD_Surg0.csv", check.names = FALSE)

#####
str(ppv_data)
str(cat_data)

dir.create("3_data_analysis/6_clustering_modeling/mfuzz/rhr_mean",recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/rhr_mean")

# -------------------- 1. 分析PPV组的NA值 --------------------
# 选择术前-4到-1天的变量
ppv_preop <- ppv_data %>% 
  select(subject_id, `day_-7_mean_rhr_1`,`day_-6_mean_rhr_1`,`day_-5_mean_rhr_1`,`day_-4_mean_rhr_1`, `day_-3_mean_rhr_1`, `day_-2_mean_rhr_1`, `day_-1_mean_rhr_1`)

# 统计每个时间点的NA数量和比例
ppv_na_by_timepoint <- ppv_preop %>%
  summarise(across(starts_with("day_"), ~sum(is.na(.)), .names = "{.col}_na"),
            across(starts_with("day_"), ~mean(is.na(.)) * 100, .names = "{.col}_pct"))

# 转换为长格式以便于可视化
ppv_na_long <- ppv_na_by_timepoint %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "metric"),
               names_pattern = "(.*)_(.*)",
               values_to = "value")

# 统计每个患者的NA数量
ppv_na_by_subject <- ppv_preop %>%
  rowwise() %>%
  mutate(na_count = sum(is.na(c_across(starts_with("day_"))))) %>%
  ungroup()

# 统计不同NA数量的患者分布
ppv_na_distribution <- ppv_na_by_subject %>%
  count(na_count) %>%
  mutate(percentage = n / sum(n) * 100)

# -------------------- 2. 分析白内障组的NA值 --------------------
# 选择术前-4到-1天的变量
cat_preop <- cat_data %>% 
  select(subject_id, `day_-7_mean_rhr_1`,`day_-6_mean_rhr_1`,`day_-5_mean_rhr_1`,`day_-4_mean_rhr_1`, `day_-3_mean_rhr_1`, `day_-2_mean_rhr_1`, `day_-1_mean_rhr_1`)

# 统计每个时间点的NA数量和比例
cat_na_by_timepoint <- cat_preop %>%
  summarise(across(starts_with("day_"), ~sum(is.na(.)), .names = "{.col}_na"),
            across(starts_with("day_"), ~mean(is.na(.)) * 100, .names = "{.col}_pct"))

# 统计每个患者的NA数量
cat_na_by_subject <- cat_preop %>%
  rowwise() %>%
  mutate(na_count = sum(is.na(c_across(starts_with("day_"))))) %>%
  ungroup()

# 统计不同NA数量的患者分布
cat_na_distribution <- cat_na_by_subject %>%
  count(na_count) %>%
  mutate(percentage = n / sum(n) * 100)

# -------------------- 3. 输出结果 --------------------
# 打印PPV组的NA值统计
cat("\n============ PPV组NA值分析 ============\n")
cat("样本量:", nrow(ppv_preop), "\n\n")

cat("每个时间点的NA值数量:\n")
ppv_na_cols <- grep("_na$", names(ppv_na_by_timepoint), value = TRUE)
for (col in ppv_na_cols) {
  day <- gsub("_na$", "", col)
  cat(sprintf("%s: %d NA (%.2f%%)\n", 
              day, 
              ppv_na_by_timepoint[[col]], 
              ppv_na_by_timepoint[[paste0(day, "_pct")]]))
}

cat("\n患者的NA值分布:\n")
print(ppv_na_distribution)

# 打印白内障组的NA值统计
cat("\n============ 白内障组NA值分析 ============\n")
cat("样本量:", nrow(cat_preop), "\n\n")

cat("每个时间点的NA值数量:\n")
cat_na_cols <- grep("_na$", names(cat_na_by_timepoint), value = TRUE)
for (col in cat_na_cols) {
  day <- gsub("_na$", "", col)
  cat(sprintf("%s: %d NA (%.2f%%)\n", 
              day, 
              cat_na_by_timepoint[[col]], 
              cat_na_by_timepoint[[paste0(day, "_pct")]]))
}

cat("\n患者的NA值分布:\n")
print(cat_na_distribution)

# -------------------- 4. 可视化NA值分布 --------------------
# 创建患者NA值热图的函数
plot_na_heatmap <- function(data, title) {
  # 准备数据
  na_matrix <- data %>%
    select(starts_with("day_")) %>%
    mutate(subject_id = row_number()) %>%
    pivot_longer(cols = starts_with("day_"), 
                 names_to = "day", 
                 values_to = "value") %>%
    mutate(is_na = is.na(value),
           day = factor(day, levels = c("day_-7_mean_rhr_1", 
                                        "day_-6_mean_rhr_1", 
                                        "day_-5_mean_rhr_1", 
                                        "day_-4_mean_rhr_1", 
                                        "day_-3_mean_rhr_1", 
                                        "day_-2_mean_rhr_1", 
                                        "day_-1_mean_rhr_1")),
           day_label = gsub("day_(.*)_mean_rhr_1", "Day \\1", day))
  
  # 绘制热图
  ggplot(na_matrix, aes(x = day_label, y = subject_id, fill = is_na)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "red"), 
                      labels = c("FALSE" = "no NA", "TRUE" = "NA")) +
    labs(title = title,
         x = "presurgery days",
         y = "Subject ID",
         fill = "Data status") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
}

# 绘制PPV组的NA热图
ppv_heatmap <- plot_na_heatmap(ppv_preop, "PPV NA data status")
print(ppv_heatmap)

# 绘制白内障组的NA热图
cat_heatmap <- plot_na_heatmap(cat_preop, "Cataract NA data status")
print(cat_heatmap)

# 保存图表
ggsave("ppv_na_heatmap.png", ppv_heatmap, width = 8, height = 6)
ggsave("cat_na_heatmap.png", cat_heatmap, width = 8, height = 8)

# -------------------- 5. 统计缺失值模式 --------------------
# 统计整体NA值比例
total_na_ppv <- sum(is.na(ppv_preop %>% select(starts_with("day_")))) / 
  (nrow(ppv_preop) * (ncol(ppv_preop) - 1)) * 100

total_na_cat <- sum(is.na(cat_preop %>% select(starts_with("day_")))) / 
  (nrow(cat_preop) * (ncol(cat_preop) - 1)) * 100

cat("\n============ 总体NA值比例 ============\n")
cat(sprintf("PPV组总体NA值比例: %.2f%%\n", total_na_ppv))
cat(sprintf("白内障组总体NA值比例: %.2f%%\n", total_na_cat))

# 检查是否有任何患者在-4到-1天中有3个或更多NA值的患者
high_na_ppv <- ppv_na_by_subject %>% 
  filter(na_count >= 3) %>% 
  nrow()

high_na_cat <- cat_na_by_subject %>% 
  filter(na_count >= 3) %>% 
  nrow()

cat("\n有3个或更多NA值的患者数量:\n")
cat(sprintf("PPV组: %d\n", high_na_ppv))
cat(sprintf("白内障组: %d\n", high_na_cat))


# 分别对PPV和白内障手术组进行处理
# 选择术前-4到-1天的数据
ppv_subset <- ppv_data[, c("subject_id", "surgery_date_rhr", 
                           "day_-4_mean_rhr_1", "day_-3_mean_rhr_1", 
                           "day_-2_mean_rhr_1", "day_-1_mean_rhr_1")]

cat_subset <- cat_data[, c("subject_id", "surgery_date_rhr", 
                           "day_-4_mean_rhr_1", "day_-3_mean_rhr_1", 
                           "day_-2_mean_rhr_1", "day_-1_mean_rhr_1")]

# -----------------------------------------------------
# 处理NA值的几种方法（分别对两组进行处理）
# -----------------------------------------------------

# # 方法1：移除包含NA的行
# ppv_complete <- na.omit(ppv_subset)
# cat_complete <- na.omit(cat_subset)
# cat("PPV组使用complete cases后剩余样本数:", nrow(ppv_complete), "\n")
# cat("白内障组使用complete cases后剩余样本数:", nrow(cat_complete), "\n")

# 方法2：用每列的均值填充NA值（推荐方法，保留更多样本）
ppv_imputed <- ppv_subset
cat_imputed <- cat_subset

# PPV组NA值填充
for(i in 3:6) {  # 第3到6列是day_-4到day_-1
  ppv_imputed[is.na(ppv_imputed[,i]), i] <- mean(ppv_imputed[,i], na.rm=TRUE)
}

# 白内障组NA值填充
for(i in 3:6) {  # 第3到6列是day_-4到day_-1
  cat_imputed[is.na(cat_imputed[,i]), i] <- mean(cat_imputed[,i], na.rm=TRUE)
}

# -----------------------------------------------------
# 准备Mfuzz聚类的数据
# -----------------------------------------------------

prep_data_for_mfuzz <- function(data) {
  # 提取术前(-4到-1)的数据
  rhr_data <- data %>%
    select(subject_id, starts_with("day_-")) %>%
    column_to_rownames("subject_id")
  
  # 转置数据以匹配Mfuzz输入格式要求
  rhr_matrix <- as.matrix(rhr_data)
  
  # 创建ExpressionSet对象
  eset <- ExpressionSet(assayData = rhr_matrix)
  
  return(eset)
}

# -----------------------------------------------------
# 为每个组分别执行Mfuzz聚类
# -----------------------------------------------------

# 1. PPV组聚类分析
# -----------------
# 准备PPV组数据
ppv_eset <- prep_data_for_mfuzz(ppv_imputed)
ppv_eset_std <- standardise(ppv_eset)

# 确定PPV组的最佳模糊系数m
ppv_m <- mestimate(ppv_eset_std)
cat("PPV组估计的最佳模糊系数 m:", ppv_m, "\n")

# 确定最佳聚类数的函数
optimal_cluster <- function(eset_std, max_c = 10, m = 1.25) {
  # 计算不同聚类数的Dmin (最小质心间距离)
  dmin_values <- numeric(max_c - 1)
  
  for (c in 2:max_c) {
    cl <- mfuzz(eset_std, c = c, m = m)
    centers <- cl$centers
    
    # 计算所有质心之间的最小距离
    min_dist <- Inf
    for (i in 1:(c-1)) {
      for (j in (i+1):c) {
        dist_ij <- sqrt(sum((centers[i,] - centers[j,])^2))
        min_dist <- min(min_dist, dist_ij)
      }
    }
    
    dmin_values[c-1] <- min_dist
  }
  
  # 绘制Dmin图
  plot(2:max_c, dmin_values, type = "b", 
       xlab = "聚类数", ylab = "最小质心间距离",
       main = "最小质心间距离vs聚类数")
  
  # 返回建议的聚类数 (Dmin急剧下降后的点)
  return(dmin_values)
}

# 运行以确定PPV组的最佳聚类数（可选步骤）
cat("计算PPV组的最佳聚类数...\n")
ppv_dmin_values <- optimal_cluster(ppv_eset_std, max_c = 8, m = ppv_m)

# 执行PPV组聚类 (假设我们选择了3个聚类)
ppv_n_clusters <- 3  # 根据optimal_cluster结果调整
ppv_cl <- mfuzz(ppv_eset_std, c = ppv_n_clusters, m = ppv_m)


# 可视化PPV组聚类结果
pdf("ppv_cluster_results.pdf", width = 10, height = 8)
mfuzz.plot(ppv_eset_std, cl = ppv_cl, mfrow = c(2, 2), 
           time.labels = c("day_-4", "day_-3", "day_-2", "day_-1"),
           new.window = FALSE)
dev.off()

# 提取PPV组的聚类成员度
ppv_membership <- ppv_cl$membership
rownames(ppv_membership) <- rownames(ppv_eset_std)
ppv_membership_df <- as.data.frame(ppv_membership)
colnames(ppv_membership_df) <- paste0("Cluster", 1:ppv_n_clusters)

# 将成员度与原始数据合并
ppv_results <- cbind(
  ppv_imputed,
  max_cluster = apply(ppv_membership, 1, which.max),
  max_membership = apply(ppv_membership, 1, max)
)

# 2. 白内障组聚类分析
# --------------------
# 准备白内障组数据
cat_eset <- prep_data_for_mfuzz(cat_imputed)
cat_eset_std <- standardise(cat_eset)

# 确定白内障组的最佳模糊系数m
cat_m <- mestimate(cat_eset_std)
cat("白内障组估计的最佳模糊系数 m:", cat_m, "\n")

# 运行以确定白内障组的最佳聚类数（可选步骤）
cat("计算白内障组的最佳聚类数...\n")
cat_dmin_values <- optimal_cluster(cat_eset_std, max_c = 8, m = cat_m)

# 执行白内障组聚类 (假设我们选择了3个聚类)
cat_n_clusters <- 3  # 根据optimal_cluster结果调整
cat_cl <- mfuzz(cat_eset_std, c = cat_n_clusters, m = cat_m)

# 可视化白内障组聚类结果
pdf("cataract_cluster_results.pdf", width = 10, height = 8)
mfuzz.plot(cat_eset_std, cl = cat_cl, mfrow = c(2, 2), 
           time.labels = c("day_-4", "day_-3", "day_-2", "day_-1"),
           new.window = FALSE)
dev.off()


# 提取白内障组的聚类成员度
cat_membership <- cat_cl$membership
rownames(cat_membership) <- rownames(cat_eset_std)
cat_membership_df <- as.data.frame(cat_membership)
colnames(cat_membership_df) <- paste0("Cluster", 1:cat_n_clusters)

# 将成员度与原始数据合并
cat_results <- cbind(
  cat_imputed,
  max_cluster = apply(cat_membership, 1, which.max),
  max_membership = apply(cat_membership, 1, max)
)

# -----------------------------------------------------
# 使用ggplot2可视化每个组的聚类结果
# -----------------------------------------------------

# 1. 可视化PPV组聚类结果
# ----------------------
ppv_cluster_summary <- ppv_results %>%
  group_by(max_cluster) %>%
  summarise(
    day_minus_4 = mean(`day_-4_mean_rhr_1`, na.rm = TRUE),
    day_minus_3 = mean(`day_-3_mean_rhr_1`, na.rm = TRUE),
    day_minus_2 = mean(`day_-2_mean_rhr_1`, na.rm = TRUE),
    day_minus_1 = mean(`day_-1_mean_rhr_1`, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = starts_with("day_"),
    names_to = "day",
    values_to = "rhr"
  ) %>%
  mutate(
    day = case_when(
      day == "day_minus_4" ~ -4,
      day == "day_minus_3" ~ -3,
      day == "day_minus_2" ~ -2,
      day == "day_minus_1" ~ -1
    )
  )

# 绘制PPV组聚类趋势图
ppv_plot <- ggplot(ppv_cluster_summary, 
                   aes(x = day, y = rhr, color = factor(max_cluster), group = max_cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1)) +
  labs(
    title = "The clustering trend of preoperative mean RHR in PPV group",
    x = "prior surgery days",
    y = "mean RHR",
    color = "cluster"
  ) +
  theme_bw()

ppv_plot

# 保存PPV组趋势图
ggsave("ppv_cluster_trends_mean_rhr.pdf", ppv_plot, width = 10, height = 6)

# 2. 可视化白内障组聚类结果
# -------------------------
cat_cluster_summary <- cat_results %>%
  group_by(max_cluster) %>%
  summarise(
    day_minus_4 = mean(`day_-4_mean_rhr_1`, na.rm = TRUE),
    day_minus_3 = mean(`day_-3_mean_rhr_1`, na.rm = TRUE),
    day_minus_2 = mean(`day_-2_mean_rhr_1`, na.rm = TRUE),
    day_minus_1 = mean(`day_-1_mean_rhr_1`, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = starts_with("day_"),
    names_to = "day",
    values_to = "rhr"
  ) %>%
  mutate(
    day = case_when(
      day == "day_minus_4" ~ -4,
      day == "day_minus_3" ~ -3,
      day == "day_minus_2" ~ -2,
      day == "day_minus_1" ~ -1
    )
  )

# 绘制白内障组聚类趋势图
cat_plot <- ggplot(cat_cluster_summary, 
                   aes(x = day, y = rhr, color = factor(max_cluster), group = max_cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1)) +
  labs(
    title = "The clustering trend of preoperative mean RHR in Cataract group",
    x = "prior surgery days",
    y = "mean RHR",
    color = "cluster"
  ) +
  theme_bw()
cat_plot
# 保存白内障组趋势图
ggsave("cataract_cluster_trends_mean_RHR.pdf", cat_plot, width = 10, height = 6)

# -----------------------------------------------------
# 输出详细的聚类成员度信息
# -----------------------------------------------------
# PPV组
ppv_details <- ppv_results %>%
  select(subject_id, max_cluster, max_membership) %>%
  arrange(max_cluster, desc(max_membership))

print("PPV组的聚类结果:")
print(ppv_details)

# 白内障组
cat_details <- cat_results %>%
  select(subject_id, max_cluster, max_membership) %>%
  arrange(max_cluster, desc(max_membership))

print("白内障组的聚类结果:")
print(cat_details)

# -----------------------------------------------------
# 保存聚类结果
# -----------------------------------------------------
write.csv(ppv_results, "ppv_cluster_results_mean_RHR.csv", row.names = FALSE)
write.csv(cat_results, "cataract_cluster_results_mean_RHR.csv", row.names = FALSE)

# -----------------------------------------------------
# 分析PPV组和白内障组的聚类模式差异
# -----------------------------------------------------
# 计算并比较两组各聚类的平均心率值
cat("\nPPV组各聚类平均心率:\n")
print(ppv_cluster_summary %>% 
        group_by(max_cluster) %>% 
        summarise(mean_rhr = mean(rhr)))

cat("\n白内障组各聚类平均心率:\n")
print(cat_cluster_summary %>% 
        group_by(max_cluster) %>% 
        summarise(mean_rhr = mean(rhr)))