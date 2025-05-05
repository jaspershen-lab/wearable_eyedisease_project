# -----------------------------------------------------
# 检查OCTA数据在各聚类中的情况
# -----------------------------------------------------
library(tidyverse)

# -----------------------------------------------------
# 1. 检查各聚类的OCTA样本量
# -----------------------------------------------------
check_octa_cluster_counts <- function(data_name, data) {
  # 确保cluster是因子
  data$cluster <- as.factor(data$cluster)
  
  # 计算每个聚类的样本量
  cluster_counts <- data %>%
    group_by(cluster) %>%
    summarise(count = n(), .groups = "drop")
  
  # 打印结果
  cat("\n========================================\n")
  cat(data_name, "中的聚类样本量:\n")
  cat("========================================\n")
  print(cluster_counts)
  
  # 返回样本量数据
  return(cluster_counts)
}

# 检查PPV组OCTA血流数据中的聚类样本量
ppv_bloodflow_counts <- check_octa_cluster_counts("PPV OCTA血流数据", ppv_bloodflow_analysis)

# 检查白内障组OCTA血流数据中的聚类样本量
cat_bloodflow_counts <- check_octa_cluster_counts("白内障 OCTA血流数据", cat_bloodflow_analysis)

# -----------------------------------------------------
# 2. 检查OCTA参数的缺失情况
# -----------------------------------------------------
check_octa_missing <- function(octa_data, param_pattern, data_name) {
  # 选择符合特定模式的OCTA参数
  octa_params <- names(octa_data)[grep(param_pattern, names(octa_data))]
  
  # 对于每个参数，计算缺失值的数量和百分比
  missing_summary <- data.frame(
    Parameter = character(),
    Missing_Count = integer(),
    Missing_Percent = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(param in octa_params) {
    missing_count <- sum(is.na(octa_data[[param]]))
    missing_percent <- missing_count / nrow(octa_data) * 100
    
    missing_summary <- missing_summary %>%
      rbind(data.frame(
        Parameter = param,
        Missing_Count = missing_count,
        Missing_Percent = missing_percent,
        stringsAsFactors = FALSE
      ))
  }
  
  # 按缺失百分比排序
  missing_summary <- missing_summary %>%
    arrange(desc(Missing_Percent))
  
  # 打印结果
  cat("\n========================================\n")
  cat(data_name, " - OCTA参数缺失情况:\n")
  cat("========================================\n")
  print(missing_summary)
  
  return(missing_summary)
}

# 检查PPV组OCTA血流参数中以"_0_21_improvement"结尾的缺失情况
ppv_bloodflow_missing <- check_octa_missing(ppv_bloodflow_analysis, "_0_21_improvement$", "PPV血流")

# 检查白内障组OCTA血流参数中以"_0_21_improvement"结尾的缺失情况
cat_bloodflow_missing <- check_octa_missing(cat_bloodflow_analysis, "_0_21_improvement$", "白内障血流")

# -----------------------------------------------------
# 3. 检查每个聚类中具体的OCTA参数值
# -----------------------------------------------------
check_octa_parameter_by_cluster <- function(octa_data, param_name, data_name) {
  # 确保参数存在
  if(!param_name %in% names(octa_data)) {
    cat("\n参数", param_name, "不存在于数据中\n")
    return(NULL)
  }
  
  # 确保cluster是因子
  octa_data$cluster <- as.factor(octa_data$cluster)
  
  # 提取ID、聚类和指定参数的值
  param_data <- octa_data %>%
    select(ID, cluster, !!sym(param_name)) %>%
    arrange(cluster, ID)
  
  # 打印每个聚类中的参数值
  cat("\n========================================\n")
  cat(data_name, " - ", param_name, "在各聚类中的值:\n")
  cat("========================================\n")
  print(param_data)
  
  # 计算每个聚类的基本统计量
  cluster_stats <- param_data %>%
    group_by(cluster) %>%
    summarise(
      count = n(),
      mean = mean(!!sym(param_name), na.rm = TRUE),
      median = median(!!sym(param_name), na.rm = TRUE),
      sd = sd(!!sym(param_name), na.rm = TRUE),
      min = min(!!sym(param_name), na.rm = TRUE),
      max = max(!!sym(param_name), na.rm = TRUE),
      missing = sum(is.na(!!sym(param_name))),
      .groups = "drop"
    )
  
  cat("\n聚类统计摘要:\n")
  print(cluster_stats)
  
  return(list(data = param_data, stats = cluster_stats))
}

# 检查PPV组中某个特定OCTA参数的情况
# 例如，检查VD_SVP参数
ppv_vd_svp_data <- check_octa_parameter_by_cluster(
  ppv_bloodflow_analysis, 
  "VD_SVP_0_21_improvement", 
  "PPV血流"
)

# 检查白内障组中的同一参数
cat_vd_svp_data <- check_octa_parameter_by_cluster(
  cat_bloodflow_analysis, 
  "VD_SVP_0_21_improvement", 
  "白内障血流"
)

# -----------------------------------------------------
# 4. 检查所有患者的OCTA数据可用情况
# -----------------------------------------------------
check_all_octa_data <- function(ppv_clusters_info, cat_clusters_info, octa_bloodflow, param_pattern) {
  # 合并PPV和白内障的聚类信息
  all_clusters <- bind_rows(
    ppv_clusters_info %>% mutate(surgery_type = "PPV"),
    cat_clusters_info %>% mutate(surgery_type = "Cataract")
  )
  
  # 提取OCTA参数名称
  octa_params <- names(octa_bloodflow)[grep(param_pattern, names(octa_bloodflow))]
  
  # 创建一个宽表格，显示每个ID的OCTA数据可用情况
  octa_availability <- all_clusters %>%
    select(ID, cluster, surgery_type)
  
  # 检查每个参数的可用情况
  for(param in octa_params) {
    # 从octa_bloodflow中提取ID和当前参数
    param_data <- octa_bloodflow %>%
      select(ID, !!sym(param))
    
    # 合并到可用性表格中
    octa_availability <- octa_availability %>%
      left_join(param_data, by = "ID") %>%
      rename(!!paste0(param, "_available") := !!sym(param))
  }
  
  # 添加总体可用性统计
  octa_availability <- octa_availability %>%
    rowwise() %>%
    mutate(
      available_params = sum(!is.na(c_across(ends_with("_available")))),
      total_params = length(octa_params),
      completeness = available_params / total_params
    ) %>%
    ungroup()
  
  # 按聚类和手术类型计算完整性摘要
  completeness_summary <- octa_availability %>%
    group_by(surgery_type, cluster) %>%
    summarise(
      patient_count = n(),
      mean_completeness = mean(completeness),
      fully_complete = sum(completeness == 1),
      partially_complete = sum(completeness > 0 & completeness < 1),
      no_data = sum(completeness == 0),
      .groups = "drop"
    )
  
  # 打印结果
  cat("\n========================================\n")
  cat("OCTA数据完整性摘要(按聚类和手术类型):\n")
  cat("========================================\n")
  print(completeness_summary)
  
  # 返回详细数据和摘要
  return(list(
    detail = octa_availability,
    summary = completeness_summary
  ))
}

# 检查所有患者的OCTA血流数据可用情况
all_octa_availability <- check_all_octa_data(
  ppv_cluster_info, 
  cat_cluster_info, 
  octa_bloodflow, 
  "_0_21_improvement$"
)

# -----------------------------------------------------
# 5. 保存结果到CSV文件
# -----------------------------------------------------
# 创建目录
dir.create("3_data_analysis/6_clustering_modeling/octa_completeness", 
           recursive = TRUE, showWarnings = FALSE)

# 保存OCTA可用性详细数据
write.csv(
  all_octa_availability$detail,
  "3_data_analysis/6_clustering_modeling/octa_completeness/octa_data_availability_detail.csv",
  row.names = FALSE
)

# 保存OCTA可用性摘要数据
write.csv(
  all_octa_availability$summary,
  "3_data_analysis/6_clustering_modeling/octa_completeness/octa_data_availability_summary.csv",
  row.names = FALSE
)
