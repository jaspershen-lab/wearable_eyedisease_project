# -----------------------------------------------------
# Analyzing OCTA improvements by clusters for PPV and Cataract groups
# Focusing on 0_21 parameters only with ggbetweenstats visualization
# -----------------------------------------------------
library(r4projects)
library(tidyverse)
library(ggplot2)
library(rstatix)      # For statistical testing
library(ggpubr)       # For creating publication-ready plots
library(corrplot)     # For correlation visualizations
library(gtsummary)    # For creating summary tables
library(ggstatsplot)  # For statistical plots with significance indicators
library(rlang)        # For proper handling of non-standard evaluation
setwd(get_project_wd())
rm(list = ls())

# -----------------------------------------------------
# 1. Load the cluster results
# -----------------------------------------------------
# Assuming these files were created in your previous clustering analysis
ppv_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/ppv_cluster_results_time_windows.csv", check.names = FALSE)
cat_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint/cataract_cluster_results_time_windows.csv", check.names = FALSE)
# ppv_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/ppv_cluster_results_all_metrics.csv", check.names = FALSE)
# cat_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/cataract_cluster_results_all_metrics.csv", check.names = FALSE)
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")


# Extract the subject_id and cluster information
ppv_cluster_info <- data.frame(
  ID = ppv_clusters$subject_id,
  cluster = ppv_clusters$max_cluster,
  membership = ppv_clusters$max_membership
)

cat_cluster_info <- data.frame(
  ID = cat_clusters$subject_id,
  cluster = cat_clusters$max_cluster,
  membership = cat_clusters$max_membership
)

# -----------------------------------------------------
# 2. Load and prepare the OCTA improvement data
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
octa_bloodflow_features <- process_all_patients(octa_bloodflow_features)
octa_thickness_features <- process_all_patients(octa_thickness_features)

#-------------------------------
# 4. 计算OCTA指标的改善值(T2-T0)
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

# -----------------------------------------------------
# 3. Merge cluster information with OCTA improvement data
# -----------------------------------------------------
# For PPV group
ppv_bloodflow_analysis <- bloodflow_improvement %>%
  inner_join(ppv_cluster_info, by = "ID") %>%
  mutate(surgery_type = "PPV")

ppv_thickness_analysis <- thickness_improvement %>%
  inner_join(ppv_cluster_info, by = "ID") %>%
  mutate(surgery_type = "PPV")

# For Cataract group
cat_bloodflow_analysis <- bloodflow_improvement %>%
  inner_join(cat_cluster_info, by = "ID") %>%
  mutate(surgery_type = "Cataract")

cat_thickness_analysis <- thickness_improvement %>%
  inner_join(cat_cluster_info, by = "ID") %>%
  mutate(surgery_type = "Cataract")

# Combine data for overall analysis if needed
all_bloodflow_analysis <- bind_rows(ppv_bloodflow_analysis, cat_bloodflow_analysis)
all_thickness_analysis <- bind_rows(ppv_thickness_analysis, cat_thickness_analysis)

# -----------------------------------------------------
# 4. Filter and analyze only 0_21 parameters
# -----------------------------------------------------
# Function to filter only 0_21 parameters
filter_0_21_params <- function(data) {
  # Get column names that contain "0_21_improvement"
  params_0_21 <- names(data)[grep("0_21_improvement$", names(data))]
  
  # Keep only ID, cluster, membership, surgery_type and 0_21 parameters
  filtered_data <- data %>%
    dplyr::select(ID, cluster, membership, surgery_type, all_of(params_0_21))
  
  return(list(
    data = filtered_data,
    params = params_0_21
  ))
}

# Apply filter to each dataset
ppv_bloodflow_0_21 <- filter_0_21_params(ppv_bloodflow_analysis)
ppv_thickness_0_21 <- filter_0_21_params(ppv_thickness_analysis)
cat_bloodflow_0_21 <- filter_0_21_params(cat_bloodflow_analysis)
cat_thickness_0_21 <- filter_0_21_params(cat_thickness_analysis)


# -----------------------------------------------------
# 4.5 Filter to focus only on specific layers for bloodflow and thickness
# -----------------------------------------------------
# Function to filter bloodflow parameters to focus on SVP, ICP, DCP, and Choroid layers
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
  
  # 保留ID, cluster, membership, surgery_type和目标参数
  filtered_data <- data %>%
    dplyr::select(ID, cluster, membership, surgery_type, all_of(params_of_interest))
  
  return(list(
    data = filtered_data,
    params = params_of_interest
  ))
}

# Function to filter thickness parameters to focus on GCL.IPL, INL, and Retina layers
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
  
  # 保留ID, cluster, membership, surgery_type和目标参数
  filtered_data <- data %>%
    dplyr::select(ID, cluster, membership, surgery_type, all_of(params_of_interest))
  
  return(list(
    data = filtered_data,
    params = params_of_interest
  ))
}

# 过滤数据并将结果重新赋值给原始变量
# 先保存原始参数列表
original_ppv_bloodflow_params <- ppv_bloodflow_0_21$params
original_cat_bloodflow_params <- cat_bloodflow_0_21$params
original_ppv_thickness_params <- ppv_thickness_0_21$params
original_cat_thickness_params <- cat_thickness_0_21$params

# 应用血流过滤器
temp_ppv_bloodflow <- filter_bloodflow_layers(ppv_bloodflow_analysis)
temp_cat_bloodflow <- filter_bloodflow_layers(cat_bloodflow_analysis)

# 应用厚度过滤器
temp_ppv_thickness <- filter_thickness_layers(ppv_thickness_analysis)
temp_cat_thickness <- filter_thickness_layers(cat_thickness_analysis)

# 如果过滤后有参数，则更新原始对象
if(length(temp_ppv_bloodflow$params) > 0) {
  ppv_bloodflow_0_21 <- temp_ppv_bloodflow
  cat("更新了PPV血流数据集，从", length(original_ppv_bloodflow_params), "个参数减少到", length(ppv_bloodflow_0_21$params), "个参数\n")
}

if(length(temp_cat_bloodflow$params) > 0) {
  cat_bloodflow_0_21 <- temp_cat_bloodflow
  cat("更新了白内障血流数据集，从", length(original_cat_bloodflow_params), "个参数减少到", length(cat_bloodflow_0_21$params), "个参数\n")
}

if(length(temp_ppv_thickness$params) > 0) {
  ppv_thickness_0_21 <- temp_ppv_thickness
  cat("更新了PPV厚度数据集，从", length(original_ppv_thickness_params), "个参数减少到", length(ppv_thickness_0_21$params), "个参数\n")
}

if(length(temp_cat_thickness$params) > 0) {
  cat_thickness_0_21 <- temp_cat_thickness
  cat("更新了白内障厚度数据集，从", length(original_cat_thickness_params), "个参数减少到", length(cat_thickness_0_21$params), "个参数\n")
}

# 清理临时变量
rm(temp_ppv_bloodflow, temp_cat_bloodflow, temp_ppv_thickness, temp_cat_thickness)
rm(original_ppv_bloodflow_params, original_cat_bloodflow_params, original_ppv_thickness_params, original_cat_thickness_params)

# -----------------------------------------------------
# 5. Create visualization using ggbetweenstats
# -----------------------------------------------------
# Function to create plots with ggbetweenstats
create_betweenstats_plot <- function(data, param_col, title_prefix, surgery_type) {
  # Convert cluster to factor for better plotting
  data$cluster <- as.factor(data$cluster)
  
  # Create a more readable parameter name
  param_name <- gsub("_0_21_improvement", "", param_col)
  
  # Create the plot with ggbetweenstats using formula approach
  p <- ggstatsplot::ggbetweenstats(
    data = data,
    x = cluster,
    y = !!sym(param_col), # Use !!sym() for proper evaluation
    type = "parametric",
    pairwise.display = "all",  # Show all pairwise comparisons, not just significant ones
    p.adjust.method = "fdr",
    effsize.type = "unbiased",
    results.subtitle = TRUE,
    xlab = "Cluster",
    ylab = paste(param_name, "Improvement (T2-T0)"),
    title = paste0(title_prefix, " - ", param_name),
    subtitle = paste("Surgery Type:", surgery_type),
    centrality.plotting = TRUE,
    centrality.type = "parametric",
    centrality.point.args = list(size = 5, color = "darkred"),
    point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                      alpha = 0.4, size = 3, stroke = 0),
    boxplot.args = list(width = 0.3, alpha = 0.2),
    violin.args = list(width = 0.5, alpha = 0.2),
    ggsignif.args = list(textsize = 3, tip_length = 0.01),
    ggtheme = ggplot2::theme_bw(),  # Use theme_bw() instead of theme_ggstatsplot
    package = "RColorBrewer",
    palette = "Dark2"
  )
  
  return(p)
}

# Create directory for plots
dir.create("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0", 
           recursive = TRUE, showWarnings = FALSE)

# Generate plots for PPV blood flow 0_21 parameters
for (param in ppv_bloodflow_0_21$params) {
  p <- create_betweenstats_plot(
    ppv_bloodflow_0_21$data, 
    param, 
    "Blood Flow Improvement", 
    "PPV"
  )
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/ppv_bloodflow_", 
           gsub("_0_21_improvement", "", param), 
           ".pdf"),
    p,
    width = 10,
    height = 8
  )
}

# Generate plots for cataract blood flow 0_21 parameters
for (param in cat_bloodflow_0_21$params) {
  p <- create_betweenstats_plot(
    cat_bloodflow_0_21$data, 
    param, 
    "Blood Flow Improvement", 
    "Cataract"
  )
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/cat_bloodflow_", 
           gsub("_0_21_improvement", "", param), 
           ".pdf"),
    p,
    width = 10,
    height = 8
  )
}

# -----------------------------------------------------
# 6. Correlation analysis for 0_21 parameters
# -----------------------------------------------------
# Function to analyze correlation between parameters
analyze_parameter_correlations <- function(data, param_cols, title, output_prefix) {
  # Select only the parameters of interest
  corr_data <- data %>%
    dplyr::select(all_of(param_cols))
  
  # Calculate correlation matrix
  corr_matrix <- cor(corr_data, use = "pairwise.complete.obs", method = "spearman")
  
  # Create the output directory if it doesn't exist
  dir.create("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0", 
             recursive = TRUE, showWarnings = FALSE)
  
  # Create correlation plot
  pdf(
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/", 
           output_prefix, 
           "_correlation_plot.pdf"),
    width = 12,
    height = 10
  )
  
  corrplot(
    corr_matrix,
    method = "color",
    type = "upper",
    order = "hclust",
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    tl.cex = 0.7,
    number.cex = 0.7,
    col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
    title = title,
    mar = c(0, 0, 2, 0)
  )
  
  dev.off()
  
  return(corr_matrix)
}

# Analyze correlations for PPV blood flow 0_21 parameters
ppv_bloodflow_corr <- analyze_parameter_correlations(
  ppv_bloodflow_0_21$data,
  ppv_bloodflow_0_21$params,
  "PPV Blood Flow 0_21 Parameter Correlations",
  "ppv_bloodflow_0_21"
)

# Analyze correlations for cataract blood flow 0_21 parameters
cat_bloodflow_corr <- analyze_parameter_correlations(
  cat_bloodflow_0_21$data,
  cat_bloodflow_0_21$params,
  "Cataract Blood Flow 0_21 Parameter Correlations",
  "cat_bloodflow_0_21"
)

# -----------------------------------------------------
# 7. Statistical analysis of differences between clusters
# -----------------------------------------------------

# Function to test for differences between clusters
test_cluster_differences <- function(data, param_cols) {
  results_list <- list()
  all_p_values <- c()
  param_names <- c()
  
  # 先收集所有p值
  for (param in param_cols) {
    # ANOVA to test for differences between clusters
    model_formula <- as.formula(paste(param, "~ cluster"))
    anova_result <- data %>%
      anova_test(model_formula)
    
    # 收集p值
    all_p_values <- c(all_p_values, anova_result$p)
    param_names <- c(param_names, param)
    
    # Add results to the list
    results_list[[paste0(param, "_anova")]] <- anova_result
  }
  
  # 多重比较校正
  adjusted_p_values <- p.adjust(all_p_values, method = "fdr")
  
  # 添加校正后的p值到结果中
  for (i in 1:length(param_names)) {
    param <- param_names[i]
    results_list[[paste0(param, "_anova_adjusted")]] <- adjusted_p_values[i]
    
    # 仅对校正后仍然显著的结果进行事后检验
    if (adjusted_p_values[i] < 0.05) {
      model_formula <- as.formula(paste(param, "~ cluster"))
      # 修改这里：移除p.adjust.method参数，或设置为"none"
      posthoc <- data %>%
        pairwise_t_test(
          model_formula,
          p.adjust.method = "none"  # 不在这一步进行校正
        )
      
      # 单独收集该参数的所有post-hoc p值并进行校正
      posthoc$p.adj <- p.adjust(posthoc$p, method = "fdr")
      
      results_list[[paste0(param, "_posthoc")]] <- posthoc
    }
  }
  
  return(results_list)
}

# Test differences for PPV blood flow 0_21 parameters
ppv_bloodflow_tests <- test_cluster_differences(
  ppv_bloodflow_0_21$data,
  ppv_bloodflow_0_21$params
)

# Test differences for cataract blood flow 0_21 parameters
cat_bloodflow_tests <- test_cluster_differences(
  cat_bloodflow_0_21$data,
  cat_bloodflow_0_21$params
)

# Create summary table of significant results
create_significant_summary <- function(test_results, params, output_file) {
  sig_results <- data.frame(
    Parameter = character(),
    P_Value = numeric(),
    P_Adjusted = numeric(),  # 添加校正后的P值列
    Significant_Pairs = character(),
    stringsAsFactors = FALSE
  )
  
  for (param in params) {
    anova_result <- test_results[[paste0(param, "_anova")]]
    adjusted_p <- test_results[[paste0(param, "_anova_adjusted")]]
    
    if (!is.null(adjusted_p) && adjusted_p < 0.05) {  # 使用校正后的P值判断显著性
      param_name <- gsub("_0_21_improvement", "", param)
      
      # Get post-hoc results if available
      posthoc <- test_results[[paste0(param, "_posthoc")]]
      sig_pairs <- ""
      
      if (!is.null(posthoc)) {
        sig_pairs_df <- posthoc %>% 
          filter(p.adj < 0.05) %>%
          mutate(pair = paste(group1, "vs", group2))
        
        if (nrow(sig_pairs_df) > 0) {
          sig_pairs <- paste(sig_pairs_df$pair, collapse = "; ")
        }
      }
      
      # Add to summary
      sig_results <- sig_results %>%
        rbind(data.frame(
          Parameter = param_name,
          P_Value = anova_result$p,
          P_Adjusted = adjusted_p,  # 添加校正后的P值
          Significant_Pairs = sig_pairs,
          stringsAsFactors = FALSE
        ))
    }
  }
  
  # 按校正后的P值排序
  sig_results <- sig_results %>%
    arrange(P_Adjusted)
  
  # Save to file
  write.csv(sig_results, output_file, row.names = FALSE)
  
  return(sig_results)
}

# Create summary tables
ppv_sig_summary <- create_significant_summary(
  ppv_bloodflow_tests,
  ppv_bloodflow_0_21$params,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/ppv_bloodflow_0_21_significant_results.csv"
)

cat_sig_summary <- create_significant_summary(
  cat_bloodflow_tests,
  cat_bloodflow_0_21$params,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/cat_bloodflow_0_21_significant_results.csv"
)

# -----------------------------------------------------
# 8. Create multi-panel visualization of top significant parameters
# -----------------------------------------------------
# Function to create a combined plot of the most significant parameters
create_top_significant_plots <- function(data, test_results, params, title_prefix, max_params = 4) {
  # Get p-values for all parameters
  p_values <- sapply(params, function(param) {
    anova_result <- test_results[[paste0(param, "_anova")]]
    if (!is.null(anova_result)) {
      return(anova_result$p)
    } else {
      return(1) # Default to 1 if result not available
    }
  })
  
  # Sort parameters by p-value
  sorted_idx <- order(p_values)
  top_params <- params[sorted_idx[1:min(max_params, length(params))]]
  
  # Create a list to store plots
  plot_list <- list()
  
  # Generate plots for top parameters
  for (i in 1:length(top_params)) {
    param <- top_params[i]
    
    # Convert cluster to factor
    data$cluster <- as.factor(data$cluster)
    
    # Create a more readable parameter name
    param_name <- gsub("_0_21_improvement", "", param)
    
    # Create plot
    p <- ggstatsplot::ggbetweenstats(
      data = data,
      x = cluster,
      y = !!sym(param),
      type = "parametric",
      pairwise.display = "all",  # Show all pairwise comparisons
      p.adjust.method = "fdr",
      effsize.type = "unbiased",
      results.subtitle = TRUE,
      xlab = "Cluster",
      ylab = "Improvement (T2-T0)",
      title = param_name,
      centrality.plotting = TRUE,
      centrality.type = "parametric",
      centrality.point.args = list(size = 5, color = "darkred"),
      point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                        alpha = 0.4, size = 3, stroke = 0),
      boxplot.args = list(width = 0.3, alpha = 0.2),
      violin.args = list(width = 0.5, alpha = 0.2),
      ggsignif.args = list(textsize = 3, tip_length = 0.01),
      ggtheme = ggplot2::theme_bw(),  # Use theme_bw
      package = "RColorBrewer",
      palette = "Dark2"
    )
    
    plot_list[[i]] <- p
  }
  
  # Combine plots into a multi-panel figure
  if (length(plot_list) > 0) {
    combined_plot <- ggpubr::ggarrange(
      plotlist = plot_list,
      ncol = 2,
      nrow = ceiling(length(plot_list) / 2),
      common.legend = TRUE,
      legend = "bottom"
    )
    
    # Add an overall title
    combined_plot <- ggpubr::annotate_figure(
      combined_plot,
      top = ggpubr::text_grob(
        paste(title_prefix, "- Top Significant Parameters"),
        face = "bold",
        size = 16
      )
    )
    
    # Create directory if it doesn't exist
    dir.create("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0", 
               recursive = TRUE, showWarnings = FALSE)
    
    # Save the combined plot
    ggsave(
      paste0("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/", 
             gsub(" ", "_", tolower(title_prefix)), 
             "_top_significant_0_21.pdf"),
      combined_plot,
      width = 14,
      height = 12
    )
    
    return(combined_plot)
  } else {
    return(NULL)
  }
}

# Create combined plots for top significant parameters
ppv_top_plots <- create_top_significant_plots(
  ppv_bloodflow_0_21$data,
  ppv_bloodflow_tests,
  ppv_bloodflow_0_21$params,
  "PPV Blood Flow"
)

cat_top_plots <- create_top_significant_plots(
  cat_bloodflow_0_21$data,
  cat_bloodflow_tests,
  cat_bloodflow_0_21$params,
  "Cataract Blood Flow"
)

# -----------------------------------------------------
# 9. Correlation between cluster membership and OCTA improvement
# -----------------------------------------------------
# Function to analyze correlation between membership and improvement
analyze_membership_correlation <- function(data, param_cols) {
  results <- data.frame(
    Parameter = character(),
    Rho = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (param in param_cols) {
    # Spearman correlation test
    corr_test <- cor.test(data$membership, data[[param]], 
                          method = "spearman", exact = FALSE)
    
    # Add results to dataframe
    results <- results %>%
      rbind(data.frame(
        Parameter = gsub("_0_21_improvement", "", param),
        Rho = corr_test$estimate,
        P_value = corr_test$p.value,
        stringsAsFactors = FALSE
      ))
  }
  
  # 添加校正后的P值
  results$P_adjusted <- p.adjust(results$P_value, method = "holm")
  
  # 按校正后的P值排序
  results <- results %>%
    arrange(P_adjusted)
  
  return(results)
}

# Analyze correlation for PPV blood flow
ppv_membership_corr <- analyze_membership_correlation(
  ppv_bloodflow_0_21$data,
  ppv_bloodflow_0_21$params
)

# Analyze correlation for cataract blood flow
cat_membership_corr <- analyze_membership_correlation(
  cat_bloodflow_0_21$data,
  cat_bloodflow_0_21$params
)

# Save correlation results
write.csv(
  ppv_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/ppv_membership_correlation_0_21.csv",
  row.names = FALSE
)

write.csv(
  cat_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/cat_membership_correlation_0_21.csv",
  row.names = FALSE
)


# Analyze correlation for PPV thickness
ppv_thickness_membership_corr <- analyze_membership_correlation(
  ppv_thickness_0_21$data,
  ppv_thickness_0_21$params
)

# Analyze correlation for cataract thickness
cat_thickness_membership_corr <- analyze_membership_correlation(
  cat_thickness_0_21$data,
  cat_thickness_0_21$params
)

# Save correlation results for thickness
write.csv(
  ppv_thickness_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1w/T1_T0/ppv_thickness_membership_correlation_0_21.csv",
  row.names = FALSE
)

write.csv(
  cat_thickness_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1w/T1_T0/cat_thickness_membership_correlation_0_21.csv",
  row.names = FALSE
)

# -----------------------------------------------------
# 10. Save processed data for further analysis
# -----------------------------------------------------
# Save the filtered datasets
saveRDS(ppv_bloodflow_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/ppv_bloodflow_0_21_analysis.rds")
saveRDS(cat_bloodflow_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/cat_bloodflow_0_21_analysis.rds")
saveRDS(ppv_thickness_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/ppv_thickness_0_21_analysis.rds")
saveRDS(cat_thickness_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/cat_thickness_0_21_analysis.rds")

# Print a summary of the analysis
cat("\n=========================================================\n")
cat("OCTA 0_21 Parameter Analysis Complete\n")
cat("=========================================================\n")

# Print summary of significant parameters
cat("\nPPV Blood Flow - Significant 0_21 Parameters (p < 0.05):\n")
if(nrow(ppv_sig_summary) > 0) {
  print(ppv_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

cat("\nCataract Blood Flow - Significant 0_21 Parameters (p < 0.05):\n")
if(nrow(cat_sig_summary) > 0) {
  print(cat_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

# Print summary of top membership correlations
cat("\nTop PPV Membership Correlations:\n")
print(head(ppv_membership_corr, 5))

cat("\nTop Cataract Membership Correlations:\n")
print(head(cat_membership_corr, 5))




# -----------------------------------------------------
# Adding Thickness Analysis Code
# -----------------------------------------------------

# 1. Generate plots for PPV thickness 0_21 parameters
for (param in ppv_thickness_0_21$params) {
  p <- create_betweenstats_plot(
    ppv_thickness_0_21$data, 
    param, 
    "Thickness Improvement", 
    "PPV"
  )
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/ppv_thickness_", 
           gsub("_0_21_improvement", "", param), 
           ".pdf"),
    p,
    width = 10,
    height = 8
  )
}

# 2. Generate plots for cataract thickness 0_21 parameters
for (param in cat_thickness_0_21$params) {
  p <- create_betweenstats_plot(
    cat_thickness_0_21$data, 
    param, 
    "Thickness Improvement", 
    "Cataract"
  )
  
  # Save the plot
  ggsave(
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/cat_thickness_", 
           gsub("_0_21_improvement", "", param), 
           ".pdf"),
    p,
    width = 10,
    height = 8
  )
}

# 3. Analyze correlations for PPV thickness 0_21 parameters
ppv_thickness_corr <- analyze_parameter_correlations(
  ppv_thickness_0_21$data,
  ppv_thickness_0_21$params,
  "PPV Thickness 0_21 Parameter Correlations",
  "ppv_thickness_0_21"
)

# 4. Analyze correlations for cataract thickness 0_21 parameters
cat_thickness_corr <- analyze_parameter_correlations(
  cat_thickness_0_21$data,
  cat_thickness_0_21$params,
  "Cataract Thickness 0_21 Parameter Correlations",
  "cat_thickness_0_21"
)

# 5. Test differences for PPV thickness 0_21 parameters
ppv_thickness_tests <- test_cluster_differences(
  ppv_thickness_0_21$data,
  ppv_thickness_0_21$params
)

# 6. Test differences for cataract thickness 0_21 parameters
cat_thickness_tests <- test_cluster_differences(
  cat_thickness_0_21$data,
  cat_thickness_0_21$params
)

# 7. Create summary tables for thickness data
ppv_thickness_sig_summary <- create_significant_summary(
  ppv_thickness_tests,
  ppv_thickness_0_21$params,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/ppv_thickness_0_21_significant_results.csv"
)

cat_thickness_sig_summary <- create_significant_summary(
  cat_thickness_tests,
  cat_thickness_0_21$params,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/cat_thickness_0_21_significant_results.csv"
)

# 8. Create combined plots for top significant thickness parameters
ppv_thickness_top_plots <- create_top_significant_plots(
  ppv_thickness_0_21$data,
  ppv_thickness_tests,
  ppv_thickness_0_21$params,
  "PPV Thickness"
)

cat_thickness_top_plots <- create_top_significant_plots(
  cat_thickness_0_21$data,
  cat_thickness_tests,
  cat_thickness_0_21$params,
  "Cataract Thickness"
)

# 9. Analyze correlation between membership and thickness improvement
ppv_thickness_membership_corr <- analyze_membership_correlation(
  ppv_thickness_0_21$data,
  ppv_thickness_0_21$params
)

cat_thickness_membership_corr <- analyze_membership_correlation(
  cat_thickness_0_21$data,
  cat_thickness_0_21$params
)

# 10. Save correlation results for thickness data
write.csv(
  ppv_thickness_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/ppv_thickness_membership_correlation_0_21.csv",
  row.names = FALSE
)

write.csv(
  cat_thickness_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/all_metrics/1m/T2_T0/cat_thickness_membership_correlation_0_21.csv",
  row.names = FALSE
)

# 11. Print summary of thickness analysis results
# Print summary of significant parameters
cat("\nPPV Blood Flow - Significant 0_21 Parameters (p < 0.05):\n")
if(nrow(ppv_sig_summary) > 0) {
  print(ppv_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

cat("\nCataract Blood Flow - Significant 0_21 Parameters (p < 0.05):\n")
if(nrow(cat_sig_summary) > 0) {
  print(cat_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

cat("\n=========================================================\n")
cat("OCTA Thickness 0_21 Parameter Analysis Complete\n")
cat("=========================================================\n")

# Print summary of significant parameters for thickness
cat("\nPPV Thickness - Significant 0_21 Parameters (p < 0.05):\n")
if(nrow(ppv_thickness_sig_summary) > 0) {
  print(ppv_thickness_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

cat("\nCataract Thickness - Significant 0_21 Parameters (p < 0.05):\n")
if(nrow(cat_thickness_sig_summary) > 0) {
  print(cat_thickness_sig_summary)
} else {
  cat("No significant parameters found.\n")
}

# Print summary of top membership correlations
cat("\nTop PPV Membership Correlations:\n")
print(head(ppv_membership_corr, 5))

cat("\nTop Cataract Membership Correlations:\n")
print(head(cat_membership_corr, 5))


# Print summary of top membership correlations for thickness
cat("\nTop PPV Thickness Membership Correlations:\n")
print(head(ppv_thickness_membership_corr, 5))

cat("\nTop Cataract Thickness Membership Correlations:\n")
print(head(cat_thickness_membership_corr, 5))

