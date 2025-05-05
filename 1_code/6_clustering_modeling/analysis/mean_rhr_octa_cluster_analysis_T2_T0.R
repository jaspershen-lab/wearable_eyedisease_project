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
ppv_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/rhr_mean/ppv_cluster_results_mean_RHR.csv", check.names = FALSE)
cat_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/rhr_mean/cataract_cluster_results_mean_RHR.csv", check.names = FALSE)
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

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
    select(ID, cluster, membership, surgery_type, all_of(params_0_21))
  
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
    p.adjust.method = "holm",
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
dir.create("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0", 
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
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/plots_0_21/ppv_bloodflow_", 
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
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/plots_0_21/cat_bloodflow_", 
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
    select(all_of(param_cols))
  
  # Calculate correlation matrix
  corr_matrix <- cor(corr_data, use = "pairwise.complete.obs", method = "spearman")
  
  # Create the output directory if it doesn't exist
  dir.create("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0", 
             recursive = TRUE, showWarnings = FALSE)
  
  # Create correlation plot
  pdf(
    paste0("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/", 
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
  
  for (param in param_cols) {
    # ANOVA to test for differences between clusters
    model_formula <- as.formula(paste(param, "~ cluster"))
    anova_result <- data %>%
      anova_test(model_formula)
    
    # Add results to the list
    results_list[[paste0(param, "_anova")]] <- anova_result
    
    # If ANOVA is significant, add post-hoc tests
    if (anova_result$p < 0.05) {
      posthoc <- data %>%
        pairwise_t_test(
          model_formula,
          p.adjust.method = "holm"
        )
      
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
    Significant_Pairs = character(),
    stringsAsFactors = FALSE
  )
  
  for (param in params) {
    anova_result <- test_results[[paste0(param, "_anova")]]
    
    if (!is.null(anova_result) && anova_result$p < 0.05) {
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
          Significant_Pairs = sig_pairs,
          stringsAsFactors = FALSE
        ))
    }
  }
  
  # Sort by p-value
  sig_results <- sig_results %>%
    arrange(P_Value)
  
  # Save to file
  write.csv(sig_results, output_file, row.names = FALSE)
  
  return(sig_results)
}

# Create summary tables
ppv_sig_summary <- create_significant_summary(
  ppv_bloodflow_tests,
  ppv_bloodflow_0_21$params,
  "3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/ppv_bloodflow_0_21_significant_results.csv"
)

cat_sig_summary <- create_significant_summary(
  cat_bloodflow_tests,
  cat_bloodflow_0_21$params,
  "3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/cat_bloodflow_0_21_significant_results.csv"
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
      p.adjust.method = "holm",
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
    dir.create("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0", 
               recursive = TRUE, showWarnings = FALSE)
    
    # Save the combined plot
    ggsave(
      paste0("3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/", 
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
  
  # Sort by absolute correlation strength
  results <- results %>%
    arrange(P_value)
  
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
  "3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/ppv_membership_correlation_0_21.csv",
  row.names = FALSE
)

write.csv(
  cat_membership_corr,
  "3_data_analysis/6_clustering_modeling/octa_analysis/T2_T0/cat_membership_correlation_0_21.csv",
  row.names = FALSE
)

# -----------------------------------------------------
# 10. Save processed data for further analysis
# -----------------------------------------------------
# Save the filtered datasets
saveRDS(ppv_bloodflow_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/ppv_bloodflow_0_21_analysis.rds")
saveRDS(cat_bloodflow_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/cat_bloodflow_0_21_analysis.rds")
saveRDS(ppv_thickness_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/ppv_thickness_0_21_analysis.rds")
saveRDS(cat_thickness_0_21, "3_data_analysis/6_clustering_modeling/octa_analysis/cat_thickness_0_21_analysis.rds")

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
