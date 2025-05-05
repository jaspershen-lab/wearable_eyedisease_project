library(tidyverse)
library(meta)
library(metafor)
library(gridExtra)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggstatsplot)  # For ggbetweenstats
library(emmeans) 

# 设置工作目录
setwd(get_project_wd())  # 取决于你的项目设置
rm(list = ls())

# 加载数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1m_prediction/daily_data_grouped/combined_all_days_groups.csv")

# 创建输出目录
dir.create("3_data_analysis/5_presurgery_analysis/wearable_prior_after_surgery_compare/surgery_group", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/wearable_prior_after_surgery_compare/surgery_group")

# 加载和处理基线信息（假设baseline_info包含diabetes_history和surgery_type）
# 如果baseline_info不在workspace中，应该先加载它
# baseline_info <- read_csv("path_to_baseline_info.csv")
combined_data <- combined_data %>%
  left_join(baseline_info, by = c("subject_id" = "ID"))

# 检查合并结果
cat("合并后，有dm_status和surgery_type的行数:", 
    sum(!is.na(combined_data$dm_status)), "/", nrow(combined_data), "\n")

# Define time periods based on day_num
combined_data <- combined_data %>%
  mutate(time_period = case_when(
    day_num >= -7 & day_num <= -1 ~ "Pre-surgery",
    day_num >= 1 & day_num <= 7 ~ "Week 1",
    day_num >= 8 & day_num <= 14 ~ "Week 2",
    day_num >= 15 & day_num <= 21 ~ "Week 3",
    day_num >= 22 & day_num <= 30 ~ "Week 4",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(time_period))

# Convert time_period to factor with correct order
combined_data$time_period <- factor(combined_data$time_period, 
                                    levels = c("Pre-surgery", "Week 1", "Week 2", "Week 3", "Week 4"))

# 创建组合分组变量
combined_data <- combined_data %>%
  mutate(
    group_label = case_when(
      dm_status == "No Diabetes" & surgery_type == 0 ~ "No Diabetes Surgery 0",
      dm_status == "Diabetes" & surgery_type == 1 ~ "Diabetes Surgery 1",
      TRUE ~ "Other"
    )
  ) %>%
  filter(group_label != "Other")  # 过滤掉不符合两个主要组的数据

# 检查数据
cat("数据总行数:", nrow(combined_data), "\n")
cat("组别分布:\n")
print(table(combined_data$group_label, useNA = "ifany"))
cat("时间段分布:\n")
print(table(combined_data$time_period, useNA = "ifany"))

# 设置主题
theme_base <- theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  )

# 定义组别颜色
group_colors <- c("#e6550d", "#3182bd")

# 定义分析变量和展示标题
variables <- list(
  list(var = "mean_rhr_1", title = "Resting Heart Rate"),
  list(var = "mean_bo", title = "Blood Oxygen (SpO2)"),
  list(var = "total_sleep", title = "Total Sleep (mins)"),
  list(var = "steps_total", title = "Total Steps")
)

# 设置组值
no_diabetes_surgery0_val <- "No Diabetes Surgery 0"
diabetes_surgery1_val <- "Diabetes Surgery 1"

# 简化后的绘图函数
plot_with_significance <- function(data, variable, title, group_color) {
  # 过滤有效数据
  plot_data <- data %>% filter(!is.na(!!sym(variable)))
  
  if (nrow(plot_data) < 10) {
    cat("警告: 数据不足以创建图表\n")
    return(NULL)
  }
  
  # 创建基础箱型图
  p <- ggplot(plot_data, aes(x = time_period, y = !!sym(variable), fill = time_period)) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    labs(
      title = title,
      x = "Time Period",
      y = title
    ) +
    scale_fill_manual(values = rep(group_color, 5)) +
    theme_base +
    theme(legend.position = "none")
  
  # 使用ggpubr的compare_means来添加显著性标记
  p <- p + 
    stat_compare_means(
      method = "t.test", 
      comparisons = list(
        c("Pre-surgery", "Week 1"),
        c("Pre-surgery", "Week 2"),
        c("Pre-surgery", "Week 3"),
        c("Pre-surgery", "Week 4")
      ),
      label = "p.signif"
    )
  
  return(p)
}

# 执行每组的绘图
for (group_val in c(no_diabetes_surgery0_val, diabetes_surgery1_val)) {
  group_color <- ifelse(group_val == no_diabetes_surgery0_val, "#e6550d", "#3182bd")
  
  cat("\n\n==========================================\n")
  cat(group_val, "的各指标随时间变化\n")
  cat("==========================================\n\n")
  
  # 过滤特定组的数据
  group_data <- combined_data %>% filter(group_label == group_val)
  
  # 为每个变量生成图表
  for (var_info in variables) {
    cat("生成", var_info$title, "的图表...\n")
    
    # 创建图表
    plot_title <- paste(var_info$title, "-", group_val)
    plot <- plot_with_significance(group_data, var_info$var, plot_title, group_color)
    
    if (!is.null(plot)) {
      # 显示图表
      print(plot)
      
      # 保存图表
      ggsave(
        filename = paste0(gsub(" ", "_", tolower(var_info$title)), "_", 
                          tolower(gsub(" ", "_", group_val)), ".pdf"),
        plot = plot,
        width = 10,
        height = 7,
        dpi = 300
      )
      
      # 另外执行统计检验并输出结果 (可选)
      stat_test <- group_data %>% 
        filter(!is.na(!!sym(var_info$var))) %>%
        pairwise_t_test(
          formula = as.formula(paste(var_info$var, "~ time_period")),
          p.adjust.method = "bonferroni"
        ) %>%
        filter(group1 == "Pre-surgery")
      
      cat("\n统计结果 -", var_info$title, ":\n")
      print(stat_test)
      cat("\n")
    }
  }
}

# 方法2：使用compare_means产生统计表格
create_stat_tables <- function() {
  # 创建目录
  dir.create("stat_tables", showWarnings = FALSE)
  
  for (group_val in c(no_diabetes_surgery0_val, diabetes_surgery1_val)) {
    # 过滤特定组的数据
    group_data <- combined_data %>% filter(group_label == group_val)
    
    # 创建统计结果表
    stat_results <- data.frame()
    
    for (var_info in variables) {
      # 跳过缺失数据过多的变量
      if (sum(!is.na(group_data[[var_info$var]])) < 10) {
        next
      }
      
      # 计算统计结果
      temp_stats <- group_data %>%
        group_by(time_period) %>%
        summarise(
          mean = mean(!!sym(var_info$var), na.rm = TRUE),
          sd = sd(!!sym(var_info$var), na.rm = TRUE),
          n = sum(!is.na(!!sym(var_info$var))),
          .groups = "drop"
        )
      
      # 计算Pre-surgery与其他时期的p值
      p_values <- compare_means(
        formula = as.formula(paste(var_info$var, "~ time_period")),
        data = group_data,
        method = "t.test",
        ref.group = "Pre-surgery",
        p.adjust.method = "bonferroni"
      )
      
      # 将p值合并到统计结果中
      temp_stats <- temp_stats %>%
        left_join(
          p_values %>% dplyr::select(time_period = group2, p.adj, p.signif),
          by = "time_period"
        )
      
      # 添加变量信息
      temp_stats$variable <- var_info$title
      
      # 合并到总结果表
      stat_results <- bind_rows(stat_results, temp_stats)
    }
    
    # 保存统计结果
    write_csv(
      stat_results,
      paste0("stat_tables/", tolower(gsub(" ", "_", group_val)), "_stats.csv")
    )
    
    # 打印统计结果
    cat("\n", group_val, "统计结果表:\n")
    print(stat_results)
  }
}

# 执行统计表格创建
create_stat_tables()

# 交互分析部分 - 使用group_label代替之前的group变量
# Load necessary libraries
library(tidyverse)
library(car)
library(emmeans)
library(ggplot2)

# Create the interaction plot with annotations
create_interaction_plot <- function(model_data, var_info, interaction_p, group_colors) {
  
  # Calculate estimated marginal means for the interaction
  interaction_model <- lm(as.formula(paste(var_info$var, "~ group_label * time_period + age + gender + bmi")), 
                          data = model_data)
  
  # Get estimated marginal means
  emm_interaction <- emmeans(interaction_model, specs = ~ group_label | time_period)
  emm_df <- as.data.frame(emm_interaction)
  
  # 已经有了组标签，不需要额外转换
  
  # Create interaction plot with annotations
  interaction_plot <- ggplot(emm_df, aes(x = time_period, y = emmean, group = group_label, color = group_label)) +
    geom_line(size = 1.2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    geom_point(size = 3) +
    labs(
      title = paste(var_info$title, "- Group by Time Interaction"),
      subtitle = paste0("Interaction effect: ", 
                        ifelse(interaction_p < 0.05, "Significant (p = ", "Non-significant (p = "), 
                        round(interaction_p, 4), ")"),
      x = "Time Period",
      y = var_info$title,
      color = "Group"
    ) +
    scale_color_manual(values = group_colors) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  return(interaction_plot)
}

# Function to create ANCOVA plots
create_ancova_plot <- function(results, title, group_color) {
  if (is.null(results)) {
    return(NULL)
  }
  
  # Extract data from results
  predicted_means <- results$predicted_means
  pairwise_comp <- results$pairwise_comp
  
  # Create plot with annotations
  p <- ggplot(predicted_means, aes(x = time_period, y = emmean)) +
    geom_bar(stat = "identity", fill = group_color, alpha = 0.7, width = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    labs(
      title = paste0(title, " (Adjusted for Age, Gender, and BMI)"),
      subtitle = "Estimated Marginal Means ± 95% CI",
      x = "Time Period",
      y = title
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  
  # Add significance markers if available
  if (!is.null(pairwise_comp) && nrow(pairwise_comp) > 0) {
    # Add significance markers for comparisons with Pre-surgery
    comp_with_baseline <- pairwise_comp %>%
      filter(grepl("Pre-surgery", contrast))
    
    if (nrow(comp_with_baseline) > 0) {
      # Calculate positions for significance markers
      y_max <- max(predicted_means$upper.CL) * 1.1
      y_positions <- seq(y_max, y_max * 1.2, length.out = nrow(comp_with_baseline))
      
      # Add markers for significant comparisons
      for (i in 1:nrow(comp_with_baseline)) {
        comp <- comp_with_baseline[i, ]
        if (comp$p.value < 0.05) {
          group2 <- gsub("Pre-surgery - ", "", comp$contrast)
          signif_marker <- ifelse(comp$p.value < 0.001, "***", 
                                  ifelse(comp$p.value < 0.01, "**", "*"))
          
          p <- p + 
            annotate("text", 
                     x = which(levels(predicted_means$time_period) == group2), 
                     y = y_positions[i], 
                     label = signif_marker, 
                     size = 6)
        }
      }
    }
  }
  
  return(p)
}

# 对所有变量创建交互图
for (var_info in variables) {
  # Filter valid data
  model_data <- combined_data %>%
    filter(!is.na(!!sym(var_info$var)) & !is.na(age) & !is.na(gender) & !is.na(bmi))
  
  if (nrow(model_data) < 50) {
    next
  }
  
  # Run interaction model
  formula_str <- paste(var_info$var, "~ group_label * time_period + age + gender + bmi")
  interaction_model <- tryCatch({
    lm(as.formula(formula_str), data = model_data)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(interaction_model)) {
    next
  }
  
  # Get interaction p-value
  anova_results <- car::Anova(interaction_model, type = 3)
  interaction_p <- anova_results["group_label:time_period", "Pr(>F)"]
  
  # Create plot
  plot <- create_interaction_plot(model_data, var_info, interaction_p, group_colors)
  
  # Save plot
  ggsave(
    filename = paste0("ancova_", gsub(" ", "_", tolower(var_info$title)), "_interaction.pdf"),
    plot = plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}
