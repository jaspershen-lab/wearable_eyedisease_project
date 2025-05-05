library(tidyverse)
library(forestplot)
library(meta)
library(metafor)
library(gridExtra)
library(gtsummary)
library(ggpubr) # For statistical tests on plots
library(dplyr)
setwd(get_project_wd())
rm(list = ls())

# First, source the 100_tools.R file to load the color scheme and theme
source("1_code/100_tools.R")

# Rest of the data preparation code remains the same
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")

dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/va_compare", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/va_compare")

# Step 1: Filter pre-surgery data (days -7 to -1)
pre_surgery_data <- combined_data %>%
  filter(as.numeric(day) >= -7 & as.numeric(day) <= -1)  # Only pre-surgery 7 days data

# Step 2: Group by subject ID and calculate pre-surgery means
patient_preop_means <- pre_surgery_data %>%
  group_by(subject_id) %>%
  summarise(
    # Select key mean metrics
    mean_hr = mean(mean_hr, na.rm = TRUE),           # Mean heart rate
    min_hr = mean(min_hr, na.rm = TRUE),  
    max_hr= mean(max_hr, na.rm = TRUE),
    mean_hr = mean(mean_hr, na.rm = TRUE),  
    mean_rhr = mean(mean_rhr_1, na.rm = TRUE),       # Mean resting heart rate
    min_rhr = mean(min_rhr_1, na.rm = TRUE), 
    max_rhr = mean(max_rhr_1, na.rm = TRUE),
    skew_rhr= mean(skew_rhr_1, na.rm=TRUE),
    kurt_rhr= mean(kurt_rhr_1, na.rm=TRUE),
    sd_rhr=mean(sd_rhr_1,na.rm=TRUE),
    mean_rhr_50 = mean(mean_rhr_50, na.rm = TRUE),       # Mean resting heart rate
    min_rhr_50 = mean(min_rhr_50, na.rm = TRUE), 
    max_rhr_50 = mean(max_rhr_50, na.rm = TRUE),
    skew_rhr_50= mean(skew_rhr_50, na.rm=TRUE),
    kurt_rhr_50= mean(kurt_rhr_50, na.rm=TRUE),
    sd_rhr_50=mean(sd_rhr_50,na.rm=TRUE),
    median_rhr = mean(median_rhr_1, na.rm = TRUE), 
    mean_bo = mean(mean_bo, na.rm = TRUE),           # Mean blood oxygen
    min_bo = mean(min_bo, na.rm = TRUE), 
    max_bo = mean(mean_bo, na.rm = TRUE), 
    cv_bo = mean(cv_bo, na.rm = TRUE), 
    total_steps = mean(steps_total, na.rm = TRUE),   # Mean total steps
    max_steps=mean(steps_max, na.rm = TRUE),
    total_sleep = mean(total_sleep, na.rm = TRUE),   # Mean total sleep
    deep_sleep = mean(deep_sleep, na.rm = TRUE), 
    light_sleep = mean(light_sleep, na.rm = TRUE),   # Mean total sleep
    dream_sleep = mean(dream_sleep, na.rm = TRUE), 
    
    # Retain other demographic variables
    dm_2 = dm_2[1],                    # Diabetes status
    pre_vision = pre_vision[1],        # Pre-surgery vision
    age = age[1],                      # Age
    gender = gender[1],                # Gender
    bmi = bmi[1],                      # BMI
    hypertension_2 = hypertension_2[1], # Hypertension status
    season = season[1]  # 获取因子标签
  ) %>%
  # Remove those with missing vision or diabetes status
  filter(!is.na(pre_vision) & !is.na(dm_2))

# Step 3: Handle missing values
# Calculate proportion of missing values
missing_values <- colSums(is.na(patient_preop_means))/nrow(patient_preop_means)
print(missing_values)

# Perform median imputation
data_imputed <- patient_preop_means
for(col in names(data_imputed)[2:(ncol(data_imputed)-6)]) {  # Skip subject_id and demographic variables
  if(sum(is.na(data_imputed[[col]])) > 0) {
    data_imputed[[col]][is.na(data_imputed[[col]])] <- median(data_imputed[[col]], na.rm = TRUE)
  }
}

# Step 4: Split data into diabetic and non-diabetic groups
diabetic_data <- data_imputed %>% filter(dm_2 == 1)
nondiabetic_data <- data_imputed %>% filter(dm_2 == 0)

# Print sample sizes
cat("Sample size of diabetic patients:", nrow(diabetic_data), "subjects\n")
cat("Sample size of non-diabetic patients:", nrow(nondiabetic_data), "subjects\n")



# Step 5: Compare wearable parameters between diabetic and non-diabetic groups
# Create comparison table
comparison_table <- data_imputed %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "DR", "Cataract"),
    gender = factor(gender, levels = c(1, 0), labels = c("Male", "Female"))
  ) %>%
  dplyr::select( mean_rhr, mean_bo, total_steps, total_sleep, 
                 pre_vision, age, gender, bmi, dm_status) %>%
  tbl_summary(
    by = dm_status,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(
      mean_rhr ~ "Mean Resting Heart Rate",
      mean_bo ~ "Mean Blood Oxygen",
      total_steps ~ "Mean Total Steps",
      total_sleep ~ "Mean Total Sleep Time",
      pre_vision ~ "Pre-surgery Vision",
      age ~ "Age",
      gender ~ "Gender",
      bmi ~ "BMI"
    )
  ) %>%
  add_p() %>%
  add_n() %>%
  modify_header(label = "**Parameter**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Comparison by DR Status**") %>%
  modify_footnote(
    all_stat_cols() ~ "Mean (SD) or n (%)",
    p.value ~ "Wilcoxon rank-sum test; Chi-square test"
  )

# Print comparison table
print(comparison_table)



# Step 6: Create box plots comparing wearable parameters between groups using the 100_tools colors
# Prepare data in long format for boxplots
# Assuming long_data is created somewhere in the code - recreate it here if missing
long_data <- data_imputed %>%
  pivot_longer(
    cols = c(mean_rhr, mean_bo, total_steps, total_sleep),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "DR", "Cataract"),
    parameter_label = case_when(
      parameter == "mean_rhr" ~ "Mean Resting Heart Rate",
      parameter == "mean_bo" ~ "Mean Blood Oxygen",
      parameter == "total_steps" ~ "Mean Total Steps",
      parameter == "total_sleep" ~ "Mean Total Sleep"
    )
  )

# Use ggplot2 instead of ggpubr for better control over colors
comparison_boxplots <- ggplot(long_data, aes(x = dm_status, y = value, fill = dm_status)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~ parameter_label, scales = "free_y") +
  # Use custom colors directly from 100_tools.R
  scale_fill_manual(values = dr_colors) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",
                     label.y.position = "top") +
  labs(
    title = "Comparison of Pre-surgery Parameters Between Diabetic and Non-diabetic Patients",
    subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
    x = "",
    y = "",
    fill = "Group"
  ) +
  theme_custom

# Display the plot
print(comparison_boxplots)

# Save the plot
ggsave("diabetic_nondiabetic_comparison_with_pvalues.pdf", 
       comparison_boxplots, width = 14, height = 10)




# Step 7: Correlation analysis for both groups using Spearman's rank correlation
# Function to calculate Spearman correlation with vision, adjusted for covariates
calculate_correlation_vision <- function(data, var_name) {
  # Basic Spearman correlation (unadjusted)
  spearman_cor <- cor(data[[var_name]], data$pre_vision, 
                      method = "spearman", use = "complete.obs")
  
  # Calculate N for reporting
  n <- sum(!is.na(data[[var_name]]) & !is.na(data$pre_vision))
  
  # Calculate p-value for Spearman correlation
  spearman_test <- cor.test(data[[var_name]], data$pre_vision, 
                            method = "spearman", exact = FALSE)
  p_value <- spearman_test$p.value
  
  # Adjusted correlation using partial correlation with covariates
  # First, create a complete cases dataset for the analysis
  vars_for_model <- c(var_name, "pre_vision", "age", "gender", "bmi")
  
  # Include diabetes status as a covariate if we're using the full dataset
  # and the column exists
  if("diabetes" %in% colnames(data) && data$Group[1] == "All Participants") {
    vars_for_model <- c(vars_for_model, "diabetes")
  }
  
  complete_data <- data[complete.cases(data[, vars_for_model]), vars_for_model]
  
  if(nrow(complete_data) > 5) { # Ensure we have enough data for the adjustment
    # Using the ppcor package for partial correlation
    # If not already installed, uncomment the line below
    # install.packages("ppcor")
    library(ppcor)
    
    # Calculate partial Spearman correlation controlling for age, gender, and bmi
    partial_result <- pcor(complete_data, method = "spearman")
    
    # Find the row and column corresponding to var_name and pre_vision
    var_index <- which(colnames(complete_data) == var_name)
    vision_index <- which(colnames(complete_data) == "pre_vision")
    
    # Extract the adjusted correlation coefficient and p-value
    adjusted_correlation <- partial_result$estimate[var_index, vision_index]
    adjusted_p_value <- partial_result$p.value[var_index, vision_index]
    adjusted_n <- nrow(complete_data)
    
    # Return both unadjusted and adjusted correlations
    return(data.frame(
      Variable = var_name,
      Unadjusted_Correlation = spearman_cor,
      Unadjusted_P_value = p_value,
      Unadjusted_N = n,
      Adjusted_Correlation = adjusted_correlation,
      Adjusted_P_value = adjusted_p_value,
      Adjusted_N = adjusted_n
    ))
  } else {
    # If not enough data for adjustment, return only unadjusted
    return(data.frame(
      Variable = var_name,
      Unadjusted_Correlation = spearman_cor,
      Unadjusted_P_value = p_value,
      Unadjusted_N = n,
      Adjusted_Correlation = NA,
      Adjusted_P_value = NA,
      Adjusted_N = NA
    ))
  }
}

# Calculate correlations for both groups
wearable_vars <- c("mean_rhr", "max_rhr","min_rhr","sd_rhr","max_bo","min_bo","mean_bo", "total_steps", "total_sleep","deep_sleeps")

# Diabetic group correlations
diabetic_correlations <- do.call(rbind, lapply(wearable_vars, function(var) {
  calculate_correlation_vision(diabetic_data, var)
}))
diabetic_correlations$Group <- "DR"

# Non-diabetic group correlations
nondiabetic_correlations <- do.call(rbind, lapply(wearable_vars, function(var) {
  calculate_correlation_vision(nondiabetic_data, var)
}))
nondiabetic_correlations$Group <- "Cataract"

# Using the provided data_imputed as combined data
combined_data_correlations <- do.call(rbind, lapply(wearable_vars, function(var) {
  calculate_correlation_vision(data_imputed, var)
}))
combined_data_correlations$Group <- "All Participants"

# Combine results from all three analyses
combined_correlations <- rbind(diabetic_correlations, nondiabetic_correlations, combined_data_correlations)

# Add human-readable variable names
combined_correlations <- combined_correlations %>%
  mutate(
    Variable_Label = case_when(
      # Variable == "mean_hr" ~ "Mean Heart Rate",
      Variable == "mean_rhr" ~ "Mean Resting Heart Rate",
      Variable == "mean_bo" ~ "Mean Blood Oxygen", 
      Variable == "total_steps" ~ "Mean Total Steps",
      Variable == "total_sleep" ~ "Mean Total Sleep"
    )
  )

# Print correlation results
print(combined_correlations)

# Visualize the results
library(ggplot2)


# 多变量回归模型（合并数据）
model_combined <- lm(pre_vision ~ mean_rhr + mean_bo + total_steps + total_sleep + 
                       age + gender + bmi + dm_2, data = data_imputed)
summary(model_combined)

# 糖尿病组多变量回归模型
model_diabetic <- lm(pre_vision ~ mean_rhr + mean_bo +total_steps + total_sleep + deep_sleep+
                       age + gender + bmi, data = diabetic_data)
summary(model_diabetic)

# 非糖尿病组多变量回归模型
model_nondiabetic <- lm(pre_vision ~ mean_rhr + mean_bo + total_steps + total_sleep + 
                          age + gender + bmi, data = nondiabetic_data)
summary(model_nondiabetic)


#-------------------------
# 1. 散点图可视化
#-------------------------
# 创建各变量与视力的散点图并添加平滑曲线
create_scatter_with_smooth <- function(data, x_var, title) {
  ggplot(data, aes_string(x = x_var, y = "pre_vision")) +
    geom_point(aes(color = factor(dm_2)), alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +  # 添加非参数平滑曲线
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
    labs(title = paste("Vision vs", title),
         x = title, y = "Vision Score",
         color = "DR Status") +
    theme_bw()
}



##########logitic_reg
# 将mean_rhr分为高心率和低心率组
# 使用中位数作为分界点
data_imputed$high_rhr <- ifelse(data_imputed$mean_rhr > 75, 1, 0)

# 查看高低心率组的分布
table(data_imputed$high_rhr)

# 进行逻辑回归：pre_vision预测高/低心率
model_logistic <- glm(high_rhr ~ pre_vision, 
                      family = binomial(link = "logit"), 
                      data = data_imputed)

# 查看模型摘要
summary(model_logistic)

# 计算效应量 - 优势比(Odds Ratio)及其置信区间
OR <- exp(coef(model_logistic))
CI <- exp(confint(model_logistic))
OR_table <- data.frame(
  OR = OR,
  Lower_CI = CI[,1],
  Upper_CI = CI[,2]
)
print(OR_table)

# 绘制视力得分与高心率概率的关系图
library(ggplot2)

# 预测不同视力得分下高心率的概率
vision_range <- seq(min(data_imputed$pre_vision, na.rm=TRUE), 
                    max(data_imputed$pre_vision, na.rm=TRUE), 
                    length.out=100)
pred_data <- data.frame(pre_vision = vision_range)
pred_data$high_rhr_prob <- predict(model_logistic, newdata=pred_data, type="response")

# 绘制概率曲线
p <- ggplot(pred_data, aes(x=pre_vision, y=high_rhr_prob)) +
  geom_line(color="blue", size=1) +
  geom_point(data=data_imputed, aes(x=pre_vision, y=high_rhr), alpha=0.5) +
  labs(title="预测概率：视力得分预测高心率",
       x="视力得分(pre_vision)",
       y="高心率概率") +
  theme_bw() +
  geom_hline(yintercept=0.5, linetype="dashed", color="red") +
  annotate("text", x=min(data_imputed$pre_vision, na.rm=TRUE), y=0.55, 
           label="临界值", color="red", hjust=0)

# 显示图表
print(p)

# 添加糖尿病状态作为协变量的模型
model_logistic_adj <- glm(high_rhr ~ pre_vision + dm_2 + age + gender + bmi, 
                          family = binomial(link = "logit"), 
                          data = data_imputed)

# 查看调整后的模型摘要
summary(model_logistic_adj)

# 计算调整后的效应量
OR_adj <- exp(coef(model_logistic_adj))
CI_adj <- exp(confint(model_logistic_adj))
OR_table_adj <- data.frame(
  OR = OR_adj,
  Lower_CI = CI_adj[,1],
  Upper_CI = CI_adj[,2]
)
print(OR_table_adj)

#---------------------------------
# 按糖尿病状态进行分层分析
#---------------------------------

# 糖尿病组分析
dm_data <- data_imputed[data_imputed$dm_2 == 1, ]
dm_data$high_rhr <- ifelse(dm_data$mean_rhr > median(dm_data$mean_rhr), 1, 0)
model_dm <- glm(high_rhr ~ pre_vision, family = binomial(link = "logit"), data = dm_data)
summary(model_dm)

# 非糖尿病组分析
non_dm_data <- data_imputed[data_imputed$dm_2 == 0, ]
non_dm_data$high_rhr <- ifelse(non_dm_data$mean_rhr > median(non_dm_data$mean_rhr), 1, 0)
model_non_dm <- glm(high_rhr ~ pre_vision, family = binomial(link = "logit"), data = non_dm_data)
summary(model_non_dm)

# 交互项分析 - 检验糖尿病状态是否调节视力与心率的关系
model_interaction <- glm(high_rhr ~ pre_vision * dm_2, 
                         family = binomial(link = "logit"), 
                         data = data_imputed)
summary(model_interaction)




######All viriables& va heatmap
# Step 7: Correlation analysis for all variables
# First, identify which variables to include in the correlation analysis
# We'll exclude subject_id, dm_2 (will be used for stratification), and categorical variables

# Define variables to include in correlation analysis
# Exclude categorical and identifier variables
exclude_vars <- c("subject_id", "season", "gender", "hypertension_2")
demo_vars <- c("age", "bmi", "pre_vision")

# Get all numeric variable names from the dataset
wearable_vars <- names(data_imputed)[sapply(data_imputed, is.numeric)]
wearable_vars <- setdiff(wearable_vars, c(exclude_vars, demo_vars, "dm_2"))

# Combine all variables for analysis
all_numeric_vars <- c(wearable_vars, demo_vars)

# Function to calculate correlation matrix
calculate_correlation_matrix <- function(data, vars) {
  # Create a matrix to store correlations
  cor_matrix <- matrix(NA, nrow = length(vars), ncol = length(vars))
  rownames(cor_matrix) <- vars
  colnames(cor_matrix) <- vars
  
  # Calculate Spearman correlations
  for (i in 1:length(vars)) {
    for (j in 1:length(vars)) {
      cor_matrix[i, j] <- cor(data[[vars[i]]], data[[vars[j]]], 
                              method = "spearman", use = "pairwise.complete.obs")
    }
  }
  
  return(cor_matrix)
}

# Calculate p-value matrix
calculate_pvalue_matrix <- function(data, vars) {
  # Create a matrix to store p-values
  p_matrix <- matrix(NA, nrow = length(vars), ncol = length(vars))
  rownames(p_matrix) <- vars
  colnames(p_matrix) <- vars
  
  # Calculate p-values for Spearman correlations
  for (i in 1:length(vars)) {
    for (j in 1:length(vars)) {
      test_result <- cor.test(data[[vars[i]]], data[[vars[j]]], 
                              method = "spearman", exact = FALSE,
                              use = "pairwise.complete.obs")
      p_matrix[i, j] <- test_result$p.value
    }
  }
  
  return(p_matrix)
}

# Calculate correlation matrices for all participants, diabetic, and non-diabetic groups
all_cor_matrix <- calculate_correlation_matrix(data_imputed, all_numeric_vars)
all_p_matrix <- calculate_pvalue_matrix(data_imputed, all_numeric_vars)

diabetic_cor_matrix <- calculate_correlation_matrix(diabetic_data, all_numeric_vars)
diabetic_p_matrix <- calculate_pvalue_matrix(diabetic_data, all_numeric_vars)

nondiabetic_cor_matrix <- calculate_correlation_matrix(nondiabetic_data, all_numeric_vars)
nondiabetic_p_matrix <- calculate_pvalue_matrix(nondiabetic_data, all_numeric_vars)

# Create nicer variable labels for plotting
variable_labels <- c(
  "mean_hr" = "Mean HR",
  "min_hr" = "Min HR",
  "max_hr" = "Max HR",
  "mean_rhr" = "Mean RHR",
  "min_rhr" = "Min RHR",
  "max_rhr" = "Max RHR",
  "skew_rhr" = "Skewness RHR",
  "kurt_rhr" = "Kurtosis RHR",
  "sd_rhr" = "SD RHR",
  "mean_rhr_50" = "Mean RHR 50",
  "min_rhr_50" = "Min RHR 50",
  "max_rhr_50" = "Max RHR 50",
  "skew_rhr_50" = "Skewness RHR 50",
  "kurt_rhr_50" = "Kurtosis RHR 50",
  "sd_rhr_50" = "SD RHR 50",
  "median_rhr" = "Median RHR",
  "mean_bo" = "Mean BO",
  "min_bo" = "Min BO",
  "max_bo" = "Max BO",
  "cv_bo" = "CV BO",
  "total_steps" = "Total Steps",
  "max_steps" = "Max Steps",
  "total_sleep" = "Total Sleep",
  "deep_sleep" = "Deep Sleep",
  "light_sleep" = "Light Sleep",
  "dream_sleep" = "REM Sleep",
  "age" = "Age",
  "bmi" = "BMI",
  "pre_vision" = "Pre-surgery Vision"
)

# Convert correlation matrices to data frames for plotting
library(reshape2)

# Function to convert matrix to long format for ggplot
matrix_to_long <- function(cor_matrix, p_matrix = NULL, sig_level = 0.05) {
  # Convert to data frame
  cor_df <- as.data.frame(cor_matrix)
  cor_df$var1 <- rownames(cor_df)
  
  # Melt to long format
  long_df <- melt(cor_df, id.vars = "var1", variable.name = "var2", value.name = "correlation")
  
  # Add significance indicator if p-value matrix is provided
  if (!is.null(p_matrix)) {
    p_df <- as.data.frame(p_matrix)
    p_df$var1 <- rownames(p_df)
    p_long <- melt(p_df, id.vars = "var1", variable.name = "var2", value.name = "p_value")
    
    # Merge correlation and p-value data frames
    long_df <- merge(long_df, p_long, by = c("var1", "var2"))
    
    # Add significance indicator
    long_df$significant <- long_df$p_value < sig_level
  }
  
  return(long_df)
}

# Convert matrices to long format
all_cor_long <- matrix_to_long(all_cor_matrix, all_p_matrix)
diabetic_cor_long <- matrix_to_long(diabetic_cor_matrix, diabetic_p_matrix)
nondiabetic_cor_long <- matrix_to_long(nondiabetic_cor_matrix, nondiabetic_p_matrix)

# Add group labels
all_cor_long$group <- "All Participants"
diabetic_cor_long$group <- "Diabetic Patients"
nondiabetic_cor_long$group <- "Non-diabetic Patients"

# Combine all data for plotting
combined_cor_long <- rbind(all_cor_long, diabetic_cor_long, nondiabetic_cor_long)

# We're not using combined_cor_long anymore in the initial heatmaps
# This section will be handled differently for the vision-specific plots

# Create correlation heatmaps
library(ggplot2)

# First, ensure the var1_label and var2_label columns exist
all_cor_long$var1_label <- factor(variable_labels[as.character(all_cor_long$var1)])
all_cor_long$var2_label <- factor(variable_labels[as.character(all_cor_long$var2)])

diabetic_cor_long$var1_label <- factor(variable_labels[as.character(diabetic_cor_long$var1)])
diabetic_cor_long$var2_label <- factor(variable_labels[as.character(diabetic_cor_long$var2)])

nondiabetic_cor_long$var1_label <- factor(variable_labels[as.character(nondiabetic_cor_long$var1)])
nondiabetic_cor_long$var2_label <- factor(variable_labels[as.character(nondiabetic_cor_long$var2)])

# Now creating combined_cor_long is not needed here; we'll work with each dataset separately
# combined_cor_long <- rbind(all_cor_long, diabetic_cor_long, nondiabetic_cor_long)

# Function to create a heatmap
create_heatmap <- function(data, title) {
  ggplot(data, aes(x = var2_label, y = var1_label, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1), name = "Correlation") +
    geom_text(aes(label = ifelse(significant, sprintf("%.2f*", correlation), 
                                 sprintf("%.2f", correlation))),
              size = 2.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = title, x = "", y = "")
}

# Full correlation heatmaps
all_heatmap <- create_heatmap(all_cor_long, "Correlation Matrix - All Participants")
diabetic_heatmap <- create_heatmap(diabetic_cor_long, "Correlation Matrix - Diabetic Patients")
nondiabetic_heatmap <- create_heatmap(nondiabetic_cor_long, "Correlation Matrix - Non-diabetic Patients")
all_heatmap
diabetic_heatmap
nondiabetic_heatmap
# Save heatmaps
ggsave("all_correlation_heatmap.pdf", all_heatmap, width = 12, height = 10)
ggsave("diabetic_correlation_heatmap.pdf", diabetic_heatmap, width = 12, height = 10)
ggsave("nondiabetic_correlation_heatmap.pdf", nondiabetic_heatmap, width = 12, height = 10)

# Create focused correlation plots with pre_vision only
# First, prepare the combined data for the vision correlations
all_cor_long$group <- "All Participants"
diabetic_cor_long$group <- "Diabetic Patients"
nondiabetic_cor_long$group <- "Non-diabetic Patients"

# Now combine them
combined_cor_long <- rbind(all_cor_long, diabetic_cor_long, nondiabetic_cor_long)

# Ensure label columns exist with correct values
combined_cor_long$var1_label <- factor(variable_labels[as.character(combined_cor_long$var1)])
combined_cor_long$var2_label <- factor(variable_labels[as.character(combined_cor_long$var2)])

# Extract correlations with pre_vision
vision_cors <- combined_cor_long %>%
  filter(var2 == "pre_vision" & var1 != "pre_vision") %>%
  arrange(group, desc(abs(correlation)))

# Create a bar plot of correlations with pre_vision
vision_plot <- ggplot(vision_cors, aes(x = reorder(var1_label, correlation), y = correlation, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = ifelse(significant, sprintf("%.2f*", correlation), 
                               sprintf("%.2f", correlation))),
            position = position_dodge(width = 0.9), hjust = -0.1, size = 3, angle = 90) +
  coord_flip() +
  scale_fill_manual(values = dr_colors) +
  theme_minimal() +
  labs(title = "Correlation of Variables with Pre-surgery Vision",
       subtitle = "* p < 0.05",
       x = "",
       y = "Spearman Correlation Coefficient",
       fill = "Group") +
  theme(legend.position = "bottom")
vision_plot
# Save the vision correlation plot
ggsave("vision_correlation_plot.pdf", vision_plot, width = 10, height = 12)

# Create a table of correlations with pre_vision
vision_table <- reshape(vision_cors, 
                        idvar = "var1_label",
                        timevar = "group",
                        direction = "wide")

# Rename columns
names(vision_table) <- gsub("correlation.", "", names(vision_table))
names(vision_table) <- gsub("p_value.", "p_", names(vision_table))
names(vision_table) <- gsub("significant.", "sig_", names(vision_table))

# Select and order columns
vision_table <- vision_table[, c("var1_label", 
                                 "All Participants", "p_All Participants", "sig_All Participants",
                                 "Diabetic Patients", "p_Diabetic Patients", "sig_Diabetic Patients",
                                 "Non-diabetic Patients", "p_Non-diabetic Patients", "sig_Non-diabetic Patients")]

# Sort by absolute correlation value in the All Participants group
vision_table <- vision_table[order(-abs(vision_table$`All Participants`)), ]

# Save correlation table
write.csv(vision_table, "vision_correlation_table.csv", row.names = FALSE)

# Create focused heatmaps showing only correlations with pre_vision 
# and selected key variables
key_vars <- c("mean_rhr", "min_rhr", "max_rhr", "sd_rhr", 
              "mean_bo", "min_bo", "max_bo", 
              "total_steps", "max_steps", 
              "total_sleep", "deep_sleep", "light_sleep", "dream_sleep",
              "age", "bmi", "pre_vision")

# Filter correlation data for key variables
key_all_cor_long <- all_cor_long %>% 
  filter(var1 %in% key_vars & var2 %in% key_vars)
key_all_cor_long$var1_label <- factor(variable_labels[as.character(key_all_cor_long$var1)])
key_all_cor_long$var2_label <- factor(variable_labels[as.character(key_all_cor_long$var2)])

key_diabetic_cor_long <- diabetic_cor_long %>% 
  filter(var1 %in% key_vars & var2 %in% key_vars)
key_diabetic_cor_long$var1_label <- factor(variable_labels[as.character(key_diabetic_cor_long$var1)])
key_diabetic_cor_long$var2_label <- factor(variable_labels[as.character(key_diabetic_cor_long$var2)])

key_nondiabetic_cor_long <- nondiabetic_cor_long %>% 
  filter(var1 %in% key_vars & var2 %in% key_vars)
key_nondiabetic_cor_long$var1_label <- factor(variable_labels[as.character(key_nondiabetic_cor_long$var1)])
key_nondiabetic_cor_long$var2_label <- factor(variable_labels[as.character(key_nondiabetic_cor_long$var2)])

# Create focused heatmaps
key_all_heatmap <- create_heatmap(key_all_cor_long, "Key Variables Correlation - All Participants")
key_diabetic_heatmap <- create_heatmap(key_diabetic_cor_long, "Key Variables Correlation - Diabetic Patients")
key_nondiabetic_heatmap <- create_heatmap(key_nondiabetic_cor_long, "Key Variables Correlation - Non-diabetic Patients")

# Save focused heatmaps
ggsave("key_all_correlation_heatmap.pdf", key_all_heatmap, width = 10, height = 8)
ggsave("key_diabetic_correlation_heatmap.pdf", key_diabetic_heatmap, width = 10, height = 8)
ggsave("key_nondiabetic_correlation_heatmap.pdf", key_nondiabetic_heatmap, width = 10, height = 8)

# Print summary of strongest correlations with pre_vision
top_correlations <- vision_cors %>%
  group_by(group) %>%
  top_n(5, abs(correlation)) %>%
  arrange(group, desc(abs(correlation)))

print("Top 5 strongest correlations with pre_vision by group:")
print(top_correlations)
