library(tidyverse)
library(forestplot)
library(meta)
library(metafor)
library(gridExtra)
library(gtsummary)
library(ggpubr) # For statistical tests on plots
library(pheatmap)
library(RColorBrewer)
library(car)        # For VIF and partial residual plots
library(corrplot)   # For correlation visualization
library(dplyr)
setwd(get_project_wd())
rm(list = ls())

# First, source the 100_tools.R file to load the color scheme and theme
source("1_code/100_tools.R")

# Load data
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Create output directory
dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis")

# Process OCTA data as in the original script
# Process OCTA blood flow data for baseline (T0)
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# Function to process blood flow data for each patient at baseline (T0)
process_patient_bloodflow <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # Use left eye data for left eye surgery, right eye data for right eye and both eyes surgery
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # Process each time point
  for(suffix in time_points) {
    # Select columns for current time point
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # Select data and rename columns
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Process each patient individually for baseline (T0)
patient_list_bloodflow <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data_bloodflow <- purrr::map(patient_list_bloodflow, process_patient_bloodflow)

# Combine results
octa_bloodflow_features <- bind_rows(processed_data_bloodflow)

# Create blood flow variables subset for baseline (T0)
bloodflow_var_T0 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),  # Select all columns ending with 0_6_T0
    -matches("PA_OuterRetina_0_6_T0"),  # Exclude these columns
    -matches("PA_PED_0_6_T0")
  )

# Process OCTA thickness data for baseline (T0)
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# Function to process thickness data for each patient at baseline (T0)
process_patient_thickness <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # Use left eye data for left eye surgery, right eye data for right eye and both eyes surgery
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # Process each time point
  for(suffix in time_points) {
    # Select columns for current time point
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # Select data and rename columns
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Process each patient individually for baseline (T0)
patient_list_thickness <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data_thickness <- purrr::map(patient_list_thickness, process_patient_thickness)

# Combine results
octa_thickness_features <- bind_rows(processed_data_thickness)

# Create thickness variables subset for baseline (T0)
thickness_var_T0 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),
    matches("Thickness_PED_0_6_T0")
  )

# Define OCTA variables for analysis
thickness_vars <- colnames(thickness_var_T0)[-1]  # Removing ID column
bloodflow_vars <- colnames(bloodflow_var_T0)[-1]  # Removing ID column

# Step 1: Filter pre-surgery data (days -7 to -1)
pre_surgery_data <- combined_data %>%
  filter(as.numeric(day) >= -7 & as.numeric(day) <= -1)  # Only pre-surgery 7 days data

# Step 2: Group by subject ID and calculate pre-surgery means
patient_preop_means <- pre_surgery_data %>%
  group_by(subject_id) %>%
  summarise(
    # Select key mean metrics
    mean_hr = mean(mean_hr, na.rm = TRUE),           # Mean heart rate
    mean_rhr = mean(mean_rhr_1, na.rm = TRUE),       # Mean resting heart rate
    mean_bo = mean(mean_bo, na.rm = TRUE),           # Mean blood oxygen
    total_steps = mean(steps_total, na.rm = TRUE),   # Mean total steps
    total_sleep = mean(total_sleep, na.rm = TRUE),   # Mean total sleep
    
    # Retain other demographic variables
    dm_2 = dplyr::first(dm_2),                    # Diabetes status
    age = dplyr::first(age),                      # Age
    gender = dplyr::first(gender),                # Gender
    bmi = dplyr::first(bmi),                      # BMI
    hypertension_2 = dplyr::first(hypertension_2) # Hypertension status
  )

# Step 3: Handle missing values
# Calculate proportion of missing values
missing_values <- colSums(is.na(patient_preop_means))/nrow(patient_preop_means)
print(missing_values)

# Perform median imputation
data_imputed <- patient_preop_means
for(col in names(data_imputed)[2:(ncol(data_imputed)-5)]) {  # Skip subject_id and demographic variables
  if(sum(is.na(data_imputed[[col]])) > 0) {
    data_imputed[[col]][is.na(data_imputed[[col]])] <- median(data_imputed[[col]], na.rm = TRUE)
  }
}

# Step 4: Merge wearable data with OCTA data
wearable_thickness_data <- data_imputed %>%
  left_join(thickness_var_T0, by = c("subject_id" = "ID"))

wearable_bloodflow_data <- data_imputed %>%
  left_join(bloodflow_var_T0, by = c("subject_id" = "ID"))

# Step 5: Split data into diabetic and non-diabetic groups
wearable_thickness_dm <- wearable_thickness_data %>% filter(dm_2 == 1)
wearable_thickness_nodm <- wearable_thickness_data %>% filter(dm_2 == 0)

wearable_bloodflow_dm <- wearable_bloodflow_data %>% filter(dm_2 == 1)
wearable_bloodflow_nodm <- wearable_bloodflow_data %>% filter(dm_2 == 0)

# Print sample sizes
cat("Sample size of diabetic patients (thickness):", nrow(wearable_thickness_dm), "subjects\n")
cat("Sample size of non-diabetic patients (thickness):", nrow(wearable_thickness_nodm), "subjects\n")
cat("Sample size of diabetic patients (blood flow):", nrow(wearable_bloodflow_dm), "subjects\n")
cat("Sample size of non-diabetic patients (blood flow):", nrow(wearable_bloodflow_nodm), "subjects\n")

# Step 6: Compare wearable parameters between diabetic and non-diabetic groups
# Create comparison table
comparison_table <- data_imputed %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "Diabetes", "No Diabetes"),
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female"))
  ) %>%
  dplyr::select(mean_hr, mean_rhr, mean_bo, total_steps, total_sleep, 
                age, gender, bmi, dm_status) %>%
  tbl_summary(
    by = dm_status,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(
      mean_hr ~ "Mean Heart Rate",
      mean_rhr ~ "Mean Resting Heart Rate",
      mean_bo ~ "Mean Blood Oxygen",
      total_steps ~ "Mean Total Steps",
      total_sleep ~ "Mean Total Sleep Time",
      age ~ "Age",
      gender ~ "Gender",
      bmi ~ "BMI"
    )
  ) %>%
  add_p() %>%
  add_n() %>%
  modify_header(label = "**Parameter**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Comparison by Diabetes Status**") %>%
  modify_footnote(
    all_stat_cols() ~ "Mean (SD) or n (%)",
    p.value ~ "Wilcoxon rank-sum test; Chi-square test"
  )

# Print comparison table
print(comparison_table)

# Save the comparison table
if(requireNamespace("gt", quietly = TRUE)) {
  comparison_table %>%
    as_gt() %>%
    gt::gtsave("diabetic_nondiabetic_comparison_table.html")
}

# Step 7: Create box plots comparing wearable parameters between groups 
# Prepare data in long format for boxplots
long_data <- data_imputed %>%
  pivot_longer(
    cols = c(mean_hr, mean_rhr, mean_bo, total_steps, total_sleep),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "Diabetes", "No Diabetes"),
    parameter_label = case_when(
      parameter == "mean_hr" ~ "Mean Heart Rate",
      parameter == "mean_rhr" ~ "Mean Resting Heart Rate",
      parameter == "mean_bo" ~ "Mean Blood Oxygen",
      parameter == "total_steps" ~ "Mean Total Steps",
      parameter == "total_sleep" ~ "Mean Total Sleep"
    )
  )

# Use ggplot2 for boxplots with diabetes_colors from 100_tools.R
comparison_boxplots <- ggplot(long_data, aes(x = dm_status, y = value, fill = dm_status)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~ parameter_label, scales = "free_y") +
  scale_fill_manual(values = diabetes_colors) +
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

# Define wearable variables for analysis
wearable_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")
covariate_vars <- c("age", "gender", "bmi")

# Step 8: Correlation analysis between wearable and OCTA variables

# Create a function to calculate correlation for one OCTA variable
calculate_correlation_octa <- function(data, octa_var, wearable_var) {
  # Remove rows with NA values
  complete_data <- data %>%
    dplyr::select(all_of(c(octa_var, wearable_var))) %>%
    na.omit()
  
  if(nrow(complete_data) < 5) {  # Require at least 5 observations
    return(data.frame(
      OCTA_Variable = octa_var,
      Wearable_Variable = wearable_var,
      Correlation = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_value = NA,
      N = nrow(complete_data)
    ))
  }
  
  # Pearson correlation coefficient
  correlation <- cor(complete_data[[octa_var]], complete_data[[wearable_var]], 
                     use = "complete.obs")
  
  # Fisher z transformation for confidence interval
  n <- nrow(complete_data)
  z <- 0.5 * log((1 + correlation) / (1 - correlation))
  se <- 1 / sqrt(n - 3)
  ci_lower <- tanh(z - 1.96 * se)
  ci_upper <- tanh(z + 1.96 * se)
  
  # Calculate p-value
  t_stat <- correlation * sqrt((n - 2) / (1 - correlation^2))
  p_value <- 2 * pt(-abs(t_stat), df = n - 2)
  
  return(data.frame(
    OCTA_Variable = octa_var,
    Wearable_Variable = wearable_var,
    Correlation = correlation,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    P_value = p_value,
    N = n
  ))
}

# Function to calculate correlations for a group of OCTA variables
calculate_all_correlations <- function(data, octa_vars, wearable_vars, group_name) {
  all_correlations <- data.frame()
  
  for(octa_var in octa_vars) {
    for(wearable_var in wearable_vars) {
      correlation_result <- calculate_correlation_octa(data, octa_var, wearable_var)
      correlation_result$Group <- group_name
      all_correlations <- rbind(all_correlations, correlation_result)
    }
  }
  
  # Add FDR correction for multiple testing
  all_correlations$P_adjusted <- p.adjust(all_correlations$P_value, method = "BH")
  
  return(all_correlations)
}

# Calculate correlations for thickness variables
thickness_dm_correlations <- calculate_all_correlations(
  wearable_thickness_dm, thickness_vars, wearable_vars, "Diabetes")

thickness_nodm_correlations <- calculate_all_correlations(
  wearable_thickness_nodm, thickness_vars, wearable_vars, "No Diabetes")

# Calculate correlations for blood flow variables
bloodflow_dm_correlations <- calculate_all_correlations(
  wearable_bloodflow_dm, bloodflow_vars, wearable_vars, "Diabetes")

bloodflow_nodm_correlations <- calculate_all_correlations(
  wearable_bloodflow_nodm, bloodflow_vars, wearable_vars, "No Diabetes")

# Combine all correlations
all_thickness_correlations <- rbind(thickness_dm_correlations, thickness_nodm_correlations)
all_bloodflow_correlations <- rbind(bloodflow_dm_correlations, bloodflow_nodm_correlations)

# Save the correlation results
write.csv(all_thickness_correlations, "thickness_correlations.csv", row.names = FALSE)
write.csv(all_bloodflow_correlations, "bloodflow_correlations.csv", row.names = FALSE)

# Print significant correlations
cat("\nSignificant thickness correlations (p < 0.05):\n")
print(all_thickness_correlations %>% 
        filter(P_value < 0.05) %>% 
        arrange(P_value))

cat("\nSignificant blood flow correlations (p < 0.05):\n")
print(all_bloodflow_correlations %>% 
        filter(P_value < 0.05) %>% 
        arrange(P_value))

# Step 9: Create correlation heatmaps

# Function to create a correlation heatmap
create_correlation_heatmap <- function(correlations, title, filename) {
  # Reshape data for heatmap
  correlation_matrix <- correlations %>%
    filter(!is.na(Correlation)) %>%
    pivot_wider(
      id_cols = OCTA_Variable,
      names_from = c(Wearable_Variable, Group),
      values_from = Correlation
    )
  
  # Extract the matrix and row names
  octa_vars <- correlation_matrix$OCTA_Variable
  correlation_matrix <- as.matrix(correlation_matrix[,-1])
  rownames(correlation_matrix) <- octa_vars
  
  # Create a p-value matrix for significance markers
  p_value_matrix <- correlations %>%
    filter(!is.na(P_value)) %>%
    pivot_wider(
      id_cols = OCTA_Variable,
      names_from = c(Wearable_Variable, Group),
      values_from = P_value
    )
  
  p_value_matrix <- as.matrix(p_value_matrix[,-1])
  rownames(p_value_matrix) <- octa_vars
  
  # Create a matrix for displaying text
  text_matrix <- matrix("", nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
  rownames(text_matrix) <- rownames(correlation_matrix)
  colnames(text_matrix) <- colnames(correlation_matrix)
  
  # Add correlation values and significance markers
  for(i in 1:nrow(correlation_matrix)) {
    for(j in 1:ncol(correlation_matrix)) {
      if(!is.na(correlation_matrix[i,j])) {
        text_val <- sprintf("%.2f", correlation_matrix[i,j])
        
        # Add significance markers
        if(!is.na(p_value_matrix[i,j])) {
          if(p_value_matrix[i,j] < 0.001) {
            text_val <- paste0(text_val, "***")
          } else if(p_value_matrix[i,j] < 0.01) {
            text_val <- paste0(text_val, "**")
          } else if(p_value_matrix[i,j] < 0.05) {
            text_val <- paste0(text_val, "*")
          }
        }
        
        text_matrix[i,j] <- text_val
      }
    }
  }
  
  # Create the heatmap
  pdf(filename, width = 14, height = 12)
  pheatmap(
    correlation_matrix,
    display_numbers = text_matrix,
    main = title,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 8,
    fontsize_col = 8,
    fontsize_number = 7,
    cluster_rows = TRUE,
    cluster_cols = TRUE
    # Removed angle_col parameter that was causing the error
  )
  dev.off()
}

# Create heatmaps
create_correlation_heatmap(
  all_thickness_correlations, 
  "Correlation Between Wearable Parameters and OCTA Thickness\n* p<0.05, ** p<0.01, *** p<0.001",
  "thickness_correlation_heatmap.pdf"
)

create_correlation_heatmap(
  all_bloodflow_correlations, 
  "Correlation Between Wearable Parameters and OCTA Blood Flow\n* p<0.05, ** p<0.01, *** p<0.001",
  "bloodflow_correlation_heatmap.pdf"
)

# Function to create a forest plot
# 创建顶级相关性数据集 - 在调用create_forest_plot之前添加
# 可以选择按相关系数绝对值大小或显著性筛选
top_thickness_correlations <- all_thickness_correlations %>%
  filter(P_value < 0.05) %>%  # 只保留显著的相关性
  arrange(desc(abs(Correlation)))  # 按相关系数绝对值排序

top_bloodflow_correlations <- all_bloodflow_correlations %>%
  filter(P_value < 0.05) %>%
  arrange(desc(abs(Correlation)))

# 如果结果太多，可以只保留前几个
if(nrow(top_thickness_correlations) > 20) {
  top_thickness_correlations <- head(top_thickness_correlations, 20)
}

if(nrow(top_bloodflow_correlations) > 20) {
  top_bloodflow_correlations <- head(top_bloodflow_correlations, 20)
}

create_forest_plot <- function(correlations, title, filename) {
  # Create more readable variable labels
  correlations <- correlations %>%
    mutate(
      Wearable_Label = case_when(
        Wearable_Variable == "mean_hr" ~ "Mean Heart Rate",
        Wearable_Variable == "mean_rhr" ~ "Mean Resting Heart Rate",
        Wearable_Variable == "mean_bo" ~ "Mean Blood Oxygen",
        Wearable_Variable == "total_steps" ~ "Mean Total Steps",
        Wearable_Variable == "total_sleep" ~ "Mean Total Sleep",
        TRUE ~ Wearable_Variable
      ),
      Significance = ifelse(P_value < 0.05, "Significant", "Not Significant"),
      Pair_Label = paste(OCTA_Variable, "-", Wearable_Label),
      Group_Pair = paste(Group, Pair_Label)
    ) %>%
    arrange(Pair_Label, desc(Group))
  
  # Create the forest plot
  p <- ggplot(correlations, 
              aes(y = Group_Pair, x = Correlation, 
                  xmin = CI_Lower, xmax = CI_Upper, 
                  color = interaction(Group, Significance))) +
    geom_point(size = 3) +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Diabetes.Significant" = "#D6604D", 
                                  "Diabetes.Not Significant" = "#FDDBC7",
                                  "No Diabetes.Significant" = "#4393C3", 
                                  "No Diabetes.Not Significant" = "#D1E5F0"),
                       name = "Group and Significance") +
    labs(
      title = title,
      subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
      x = "Correlation Coefficient",
      y = ""
    ) +
    theme_custom +
    theme(
      axis.text.y = element_text(size = 8),
      plot.margin = unit(c(1, 3, 1, 1), "cm")  # Fixed the margin syntax
    ) +
    scale_x_continuous(limits = c(min(correlations$CI_Lower) - 0.1, 
                                  max(correlations$CI_Upper) + 0.3))
  
  # Add text annotations for correlation values
  p <- p + 
    geom_text(
      aes(x = max(CI_Upper) + 0.15, 
          label = sprintf("%.2f [%.2f, %.2f]%s", 
                          Correlation, 
                          CI_Lower, 
                          CI_Upper,
                          ifelse(P_value < 0.001, "***", 
                                 ifelse(P_value < 0.01, "**", 
                                        ifelse(P_value < 0.05, "*", ""))))),
      hjust = 0, size = 3
    )
  
  # Save and return the plot
  ggsave(filename, p, width = 12, height = length(unique(correlations$Pair_Label)) * 0.8)
  return(p)
}

# Create forest plots
thickness_forest <- create_forest_plot(
  top_thickness_correlations,
  "Correlation Between Wearable Parameters and OCTA Thickness",
  "thickness_correlation_forest.pdf"
)

bloodflow_forest <- create_forest_plot(
  top_bloodflow_correlations,
  "Correlation Between Wearable Parameters and OCTA Blood Flow",
  "bloodflow_correlation_forest.pdf"
)

# Step 11: Create scatter plots for top correlations

# Function to create a scatter plot for a correlation pair
create_scatter_plot <- function(data_dm, data_nodm, octa_var, wearable_var, 
                                dm_corr, nodm_corr, filename) {
  # Prepare combined data
  dm_data <- data_dm %>% 
    dplyr::select(subject_id, all_of(c(octa_var, wearable_var))) %>%
    mutate(Group = "Diabetes")
  
  nodm_data <- data_nodm %>% 
    dplyr::select(subject_id, all_of(c(octa_var, wearable_var))) %>%
    mutate(Group = "No Diabetes")
  
  combined_data <- bind_rows(dm_data, nodm_data)
  
  # Find the correlation values
  dm_r <- dm_corr$Correlation[dm_corr$OCTA_Variable == octa_var & 
                                dm_corr$Wearable_Variable == wearable_var]
  dm_p <- dm_corr$P_value[dm_corr$OCTA_Variable == octa_var & 
                            dm_corr$Wearable_Variable == wearable_var]
  
  nodm_r <- nodm_corr$Correlation[nodm_corr$OCTA_Variable == octa_var & 
                                    nodm_corr$Wearable_Variable == wearable_var]
  nodm_p <- nodm_corr$P_value[nodm_corr$OCTA_Variable == octa_var & 
                                nodm_corr$Wearable_Variable == wearable_var]
  
  # Create a prettier wearable variable label
  wearable_label <- case_when(
    wearable_var == "mean_hr" ~ "Mean Heart Rate",
    wearable_var == "mean_rhr" ~ "Mean Resting Heart Rate",
    wearable_var == "mean_bo" ~ "Mean Blood Oxygen",
    wearable_var == "total_steps" ~ "Mean Total Steps",
    wearable_var == "total_sleep" ~ "Mean Total Sleep",
    TRUE ~ wearable_var
  )
  
  # Create the scatter plot
  p <- ggplot(combined_data, aes(x = .data[[wearable_var]], 
                                 y = .data[[octa_var]], 
                                 color = Group))+
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(group = Group)) +
    scale_color_manual(values = c("Diabetes" = "#D6604D", "No Diabetes" = "#4393C3")) +
    labs(
      title = paste(octa_var, "vs.", wearable_label),
      subtitle = sprintf("Diabetes: r = %.2f, p = %.3f; No Diabetes: r = %.2f, p = %.3f", 
                         dm_r, dm_p, nodm_r, nodm_p),
      x = wearable_label,
      y = octa_var
    ) +
    theme_custom
  
  # Save the plot
  ggsave(filename, p, width = 8, height = 6)
  return(p)
}

# 首先检查每个组有多少显著相关记录
diabetes_count <- sum(top_thickness_correlations$Group == "Diabetes")
no_diabetes_count <- sum(top_thickness_correlations$Group == "No Diabetes")

cat("Number of significant correlations in Diabetes group:", diabetes_count, "\n")
cat("Number of significant correlations in No Diabetes group:", no_diabetes_count, "\n")

# 针对实际存在的组执行循环 - 厚度变量
if(diabetes_count > 0) {
  # 处理糖尿病组
  cat("Creating thickness scatter plots for Diabetes group\n")
  for(i in 1:min(5, diabetes_count)) {
    pair <- top_thickness_correlations %>%
      filter(Group == "Diabetes") %>%
      arrange(desc(abs(Correlation))) %>%
      slice(i)
    
    octa_var <- pair$OCTA_Variable
    wearable_var <- pair$Wearable_Variable
    
    cat("  Plot", i, ":", octa_var, "vs", wearable_var, "\n")
    
    create_scatter_plot(
      wearable_thickness_dm,
      wearable_thickness_nodm,
      octa_var,
      wearable_var,
      thickness_dm_correlations,
      thickness_nodm_correlations,
      paste0("scatter_thickness_diabetes_", gsub("[^a-zA-Z0-9]", "_", octa_var), "_", wearable_var, ".pdf")
    )
  }
}

if(no_diabetes_count > 0) {
  # 处理非糖尿病组
  cat("Creating thickness scatter plots for No Diabetes group\n")
  for(i in 1:min(5, no_diabetes_count)) {
    pair <- top_thickness_correlations %>%
      filter(Group == "No Diabetes") %>%
      arrange(desc(abs(Correlation))) %>%
      slice(i)
    
    octa_var <- pair$OCTA_Variable
    wearable_var <- pair$Wearable_Variable
    
    cat("  Plot", i, ":", octa_var, "vs", wearable_var, "\n")
    
    create_scatter_plot(
      wearable_thickness_dm,
      wearable_thickness_nodm,
      octa_var,
      wearable_var,
      thickness_dm_correlations,
      thickness_nodm_correlations,
      paste0("scatter_thickness_nodiabetes_", gsub("[^a-zA-Z0-9]", "_", octa_var), "_", wearable_var, ".pdf")
    )
  }
}

# 对血流变量重复同样的过程
diabetes_count_bf <- sum(top_bloodflow_correlations$Group == "Diabetes")
no_diabetes_count_bf <- sum(top_bloodflow_correlations$Group == "No Diabetes")

cat("\nNumber of significant correlations in Diabetes group (blood flow):", diabetes_count_bf, "\n")
cat("Number of significant correlations in No Diabetes group (blood flow):", no_diabetes_count_bf, "\n")

# 针对实际存在的组执行循环 - 血流变量
if(diabetes_count_bf > 0) {
  # 处理糖尿病组
  cat("Creating blood flow scatter plots for Diabetes group\n")
  for(i in 1:min(5, diabetes_count_bf)) {
    pair <- top_bloodflow_correlations %>%
      filter(Group == "Diabetes") %>%
      arrange(desc(abs(Correlation))) %>%
      slice(i)
    
    octa_var <- pair$OCTA_Variable
    wearable_var <- pair$Wearable_Variable
    
    cat("  Plot", i, ":", octa_var, "vs", wearable_var, "\n")
    
    create_scatter_plot(
      wearable_bloodflow_dm,
      wearable_bloodflow_nodm,
      octa_var,
      wearable_var,
      bloodflow_dm_correlations,
      bloodflow_nodm_correlations,
      paste0("scatter_bloodflow_diabetes_", gsub("[^a-zA-Z0-9]", "_", octa_var), "_", wearable_var, ".pdf")
    )
  }
}

if(no_diabetes_count_bf > 0) {
  # 处理非糖尿病组
  cat("Creating blood flow scatter plots for No Diabetes group\n")
  for(i in 1:min(5, no_diabetes_count_bf)) {
    pair <- top_bloodflow_correlations %>%
      filter(Group == "No Diabetes") %>%
      arrange(desc(abs(Correlation))) %>%
      slice(i)
    
    octa_var <- pair$OCTA_Variable
    wearable_var <- pair$Wearable_Variable
    
    cat("  Plot", i, ":", octa_var, "vs", wearable_var, "\n")
    
    create_scatter_plot(
      wearable_bloodflow_dm,
      wearable_bloodflow_nodm,
      octa_var,
      wearable_var,
      bloodflow_dm_correlations,
      bloodflow_nodm_correlations,
      paste0("scatter_bloodflow_nodiabetes_", gsub("[^a-zA-Z0-9]", "_", octa_var), "_", wearable_var, ".pdf")
    )
  }
}

# Step 12: Multiple regression models for key OCTA variables

# Function to build multiple regression models
build_regression_models <- function(data, octa_var, wearable_vars, covariates, group_name) {
  # Prepare formula
  formula_str <- paste(octa_var, "~", paste(c(wearable_vars, covariates), collapse = " + "))
  
  # Create regression model
  model <- tryCatch({
    lm(as.formula(formula_str), data = data)
  }, error = function(e) {
    cat("Error in creating model for", octa_var, "in", group_name, "group:", e$message, "\n")
    return(NULL)
  })
  
  # If model was successfully created
  if(!is.null(model)) {
    # Extract coefficients and confidence intervals
    coef_table <- summary(model)$coefficients
    conf_int <- confint(model)
    
    # Combine into a data frame
    result <- data.frame(
      OCTA_Variable = octa_var,
      Variable = rownames(coef_table),
      Coefficient = coef_table[, "Estimate"],
      Std_Error = coef_table[, "Std. Error"],
      t_value = coef_table[, "t value"],
      P_value = coef_table[, "Pr(>|t|)"],
      Lower_CI = conf_int[, "2.5 %"],
      Upper_CI = conf_int[, "97.5 %"],
      Group = group_name,
      stringsAsFactors = FALSE
    ) %>% filter(Variable != "(Intercept)")
    
    # Add model performance metrics
    model_summary <- summary(model)
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared
    
    # Calculate RMSE
    rmse <- sqrt(mean(model_summary$residuals^2))
    
    # Return results and model
    return(list(
      coef_data = result,
      model = model,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      rmse = rmse
    ))
  } else {
    return(NULL)
  }
}

# 修改后的函数用于选择顶级OCTA变量进行回归建模
select_top_octa_vars <- function(correlations, n = 5) {
  # 如果不是数据框，转换为数据框
  corr_df <- as.data.frame(correlations)
  
  # 过滤掉NA值
  corr_df <- corr_df[!is.na(corr_df$Correlation), ]
  
  # 添加绝对相关性列
  corr_df$abs_corr <- abs(corr_df$Correlation)
  
  # 为每个OCTA变量计算最大绝对相关性
  result <- aggregate(abs_corr ~ OCTA_Variable, data = corr_df, FUN = max)
  
  # 按降序排列相关性
  result <- result[order(-result$abs_corr), ]
  
  # 选择前n个变量
  if(nrow(result) > n) {
    top_vars <- result$OCTA_Variable[1:n]
  } else {
    top_vars <- result$OCTA_Variable
  }
  
  return(top_vars)
}

# 选择顶级变量
top_thickness_vars <- select_top_octa_vars(all_thickness_correlations)
top_bloodflow_vars <- select_top_octa_vars(all_bloodflow_correlations)

# 为顶级厚度变量建立模型
thickness_dm_models <- list()
thickness_nodm_models <- list()
thickness_dm_results <- data.frame()
thickness_nodm_results <- data.frame()
thickness_model_performance <- data.frame()

for(var in top_thickness_vars) {
  # 糖尿病组
  dm_result <- build_regression_models(
    wearable_thickness_dm, var, wearable_vars, covariate_vars, "Diabetes"
  )
  
  if(!is.null(dm_result)) {
    thickness_dm_models[[var]] <- dm_result$model
    thickness_dm_results <- rbind(thickness_dm_results, dm_result$coef_data)
    
    # 添加性能指标
    thickness_model_performance <- rbind(
      thickness_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "Diabetes",
        R_squared = dm_result$r_squared,
        Adj_R_squared = dm_result$adj_r_squared,
        RMSE = dm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # 非糖尿病组
  nodm_result <- build_regression_models(
    wearable_thickness_nodm, var, wearable_vars, covariate_vars, "No Diabetes"
  )
  
  if(!is.null(nodm_result)) {
    thickness_nodm_models[[var]] <- nodm_result$model
    thickness_nodm_results <- rbind(thickness_nodm_results, nodm_result$coef_data)
    
    # 添加性能指标
    thickness_model_performance <- rbind(
      thickness_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "No Diabetes",
        R_squared = nodm_result$r_squared,
        Adj_R_squared = nodm_result$adj_r_squared,
        RMSE = nodm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }
}

# 为顶级血流变量建立模型
bloodflow_dm_models <- list()
bloodflow_nodm_models <- list()
bloodflow_dm_results <- data.frame()
bloodflow_nodm_results <- data.frame()
bloodflow_model_performance <- data.frame()

for(var in top_bloodflow_vars) {
  # 糖尿病组
  dm_result <- build_regression_models(
    wearable_bloodflow_dm, var, wearable_vars, covariate_vars, "Diabetes"
  )
  
  if(!is.null(dm_result)) {
    bloodflow_dm_models[[var]] <- dm_result$model
    bloodflow_dm_results <- rbind(bloodflow_dm_results, dm_result$coef_data)
    
    # 添加性能指标
    bloodflow_model_performance <- rbind(
      bloodflow_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "Diabetes",
        R_squared = dm_result$r_squared,
        Adj_R_squared = dm_result$adj_r_squared,
        RMSE = dm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # 非糖尿病组
  nodm_result <- build_regression_models(
    wearable_bloodflow_nodm, var, wearable_vars, covariate_vars, "No Diabetes"
  )
  
  if(!is.null(nodm_result)) {
    bloodflow_nodm_models[[var]] <- nodm_result$model
    bloodflow_nodm_results <- rbind(bloodflow_nodm_results, nodm_result$coef_data)
    
    # 添加性能指标
    bloodflow_model_performance <- rbind(
      bloodflow_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "No Diabetes",
        R_squared = nodm_result$r_squared,
        Adj_R_squared = nodm_result$adj_r_squared,
        RMSE = nodm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }
}

# 合并所有回归结果
all_thickness_reg_results <- rbind(thickness_dm_results, thickness_nodm_results)
all_bloodflow_reg_results <- rbind(bloodflow_dm_results, bloodflow_nodm_results)

# 保存回归结果
write.csv(all_thickness_reg_results, "thickness_regression_results.csv", row.names = FALSE)
write.csv(all_bloodflow_reg_results, "bloodflow_regression_results.csv", row.names = FALSE)
write.csv(rbind(thickness_model_performance, bloodflow_model_performance), 
          "model_performance.csv", row.names = FALSE)


# Step 13: Create forest plots for regression coefficients

# Function to create a forest plot for regression coefficients
create_reg_forest_plot <- function(reg_results, title, filename) {
  # Focus on wearable variables only
  wearable_results <- reg_results %>%
    filter(Variable %in% wearable_vars) %>%
    mutate(
      Wearable_Label = case_when(
        Variable == "mean_hr" ~ "Mean Heart Rate",
        Variable == "mean_rhr" ~ "Mean Resting Heart Rate",
        Variable == "mean_bo" ~ "Mean Blood Oxygen",
        Variable == "total_steps" ~ "Mean Total Steps",
        Variable == "total_sleep" ~ "Mean Total Sleep",
        TRUE ~ Variable
      ),
      Significance = ifelse(P_value < 0.05, "Significant", "Not Significant"),
      Pair_Label = paste(OCTA_Variable, "-", Wearable_Label),
      Group_Pair = paste(Group, Pair_Label)
    ) %>%
    arrange(Pair_Label, desc(Group))
  
  # Create the forest plot
  p <- ggplot(wearable_results, 
              aes(y = Group_Pair, x = Coefficient, 
                  xmin = Lower_CI, xmax = Upper_CI, 
                  color = interaction(Group, Significance))) +
    geom_point(size = 3) +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Diabetes.Significant" = "#D6604D", 
                                  "Diabetes.Not Significant" = "#FDDBC7",
                                  "No Diabetes.Significant" = "#4393C3", 
                                  "No Diabetes.Not Significant" = "#D1E5F0"),
                       name = "Group and Significance") +
    labs(
      title = title,
      subtitle = "Regression Coefficients with 95% CI\n* p<0.05, ** p<0.01, *** p<0.001",
      x = "Regression Coefficient",
      y = ""
    ) +
    theme_custom +
    theme(
      axis.text.y = element_text(size = 8),
      plot.margin = unit(c(1, 3, 1, 1), "cm")  # 修改了这里，使用unit函数
    ) +
    scale_x_continuous(limits = c(min(wearable_results$Lower_CI) - 0.1, 
                                  max(wearable_results$Upper_CI) + 0.3))
  
  # Add text annotations for coefficient values
  p <- p + 
    geom_text(
      aes(x = max(Upper_CI) + 0.15, 
          label = sprintf("%.2f [%.2f, %.2f]%s", 
                          Coefficient, 
                          Lower_CI, 
                          Upper_CI,
                          ifelse(P_value < 0.001, "***", 
                                 ifelse(P_value < 0.01, "**", 
                                        ifelse(P_value < 0.05, "*", ""))))),
      hjust = 0, size = 3
    )
  
  # Save and return the plot
  ggsave(filename, p, width = 12, height = length(unique(wearable_results$Pair_Label)) * 0.8)
  return(p)
}

# Create regression forest plots
thickness_reg_forest <- create_reg_forest_plot(
  all_thickness_reg_results,
  "Wearable Parameters' Association with OCTA Thickness",
  "thickness_regression_forest.pdf"
)

bloodflow_reg_forest <- create_reg_forest_plot(
  all_bloodflow_reg_results,
  "Wearable Parameters' Association with OCTA Blood Flow",
  "bloodflow_regression_forest.pdf"
)

# Step 14: Model performance comparison between groups

# Function to create performance comparison plot
create_performance_plot <- function(performance_data, metric = "R_squared", title, filename) {
  # Prepare data
  plot_data <- performance_data %>%
    arrange(desc(!!sym(metric)))
  
  # Nicer metric label
  metric_label <- case_when(
    metric == "R_squared" ~ "R²",
    metric == "Adj_R_squared" ~ "Adjusted R²",
    metric == "RMSE" ~ "RMSE",
    TRUE ~ metric
  )
  
  # Create bar plot
  p <- ggplot(plot_data, aes(x = reorder(OCTA_Variable, !!sym(metric)), 
                             y = !!sym(metric), 
                             fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Diabetes" = "#D6604D", "No Diabetes" = "#4393C3")) +
    labs(
      title = title,
      subtitle = paste("Comparison of", metric_label, "Between Diabetic and Non-diabetic Groups"),
      x = "OCTA Variables",
      y = metric_label,
      fill = "Group"
    ) +
    theme_custom +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Save and return the plot
  ggsave(filename, p, width = 10, height = 6)
  return(p)
}

# Create performance comparison plots
thickness_perf_r2 <- create_performance_plot(
  thickness_model_performance,
  "R_squared",
  "OCTA Thickness Model Performance",
  "thickness_model_r2_comparison.pdf"
)

thickness_perf_rmse <- create_performance_plot(
  thickness_model_performance,
  "RMSE",
  "OCTA Thickness Model Performance",
  "thickness_model_rmse_comparison.pdf"
)

bloodflow_perf_r2 <- create_performance_plot(
  bloodflow_model_performance,
  "R_squared",
  "OCTA Blood Flow Model Performance",
  "bloodflow_model_r2_comparison.pdf"
)

bloodflow_perf_rmse <- create_performance_plot(
  bloodflow_model_performance,
  "RMSE",
  "OCTA Blood Flow Model Performance",
  "bloodflow_model_rmse_comparison.pdf"
)

# Step 15: Predicted vs Actual plots for best models

# Function to find best model
find_best_model <- function(performance_data, group, type = "thickness") {
  # Find variable with highest R²
  best_var <- performance_data %>%
    filter(Group == group) %>%
    arrange(desc(R_squared)) %>%
    slice(1) %>%
    pull(OCTA_Variable)
  
  if(type == "thickness") {
    if(group == "Diabetes") {
      return(list(variable = best_var, model = thickness_dm_models[[best_var]]))
    } else {
      return(list(variable = best_var, model = thickness_nodm_models[[best_var]]))
    }
  } else {
    if(group == "Diabetes") {
      return(list(variable = best_var, model = bloodflow_dm_models[[best_var]]))
    } else {
      return(list(variable = best_var, model = bloodflow_nodm_models[[best_var]]))
    }
  }
}

# Function to create predicted vs actual plot
create_pred_actual_plot <- function(model, data, var, group, title, filename) {
  # Calculate predictions
  data$predicted <- predict(model, newdata = data)
  
  # Create plot
  p <- ggplot(data, aes(x = predicted, y = .data[[var]])) +
    geom_point(alpha = 0.7, color = ifelse(group == "Diabetes", "#D6604D", "#4393C3")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = title,
      subtitle = paste(group, "Group"),
      x = "Predicted Values",
      y = "Actual Values"
    ) +
    annotate("text", x = min(data$predicted, na.rm = TRUE), 
             y = max(data[[var]], na.rm = TRUE),
             label = sprintf("R² = %.3f\nRMSE = %.3f", 
                             summary(model)$r.squared,
                             sqrt(mean(resid(model)^2))),
             hjust = 0, vjust = 1) +
    theme_custom
  
  # Save and return the plot
  ggsave(filename, p, width = 8, height = 6)
  return(p)
}

# Create predicted vs actual plots for best models
# Thickness - Diabetic
best_thickness_dm <- find_best_model(thickness_model_performance, "Diabetes", "thickness")
if(!is.null(best_thickness_dm$model)) {
  thickness_dm_pred_plot <- create_pred_actual_plot(
    best_thickness_dm$model,
    wearable_thickness_dm,
    best_thickness_dm$variable,
    "Diabetes",
    paste("Predicted vs Actual", best_thickness_dm$variable),
    "thickness_dm_pred_actual.pdf"
  )
}

# Thickness - Non-diabetic
best_thickness_nodm <- find_best_model(thickness_model_performance, "No Diabetes", "thickness")
if(!is.null(best_thickness_nodm$model)) {
  thickness_nodm_pred_plot <- create_pred_actual_plot(
    best_thickness_nodm$model,
    wearable_thickness_nodm,
    best_thickness_nodm$variable,
    "No Diabetes",
    paste("Predicted vs Actual", best_thickness_nodm$variable),
    "thickness_nodm_pred_actual.pdf"
  )
}

# Blood flow - Diabetic
best_bloodflow_dm <- find_best_model(bloodflow_model_performance, "Diabetes", "bloodflow")
if(!is.null(best_bloodflow_dm$model)) {
  bloodflow_dm_pred_plot <- create_pred_actual_plot(
    best_bloodflow_dm$model,
    wearable_bloodflow_dm,
    best_bloodflow_dm$variable,
    "Diabetes",
    paste("Predicted vs Actual", best_bloodflow_dm$variable),
    "bloodflow_dm_pred_actual.pdf"
  )
}

# Blood flow - Non-diabetic
best_bloodflow_nodm <- find_best_model(bloodflow_model_performance, "No Diabetes", "bloodflow")
if(!is.null(best_bloodflow_nodm$model)) {
  bloodflow_nodm_pred_plot <- create_pred_actual_plot(
    best_bloodflow_nodm$model,
    wearable_bloodflow_nodm,
    best_bloodflow_nodm$variable,
    "No Diabetes",
    paste("Predicted vs Actual", best_bloodflow_nodm$variable),
    "bloodflow_nodm_pred_actual.pdf"
  )
}

# Step 16: Summary of results

# Significant correlations
cat("\n\n=== SUMMARY OF RESULTS ===\n\n")

cat("Top significant OCTA thickness correlations with wearable parameters:\n")
print(all_thickness_correlations %>% 
        filter(P_value < 0.05) %>% 
        arrange(P_value) %>%
        head(10) %>%
        dplyr::select(Group, OCTA_Variable, Wearable_Variable, Correlation, P_value))

cat("\nTop significant OCTA blood flow correlations with wearable parameters:\n")
print(all_bloodflow_correlations %>% 
        filter(P_value < 0.05) %>% 
        arrange(P_value) %>%
        head(10) %>%
        dplyr::select(Group, OCTA_Variable, Wearable_Variable, Correlation, P_value))

# Significant regression coefficients
cat("\nTop significant regression coefficients for OCTA thickness models:\n")
print(all_thickness_reg_results %>% 
        filter(P_value < 0.05, Variable %in% wearable_vars) %>% 
        arrange(P_value) %>%
        head(10) %>%
        dplyr::select(Group, OCTA_Variable, Variable, Coefficient, P_value))

cat("\nTop significant regression coefficients for OCTA blood flow models:\n")
print(all_bloodflow_reg_results %>% 
        filter(P_value < 0.05, Variable %in% wearable_vars) %>% 
        arrange(P_value) %>%
        head(10) %>%
        dplyr::select(Group, OCTA_Variable, Variable, Coefficient, P_value))

# Model performance summary
cat("\nBest performing models for OCTA thickness:\n")
print(thickness_model_performance %>% 
        group_by(Group) %>%
        arrange(desc(R_squared)) %>%
        slice(1) %>%
        ungroup() %>%
        dplyr::select(Group, OCTA_Variable, R_squared, RMSE))

cat("\nBest performing models for OCTA blood flow:\n")
print(bloodflow_model_performance %>% 
        group_by(Group) %>%
        arrange(desc(R_squared)) %>%
        slice(1) %>%
        ungroup() %>%
        dplyr::select(Group, OCTA_Variable, R_squared, RMSE))
