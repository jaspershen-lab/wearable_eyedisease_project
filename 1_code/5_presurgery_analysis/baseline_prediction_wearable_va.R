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
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")

dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/va", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/va")

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
    median_rhr = mean(median_rhr_1, na.rm = TRUE), 
    mean_bo = mean(mean_bo, na.rm = TRUE),           # Mean blood oxygen
    min_bo = mean(min_bo, na.rm = TRUE), 
    max_bo = mean(mean_bo, na.rm = TRUE), 
    total_steps = mean(steps_total, na.rm = TRUE),   # Mean total steps
    max_steps=mean(steps_max, na.rm = TRUE),
    total_sleep = mean(total_sleep, na.rm = TRUE),   # Mean total sleep
    deep_sleep = mean(deep_sleep, na.rm = TRUE), 
    
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
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female"))
  ) %>%
  dplyr::select(mean_hr, mean_rhr, mean_bo, total_steps, total_sleep, 
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
      mean_hr ~ "Mean Heart Rate",
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
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Comparison by Diabetes Status**") %>%
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
    cols = c(mean_hr, mean_rhr, mean_bo, total_steps, total_sleep),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "DR", "Cataract"),
    parameter_label = case_when(
      parameter == "mean_hr" ~ "Mean Heart Rate",
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

# If you want to add p-values, first define the p_values dataframe correctly:
# Create a separate data frame for p-values
calculate_p_values <- function() {
  p_values_df <- data.frame(parameter = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each parameter
  for(param in unique(long_data$parameter)) {
    # Extract data for this parameter
    param_data <- long_data[long_data$parameter == param, ]
    
    # Only perform test if both groups have data
    if(length(unique(param_data$dm_status)) > 1) {
      # Perform Wilcoxon rank sum test
      test_result <- wilcox.test(value ~ dm_status, data = param_data)
      
      # Add result to data frame
      p_values_df <- rbind(p_values_df, 
                           data.frame(parameter = param, 
                                      p_value = test_result$p.value,
                                      stringsAsFactors = FALSE))
    }
  }
  
  # Format p-values for display
  p_values_df$parameter_label <- case_when(
    p_values_df$parameter == "mean_hr" ~ "Mean Heart Rate",
    p_values_df$parameter == "mean_rhr" ~ "Mean Resting Heart Rate",
    p_values_df$parameter == "mean_bo" ~ "Mean Blood Oxygen",
    p_values_df$parameter == "total_steps" ~ "Mean Total Steps",
    p_values_df$parameter == "total_sleep" ~ "Mean Total Sleep",
    p_values_df$parameter == "pre_vision" ~ "Pre-surgery Vision"
  )
  
  p_values_df$p_text <- sprintf("p = %.3f", p_values_df$p_value)
  p_values_df$p_text <- ifelse(p_values_df$p_value < 0.05, 
                               paste0(p_values_df$p_text, "*"), 
                               p_values_df$p_text)
  
  return(p_values_df)
}

# Calculate p-values
p_values <- calculate_p_values()
print(p_values)

# Now create boxplots with p-values if the function ran successfully
if(exists("p_values") && nrow(p_values) > 0) {
  comparison_boxplots_with_p <- ggplot(long_data, aes(x = dm_status, y = value, fill = dm_status)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ parameter_label, scales = "free_y") +
    # Use custom colors directly from 100_tools.R
    scale_fill_manual(values = diabetes_colors) +
    labs(
      title = "Comparison of Pre-surgery Parameters Between Diabetic and Non-diabetic Patients",
      subtitle = "* indicates p < 0.05",
      x = "",
      y = "",
      fill = "Group"
    ) +
    geom_text(data = p_values, aes(x = 1.5, y = Inf, label = p_text), 
              vjust = -0.5, hjust = 0.5, inherit.aes = FALSE) +
    # MODIFIED: Use theme_custom from 100_tools.R
    theme_custom
  
  print(comparison_boxplots_with_p)
  ggsave("diabetic_nondiabetic_comparison_boxplots_with_p.pdf", comparison_boxplots_with_p, width = 14, height = 10)
}

# Step 7: Correlation analysis for both groups
# Function to calculate correlation with vision
calculate_correlation_vision <- function(data, var_name) {
  # Pearson correlation coefficient
  correlation <- cor(data[[var_name]], data$pre_vision, use = "complete.obs")
  
  # Fisher z transformation for confidence interval
  n <- sum(!is.na(data[[var_name]]) & !is.na(data$pre_vision))
  z <- 0.5 * log((1 + correlation) / (1 - correlation))
  se <- 1 / sqrt(n - 3)
  ci_lower <- tanh(z - 1.96 * se)
  ci_upper <- tanh(z + 1.96 * se)
  
  # Calculate p-value
  t_stat <- correlation * sqrt((n - 2) / (1 - correlation^2))
  p_value <- 2 * pt(-abs(t_stat), df = n - 2)
  
  return(data.frame(
    Variable = var_name,
    Correlation = correlation,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    P_value = p_value,
    N = n
  ))
}

# Calculate correlations for both groups
wearable_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")

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

# Combine results
combined_correlations <- rbind(diabetic_correlations, nondiabetic_correlations)

# Add human-readable variable names
combined_correlations <- combined_correlations %>%
  mutate(
    Variable_Label = case_when(
      Variable == "mean_hr" ~ "Mean Heart Rate",
      Variable == "mean_rhr" ~ "Mean Resting Heart Rate",
      Variable == "mean_bo" ~ "Mean Blood Oxygen",
      Variable == "total_steps" ~ "Mean Total Steps",
      Variable == "total_sleep" ~ "Mean Total Sleep"
    )
  )

# Print correlation results
print(combined_correlations)

# Step 8: Create forest plot comparing correlations between groups
# Prepare forest plot data
combined_forest_data <- combined_correlations %>%
  mutate(
    Significance = ifelse(P_value < 0.05, "Significant", "Not Significant"),
    Correlation_Direction = ifelse(Correlation > 0, "Positive", "Negative"),
    Group_Variable = paste(Group, Variable_Label)
  ) %>%
  arrange(Variable_Label, desc(Group))

# Function to create lighter version of a color for non-significant results
lighten_color <- function(color, factor = 0.7) {
  rgb_col <- col2rgb(color) / 255
  lighter <- rgb_col * factor + (1 - factor)
  rgb(lighter[1], lighter[2], lighter[3])
}

# Create the forest plot with 100_tools colors
comparison_forest_plot <- ggplot(combined_forest_data, 
                                 aes(y = Group_Variable, x = Correlation, 
                                     xmin = CI_Lower, xmax = CI_Upper, 
                                     color = interaction(Group, Significance))) +
  geom_point(size = 4) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # 直接定义颜色，只有基本分组用source的颜色
  scale_color_manual(values = c("Diabetes.Significant" = "#D6604D", 
                                "Diabetes.Not Significant" = "#FDDBC7",
                                "Cataract.Significant" = "#4393C3", 
                                "Cataract.Not Significant" = "#D1E5F0"),
                     name = "Group and Significance") +
  labs(
    title = "Correlation Between Pre-surgery Wearable Parameters and Vision",
    subtitle = "Comparison Between Diabetic and Non-diabetic Patients",
    x = "Correlation Coefficient",
    y = ""
  ) +
  # MODIFIED: Use theme_custom from 100_tools.R
  theme_custom +
  annotate("text", x = max(combined_forest_data$CI_Upper) + 0.1, 
           y = 1:nrow(combined_forest_data), 
           label = sprintf("%.3f [%.3f, %.3f]", 
                           combined_forest_data$Correlation, 
                           combined_forest_data$CI_Lower, 
                           combined_forest_data$CI_Upper), 
           hjust = 0, size = 3)

# Display forest plot
print(comparison_forest_plot)

# Save forest plot
ggsave("diabetic_nondiabetic_correlation_comparison.pdf", comparison_forest_plot, width = 12, height = 8)





# ===============================================================
# NEW CODE: RCS Analysis for Diabetic and Non-Diabetic Groups
# ===============================================================

# Step 4: Convert pre_vision to binary
# Calculate the median vision to use as cutoff
vision_median <- median(data_imputed$pre_vision, na.rm = TRUE)
cat("Median pre-surgery vision:", vision_median, "\n")

# Create binary vision status (1 = good vision (above median), 0 = poor vision (below or equal to median))
data_imputed$vision_status <- ifelse(data_imputed$pre_vision > 0.5, 1, 0)

# Print counts for each group
vision_table <- table(data_imputed$vision_status)
cat("Vision status distribution:\n")
cat("Poor vision (0):", vision_table[1], "subjects\n")
cat("Good vision (1):", vision_table[2], "subjects\n")

# Step 5: Split data into diabetic and non-diabetic groups
diabetic_data <- data_imputed %>% filter(dm_2 == 1)
nondiabetic_data <- data_imputed %>% filter(dm_2 == 0)

# Print sample sizes
cat("Sample size of diabetic patients:", nrow(diabetic_data), "subjects\n")
cat("Sample size of non-diabetic patients:", nrow(nondiabetic_data), "subjects\n")

# Define the wearable variables
wearable_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")

# Create human-readable variable labels
var_labels <- c(
  mean_hr = "Mean Heart Rate",
  mean_rhr = "Mean Resting Heart Rate",
  mean_bo = "Mean Blood Oxygen",
  total_steps = "Mean Total Steps",
  total_sleep = "Mean Total Sleep"
)

# ===============================================================
# RCS Logistic Regression Analysis for Each Variable
# ===============================================================

# Setup datadist for diabetic and nondiabetic groups
dd_diabetic <- datadist(diabetic_data)
dd_nondiabetic <- datadist(nondiabetic_data)

# Number of knots for RCS
n_knots <-3

# Storage for results
diabetic_lrm_results <- list()
nondiabetic_lrm_results <- list()
combined_or_plots <- list()

# Loop through each wearable variable
for (var in wearable_vars) {
  var_label <- var_labels[var]
  cat("\nAnalyzing", var_label, "...\n")
  
  # ---- Diabetic Group Analysis ----
  options(datadist = "dd_diabetic")
  
  # Set reference point to median of the variable
  var_median <- median(diabetic_data[[var]], na.rm = TRUE)
  dd_diabetic$limits[[var]][2] <- var_median
  options(datadist = "dd_diabetic")
  
  # Fit logistic regression model with RCS for diabetic group
  formula_str <- paste("vision_status ~ rcs(", var, ", ", n_knots, ") + age + gender + bmi")
  diabetic_fit <- lrm(as.formula(formula_str), data = diabetic_data)
  
  # Get ANOVA results
  diabetic_anova <- anova(diabetic_fit)
  
  # Extract p-values - safely handle missing rows
  diabetic_nonlinear_p <- if("Nonlinear" %in% rownames(diabetic_anova)) {
    diabetic_anova["Nonlinear", "P"]
  } else {
    NA
  }
  
  # Generate OR predictions
  diabetic_or <- Predict(diabetic_fit, name = var, fun = exp, ref.zero = TRUE)
  diabetic_or_df <- as.data.frame(diabetic_or)
  diabetic_or_df$group <- "DR"
  
  # Create OR plot for diabetic group
  diabetic_or_plot <- ggplot(diabetic_or_df, aes(x = .data[[var]], y = yhat)) +
    geom_line(size = 1, color = "#D6604D") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#D6604D") +
    geom_hline(yintercept = 1, linetype = 2, size = 1) +
    labs(
      title = paste(var_label, "vs. Vision Status - Diabetes Group"),
      subtitle = paste("Non-linear P =", format(diabetic_nonlinear_p, digits = 3)),
      x = var_label,
      y = "OR (95% CI)"
    ) +
    theme_custom
  
  # Save results for diabetic group
  diabetic_lrm_results[[var]] <- list(
    fit = diabetic_fit,
    anova = diabetic_anova,
    or_predictions = diabetic_or_df,
    or_plot = diabetic_or_plot
  )
  
  # Print and save plot
  print(diabetic_or_plot)
  ggsave(paste0("or_", var, "_diabetic.pdf"), diabetic_or_plot, width = 8, height = 6)
  
  # ---- Non-Diabetic Group Analysis ----
  options(datadist = "dd_nondiabetic")
  
  # Set reference point to median of the variable
  var_median <- median(nondiabetic_data[[var]], na.rm = TRUE)
  dd_nondiabetic$limits[[var]][2] <- var_median
  options(datadist = "dd_nondiabetic")
  
  # Fit logistic regression model with RCS for non-diabetic group
  nondiabetic_fit <- lrm(as.formula(formula_str), data = nondiabetic_data)
  
  # Get ANOVA results
  nondiabetic_anova <- anova(nondiabetic_fit)
  
  # Extract p-values - safely handle missing rows
  nondiabetic_nonlinear_p <- if("Nonlinear" %in% rownames(nondiabetic_anova)) {
    nondiabetic_anova["Nonlinear", "P"]
  } else {
    NA
  }
  
  # Generate OR predictions
  nondiabetic_or <- Predict(nondiabetic_fit, name = var, fun = exp, ref.zero = TRUE)
  nondiabetic_or_df <- as.data.frame(nondiabetic_or)
  nondiabetic_or_df$group <- "Cataract"
  
  # Create OR plot for non-diabetic group
  nondiabetic_or_plot <- ggplot(nondiabetic_or_df, aes(x = .data[[var]], y = yhat)) +
    geom_line(size = 1, color = "#4393C3") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#4393C3") +
    geom_hline(yintercept = 1, linetype = 2, size = 1) +
    labs(
      title = paste(var_label, "vs. Vision Status - Cataract Group"),
      subtitle = paste("Non-linear P =", format(nondiabetic_nonlinear_p, digits = 3)),
      x = var_label,
      y = "OR (95% CI)"
    ) +
    theme_custom
  
  # Save results for non-diabetic group
  nondiabetic_lrm_results[[var]] <- list(
    fit = nondiabetic_fit,
    anova = nondiabetic_anova,
    or_predictions = nondiabetic_or_df,
    or_plot = nondiabetic_or_plot
  )
  
  # Print and save plot
  print(nondiabetic_or_plot)
  ggsave(paste0("or_", var, "_nondiabetic.pdf"), nondiabetic_or_plot, width = 8, height = 6)

}


# ===============================================================
# Create summary table of RCS analysis results
# ===============================================================

# Prepare data for summary table
lrm_summary <- data.frame(
  Variable = character(),
  Group = character(),
  NonlinearP = numeric(),
  OverallP = numeric(),
  stringsAsFactors = FALSE
)

for (var in wearable_vars) {
  var_label <- var_labels[var]
  
  # Extract p-values for diabetic group
  diabetic_anova <- diabetic_lrm_results[[var]]$anova
  
  diabetic_nonlinear_p <- if("Nonlinear" %in% rownames(diabetic_anova)) {
    diabetic_anova["Nonlinear", "P"]
  } else {
    NA
  }
  
  # For overall effect, check different possible row names
  diabetic_overall_p <- NA
  if(paste0("rcs(", var, ", ", n_knots, ")") %in% rownames(diabetic_anova)) {
    diabetic_overall_p <- diabetic_anova[paste0("rcs(", var, ", ", n_knots, ")"), "P"]
  } else if(var %in% rownames(diabetic_anova)) {
    diabetic_overall_p <- diabetic_anova[var, "P"]
  } else if("TOTAL" %in% rownames(diabetic_anova)) {
    diabetic_overall_p <- diabetic_anova["TOTAL", "P"]
  }
  
  # Extract p-values for non-diabetic group
  nondiabetic_anova <- nondiabetic_lrm_results[[var]]$anova
  
  nondiabetic_nonlinear_p <- if("Nonlinear" %in% rownames(nondiabetic_anova)) {
    nondiabetic_anova["Nonlinear", "P"]
  } else {
    NA
  }
  
  # For overall effect, check different possible row names
  nondiabetic_overall_p <- NA
  if(paste0("rcs(", var, ", ", n_knots, ")") %in% rownames(nondiabetic_anova)) {
    nondiabetic_overall_p <- nondiabetic_anova[paste0("rcs(", var, ", ", n_knots, ")"), "P"]
  } else if(var %in% rownames(nondiabetic_anova)) {
    nondiabetic_overall_p <- nondiabetic_anova[var, "P"]
  } else if("TOTAL" %in% rownames(nondiabetic_anova)) {
    nondiabetic_overall_p <- nondiabetic_anova["TOTAL", "P"]
  }
  
  # Add to summary table
  lrm_summary <- rbind(
    lrm_summary,
    data.frame(
      Variable = var_label,
      Group = "DR",
      NonlinearP = diabetic_nonlinear_p,
      OverallP = diabetic_overall_p
    ),
    data.frame(
      Variable = var_label,
      Group = "Cataract",
      NonlinearP = nondiabetic_nonlinear_p,
      OverallP = nondiabetic_overall_p
    )
  )
}

# Print the summary table
print(lrm_summary)




# Step 9: Multiple regression models for both groups
# Create separate regression models
diabetic_model <- lm(
  pre_vision ~ mean_hr + mean_rhr + mean_bo + total_steps  + age + gender + bmi,
  data = diabetic_data
)

diabetic_model <- glm(
  vision_status ~ mean_hr + mean_rhr + mean_bo + total_steps + total_sleep + age + gender + bmi,
  data = data_imputed
)

print(summary(diabetic_model))

nondiabetic_model <- lm(
  pre_vision ~ mean_hr + mean_rhr + mean_bo + total_steps + total_sleep + age + gender + bmi,
  data = nondiabetic_data
)

print(summary(nondiabetic_model))


combined_model <- lm(
  pre_vision ~ mean_hr + mean_rhr + mean_bo + total_steps + total_sleep + age + gender + bmi + dm_2,
  data = data_imputed  # 使用已经包含两组患者的完整数据集
)

summary(combined_model)



# Extract coefficients and confidence intervals
# For diabetic patients
diabetic_coef <- coef(diabetic_model)
diabetic_ci <- confint(diabetic_model)
diabetic_p <- summary(diabetic_model)$coefficients[, "Pr(>|t|)"]

diabetic_coef_data <- data.frame(
  Variable = names(diabetic_coef),
  Coefficient = diabetic_coef,
  Lower_CI = diabetic_ci[, "2.5 %"],
  Upper_CI = diabetic_ci[, "97.5 %"],
  P_value = diabetic_p,
  Group = "Diabetic"
) %>% filter(Variable != "(Intercept)")

# For non-diabetic patients
nondiabetic_coef <- coef(nondiabetic_model)
nondiabetic_ci <- confint(nondiabetic_model)
nondiabetic_p <- summary(nondiabetic_model)$coefficients[, "Pr(>|t|)"]

nondiabetic_coef_data <- data.frame(
  Variable = names(nondiabetic_coef),
  Coefficient = nondiabetic_coef,
  Lower_CI = nondiabetic_ci[, "2.5 %"],
  Upper_CI = nondiabetic_ci[, "97.5 %"],
  P_value = nondiabetic_p,
  Group = "Non-diabetic"
) %>% filter(Variable != "(Intercept)")

# Combine coefficient data
combined_coef_data <- rbind(diabetic_coef_data, nondiabetic_coef_data) %>%
  mutate(
    Variable_Label = case_when(
      Variable == "mean_hr" ~ "Mean Heart Rate",
      Variable == "mean_rhr" ~ "Mean Resting Heart Rate",
      Variable == "mean_bo" ~ "Mean Blood Oxygen",
      Variable == "total_steps" ~ "Mean Total Steps",
      Variable == "total_sleep" ~ "Mean Total Sleep",
      Variable == "age" ~ "Age",
      Variable == "gender" ~ "Gender",
      Variable == "bmi" ~ "BMI"
    ),
    Significance = ifelse(P_value < 0.05, "Significant", "Not Significant"),
    Group_Variable = paste(Group, Variable_Label)
  ) %>%
  arrange(Variable_Label, desc(Group))

# Create the forest plot for regression coefficients with 100_tools colors
# 重新定义显著性
combined_coef_data$Significance <- "Not Significant"
combined_coef_data$Significance[combined_coef_data$P_value < 0.05] <- "Significant"

# 确保组名一致
combined_coef_data$Group_for_color <- combined_coef_data$Group
combined_coef_data$Group_for_color[combined_coef_data$Group_for_color == "Diabetic"] <- "DR"
combined_coef_data$Group_for_color[combined_coef_data$Group_for_color == "Non-diabetic"] <- "Cataract"

# 创建正确的交互变量
combined_coef_data$Group_Sig <- paste(combined_coef_data$Group_for_color, combined_coef_data$Significance, sep=".")
# 然后在ggplot中使用这个新变量
coef_comparison_plot <- ggplot(combined_coef_data, 
                               aes(y = Group_Variable, x = Coefficient, 
                                   xmin = Lower_CI, xmax = Upper_CI, 
                                 color = Group_Sig)) +  # 使用新创建的变量
  geom_point(size = 4) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # MODIFIED: Use diabetes_colors from 100_tools.R
  scale_color_manual(values = c("DR.Significant" = "#D6604D", 
                                "DR.Not Significant" = "#FDDBC7",
                                "Cataract.Significant" = "#4393C3", 
                                "Cataract.Not Significant" = "#D1E5F0"),
                     name = "Group and Significance") +
  labs(
    title = "Association Between Pre-surgery Parameters and Vision",
    subtitle = "Regression Coefficients Comparison Between Diabetic and Non-diabetic Patients",
    x = "Regression Coefficient",
    y = ""
  ) +
  # MODIFIED: Use theme_custom from 100_tools.R
  theme_custom +
  annotate("text", x = max(abs(c(combined_coef_data$Lower_CI, combined_coef_data$Upper_CI))) * 1.2, 
           y = 1:nrow(combined_coef_data), 
           label = sprintf("%.3f [%.3f, %.3f]", 
                           combined_coef_data$Coefficient, 
                           combined_coef_data$Lower_CI, 
                           combined_coef_data$Upper_CI), 
           hjust = 0, size = 3)

# Display coefficient comparison plot
print(coef_comparison_plot)

# Save coefficient comparison plot
ggsave("diabetic_nondiabetic_regression_comparison.pdf", coef_comparison_plot, width = 12, height = 8)

# Step 10: Scatter plots comparing relationships
# Create scatter plots with group comparison
scatter_comparison_plots <- list()

for (i in seq_along(wearable_vars)) {
  var <- wearable_vars[i]
  var_label <- case_when(
    var == "mean_hr" ~ "Mean Heart Rate",
    var == "mean_rhr" ~ "Mean Resting Heart Rate",
    var == "mean_bo" ~ "Mean Blood Oxygen",
    var == "total_steps" ~ "Mean Total Steps",
    var == "total_sleep" ~ "Mean Total Sleep"
  )
  
  # Get correlations for both groups
  diabetic_r <- diabetic_correlations$Correlation[diabetic_correlations$Variable == var]
  diabetic_p <- diabetic_correlations$P_value[diabetic_correlations$Variable == var]
  nondiabetic_r <- nondiabetic_correlations$Correlation[nondiabetic_correlations$Variable == var]
  nondiabetic_p <- nondiabetic_correlations$P_value[nondiabetic_correlations$Variable == var]
  
  p <- ggplot(data_imputed, aes_string(x = var, y = "pre_vision", color = "factor(dm_2)")) +
    # 添加 group 参数确保平滑线颜色与点颜色一致
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(color = factor(dm_2), linetype = factor(dm_2), group = factor(dm_2))) +
    # 使用硬编码的颜色值
    scale_color_manual(values = c("0" = "#3182bd",  # Cataract (蓝色)
                                  "1" = "#e6550d"),  # Diabetes (橙色)
                       labels = c("Cataract", "DR"),
                       name = "Group") +
    scale_linetype_manual(values = c("0" = "dashed", "1" = "solid"),
                          labels = c("Cataract", "DR"),
                          name = "Group") +
    labs(
      title = paste(var_label, "vs. Pre-surgery Vision"),
      subtitle = sprintf("DR: r = %.3f, p = %.4f; Cataract: r = %.3f, p = %.4f", 
                         diabetic_r, diabetic_p, nondiabetic_r, nondiabetic_p),
      x = var_label,
      y = "Pre-surgery Vision"
    ) +
    # MODIFIED: Use theme_custom from 100_tools.R
    theme_custom
  
  scatter_comparison_plots[[var]] <- p
  print(p)
  ggsave(paste0("comparison_scatter_", var, "_vs_vision.pdf"), p, width = 8, height = 6)
}

# Step 10 (improved version): Multivariate regression scatter plots
# Calculate R² values for models
diabetic_r2 <- summary(diabetic_model)$r.squared
nondiabetic_r2 <- summary(nondiabetic_model)$r.squared

# Create prediction data
diabetic_data$predicted_vision <- predict(diabetic_model, newdata = diabetic_data)
nondiabetic_data$predicted_vision <- predict(nondiabetic_model, newdata = nondiabetic_data)

# Combine data for plotting
data_with_predictions <- rbind(
  diabetic_data %>% 
    dplyr::select(subject_id, pre_vision, predicted_vision, mean_hr, mean_rhr, mean_bo, 
                  total_steps, total_sleep, dm_2) %>%
    mutate(group = "DR"),
  nondiabetic_data %>% 
    dplyr::select(subject_id, pre_vision, predicted_vision, mean_hr, mean_rhr, mean_bo, 
                  total_steps, total_sleep, dm_2) %>%
    mutate(group = "Cataract")
)

# Create improved scatter plots
improved_scatter_plots <- list()

for (i in seq_along(wearable_vars)) {
  var <- wearable_vars[i]
  var_label <- case_when(
    var == "mean_hr" ~ "Mean Heart Rate",
    var == "mean_rhr" ~ "Mean Resting Heart Rate",
    var == "mean_bo" ~ "Mean Blood Oxygen",
    var == "total_steps" ~ "Mean Total Steps",
    var == "total_sleep" ~ "Mean Total Sleep"
  )
  
  # Get regression coefficients
  diabetic_coef_val <- diabetic_coef_data$Coefficient[diabetic_coef_data$Variable == var]
  nondiabetic_coef_val <- nondiabetic_coef_data$Coefficient[nondiabetic_coef_data$Variable == var]
  
  # 为每个组分配不同颜色
  p <- ggplot(data_with_predictions, aes_string(x = var, y = "pre_vision", color = "group", group = "group")) +
    # 原始数据点
    geom_point(alpha = 0.7) +
    # 多元回归拟合线
    geom_line(aes(y = predicted_vision), linewidth = 1) +
    # 使用硬编码的颜色
    scale_color_manual(values = c("DR" = "#e6550d", 
                                  "Cataract" = "#3182bd"), 
                       name = "Group") +
    labs(
      title = paste(var_label, "vs. Pre-surgery Vision (Multivariate Model)"),
      subtitle = sprintf("DR coef: %.3f, Cataract coef: %.3f", 
                         diabetic_coef_val, nondiabetic_coef_val),
      x = var_label,
      y = "Pre-surgery Vision"
    ) +
    # Add R² information
    annotate("text", x = min(data_with_predictions[[var]], na.rm = TRUE), 
             y = max(data_with_predictions$pre_vision, na.rm = TRUE),
             label = sprintf("DR model R² = %.3f\nCataract model R² = %.3f", 
                             diabetic_r2, nondiabetic_r2),
             hjust = 0, vjust = 1) +
    # MODIFIED: Use theme_custom from 100_tools.R
    theme_custom
  
  improved_scatter_plots[[var]] <- p
  print(p)
  ggsave(paste0("multivariate_scatter_", var, "_vs_vision.pdf"), p, width = 8, height = 6)
}

# Create partial residual plots if car package is available
if(requireNamespace("car", quietly = TRUE)) {
  library(car)
  
  for (i in seq_along(wearable_vars)) {
    var <- wearable_vars[i]
    var_label <- case_when(
      var == "mean_hr" ~ "Mean Heart Rate",
      var == "mean_rhr" ~ "Mean Resting Heart Rate",
      var == "mean_bo" ~ "Mean Blood Oxygen",
      var == "total_steps" ~ "Mean Total Steps",
      var == "total_sleep" ~ "Mean Total Sleep"
    )
    
    # Diabetic group partial residuals
    pdf(paste0("partial_residual_diabetic_", var, ".pdf"), width = 8, height = 6)
    avPlot(diabetic_model, variable = var, 
           main = paste("Partial Residual Plot for", var_label, "- DR Group"))
    dev.off()
    
    # Non-diabetic group partial residuals
    pdf(paste0("partial_residual_nondiabetic_", var, ".pdf"), width = 8, height = 6)
    avPlot(nondiabetic_model, variable = var,
           main = paste("Partial Residual Plot for", var_label, "- Cataract Group"))
    dev.off()
  }
}

# 创建实际值vs预测值的图表
pred_vs_actual <- ggplot(data_with_predictions, aes(x = predicted_vision, y = pre_vision, color = group)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # 使用硬编码的颜色
  scale_color_manual(values = c("DR" = "#e6550d", 
                                "Cataract" = "#3182bd"), 
                     name = "Group") +
  labs(
    title = "Predicted vs. Actual Pre-surgery Vision",
    subtitle = sprintf("DR model R² = %.3f, Cataract model R² = %.3f", 
                       diabetic_r2, nondiabetic_r2),
    x = "Predicted Vision",
    y = "Actual Vision"
  ) +
  # MODIFIED: Use theme_custom from 100_tools.R
  theme_custom

print(pred_vs_actual)
ggsave("predicted_vs_actual_vision.pdf", pred_vs_actual, width = 8, height = 8)

# Step 11: Model performance comparison
# Calculate performance metrics for both models
diabetic_rmse <- sqrt(mean(summary(diabetic_model)$residuals^2))
nondiabetic_rmse <- sqrt(mean(summary(nondiabetic_model)$residuals^2))

# Create performance comparison table
performance_comparison <- data.frame(
  Group = c("DR", "Cataract"),
  R_squared = c(diabetic_r2, nondiabetic_r2),
  Adjusted_R_squared = c(summary(diabetic_model)$adj.r.squared, summary(nondiabetic_model)$adj.r.squared),
  RMSE = c(diabetic_rmse, nondiabetic_rmse)
)

# Print performance comparison
print(performance_comparison)