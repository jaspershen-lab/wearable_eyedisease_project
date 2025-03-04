library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)

# Load the RHR data
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/time_period_rhr_results.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/daily_bo_result.rda")
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/daily_sleep_result.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result.rda")


# Read baseline info
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow<- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness<- read.csv("2_data/analysis_data/octa_data_thickness.csv")

####disease data
disease_data <- baseline_info %>%
  mutate(
    cataract_2 = case_when(
      cataract == 1 ~ 0,  
      cataract %in% c(2, 3, 4) ~ 1, 
      TRUE ~ NA_real_
    ),
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  
      diabetes_history == 2 ~ 0, 
      TRUE ~ NA_real_
    ),
    hypertension_2 = case_when(
      hypertension_history == 1 ~ 1,  
      hypertension_history == 2 ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(ID, cataract_2, dm_2, hypertension_2)

# Process vision data
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,  # 右眼手术
      surgery_eye_1 == 1 ~ os_corrected_bas,  # 左眼手术
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # 双眼手术取平均
      TRUE ~ NA_real_
    ),
    post_vision= case_when(
      surgery_eye_1 == 0 ~ od_corrected_1w,   # 右眼手术后
      surgery_eye_1 == 1 ~ os_corrected_1w,   # 左眼手术后
      surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # 双眼手术后取平均
      TRUE ~ NA_real_
    ),
    vision_improvement = post_vision - pre_vision,
    vision_improved = if_else(vision_improvement >= 0, 1, 0),  # 添加二分类指标
    vision_improved_factor = factor(vision_improved, 
                                    levels = c(0, 1), 
                                    labels = c("NoImprovement", "Improved"))  # 添加因子版本
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision, 
                vision_improvement, vision_improved, vision_improved_factor,
                age, gender)



# OCTA features
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# 来处理单个患者的数据
process_patient <- function(patient_data) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术选择左眼数据，右眼和双眼手术选择右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"

  # 选择列
  cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
  cols_to_keep <- cols_to_keep[grep("T0$", cols_to_keep)]

  # 选择数据并重命名
  result <- patient_data %>%
    dplyr::select("ID", all_of(cols_to_keep)) %>%
    rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)

  return(result)
}

# 对每个患者分别处理
patient_list <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data <- purrr::map(patient_list, process_patient)

# 合并结果
octa_bloodflow_features <- bind_rows(processed_data)

# 验证结果
print(ncol(octa_bloodflow_features))
print(head(octa_bloodflow_features))
print(names(octa_bloodflow_features))

bloodflow_var <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # 保留ID列
    matches(".*_0_6_T0$"),  # 选择所有以 0_6_T0 结尾的列
    -"PA_OuterRetina_0_6_T0",  # 直接用列名排除
    -"PA_PED_0_6_T0"
  )

# 结果包含哪些变量
colnames(bloodflow_var)


octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# 来处理单个患者的数据
process_patient <- function(patient_data) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术选择左眼数据，右眼和双眼手术选择右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  # 选择列
  cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
  cols_to_keep <- cols_to_keep[grep("T0$", cols_to_keep)]
  
  # 选择数据并重命名
  result <- patient_data %>%
    dplyr::select("ID", all_of(cols_to_keep)) %>%
    rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
  
  return(result)
}

# 对每个患者分别处理
patient_list <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data <- purrr::map(patient_list, process_patient)

# 合并结果
octa_thickness_features <- bind_rows(processed_data)

# 验证结果
print(ncol(octa_thickness_features))
print(head(octa_thickness_features))
print(names(octa_thickness_features))

thickness_var <- octa_thickness_features %>%
  dplyr::select(
    ID,  # 保留ID列
    matches(".*_0_6_T0$"),
    "Thickness_PED_0_6_T0"
  )

# 结果包含哪些变量
colnames(thickness_var)



#Extract RHR features from day -1 (1 days before surgery), only using RHR_1
rhr_features <- daily_rhr_result %>%
  dplyr::select(
    subject_id,
    # Mean RHR
    rhr_mean_1 = matches("day_-1_mean_rhr_1$"),
    # Minimum RHR
    rhr_min_1 = matches("day_-1_min_rhr_1$"),
    # Maximum RHR
    rhr_max_1 = matches("day_-1_max_rhr_1$"),
    # Standard deviation
    rhr_sd_1 = matches("day_-1_sd_rhr_1$"),
    # Coefficient of variation
    rhr_median_1 = matches("day_-1_median_rhr_1$"),
    # Interquartile range
    rhr_iqr_1 = matches("day_-1_iqr_rhr_1$")
    # rhr_skew_1 = matches("day_-1_skew_rhr_1$"),
    # rhr_kurt_1 = matches("day_-1_kurt_rhr_1$")
  )


bo_features <- daily_bo_result %>%
  dplyr::select(
    subject_id,
    # Mean RHR
    bo_mean = matches("day_-1_mean_bo$"),
    # Minimum RHR
    bo_min = matches("day_-1_min_bo$"),
    # Maximum RHR
    bo_max = matches("day_-1_max_bo$"),
    # Standard deviation
    bo_sd = matches("day_-1_sd_bo$"),
    # Coefficient of variation
    bo_median = matches("day_-1_median_bo$"),
    # Interquartile range
    bo_iqr = matches("day_-1_iqr_bo$")
    # bo_skew = matches("day_-1_skew_bo$"),
    # bo_kurt = matches("day_-1_kurt_bo$")
  )

sleep_features <- daily_sleep_result %>%
  dplyr::select(
    subject_id,
    deep_sleep = matches("day_-1_deep_sleep$"),
    light_sleep = matches("day_-1_light_sleep$"),
    total_sleep = matches("day_-1_total_sleep$"),
    dream_sleep = matches("day_-1_dream_sleep$"),
    awake = matches("day_-1_awake$"),
    daytime_sleep = matches("day_-1_daytime_sleep$"))
    

steps_features <- daily_steps_result %>%
  dplyr::select(
    subject_id,
    steps_total = matches("day_-1_steps_total$"),
    steps_mean = matches("day_-1_steps_mean$"),
    steps_median = matches("day_-1_steps_median$"),
    steps_max = matches("day_-1_steps_max$"))

#Combine RHR features with vision data
model_data <- rhr_features %>%
  inner_join(disease_data, by = c("subject_id" = "ID")) %>%
  inner_join(vision_data, by = c("subject_id" = "ID")) %>%
  inner_join(bo_features, by = c("subject_id" = "subject_id")) %>%
  inner_join(sleep_features, by = c("subject_id" = "subject_id")) %>%
  inner_join(steps_features, by = c("subject_id" = "subject_id"))%>%
  inner_join(bloodflow_var, by = c("subject_id" = "ID"))%>%
  inner_join(thickness_var, by = c("subject_id" = "ID"))


# Convert gender to factor
model_data$gender <- as.factor(model_data$gender)
model_data$cataract_2 <- as.factor(model_data$cataract_2)
model_data$dm_2  <- as.factor(model_data$dm_2 )
model_data$hypertension_2  <- as.factor(model_data$hypertension_2 )

print(summary(model_data))



######imputation
# First remove rows where vision_improvement is NA (our target variable)
model_data_clean <- model_data[!is.na(model_data$vision_improvement), ]

# Select features for imputation
features_for_imputation <- model_data_clean %>%
  dplyr::select(
    # RHR features
    starts_with("rhr_"),
    # BO features
    starts_with("bo_"),
    # Sleep features
    deep_sleep, light_sleep, total_sleep, dream_sleep,
    awake, daytime_sleep,
    # Steps features
    starts_with("steps_"),
    # Demographics and medical history
    age, gender, cataract_2, dm_2, hypertension_2, pre_vision,
    # OCTA features
    matches("(SVD|PA|VD).*_0_6_T0"),  # 选择所有 bloodflow 0_6_T0 变量
    matches("Thickness.*_0_6_T0"),     # 选择所有 thickness 0_6_T0 变量
    # Target variable
    vision_improvement
  )

# Print missing values summary before imputation
cat("\nMissing values before imputation:\n")
print(colSums(is.na(features_for_imputation)))

# Ensure categorical variables are properly coded as factors
features_for_imputation$gender <- as.factor(features_for_imputation$gender)
features_for_imputation$cataract_2 <- as.factor(features_for_imputation$cataract_2)
features_for_imputation$dm_2 <- as.factor(features_for_imputation$dm_2)
features_for_imputation$hypertension_2 <- as.factor(features_for_imputation$hypertension_2)


# 创建一个列表来存储wearable数据异常值信息
all_outliers <- list()

# 创建一个函数来识别和处理异常值，但使用更保守的阈值
detect_and_handle_outliers <- function(data, patterns, z_threshold = 3, replace_with_na = TRUE) {
  all_outliers <- list()
  
  for(col in colnames(data)) {
    # 检查列名是否匹配模式
    is_matched <- any(sapply(patterns, function(pattern) grepl(pattern, col)))
    
    if(is_matched && is.numeric(data[[col]])) {
      # 获取非NA值
      valid_values <- data[[col]][!is.na(data[[col]])]
      
      if(length(valid_values) > 5) {  # 确保有足够的数据计算z分数
        # 计算z-scores
        valid_mean <- mean(valid_values, na.rm = TRUE)
        valid_sd <- sd(valid_values, na.rm = TRUE)
        
        if(valid_sd > 0) {  # 避免除以0
          z_scores <- (data[[col]] - valid_mean) / valid_sd
          
          # 找出超过阈值的值
          outliers <- which(abs(z_scores) > z_threshold & !is.na(z_scores))
          
          # 存储和处理异常值信息
          if(length(outliers) > 0) {
            all_outliers[[col]] <- list(
              indices = outliers,
              values = data[[col]][outliers],
              z_scores = z_scores[outliers]
            )
            
            # 根据参数决定是否替换为NA
            if(replace_with_na) {
              data[outliers, col] <- NA
            }
            
            # 打印该变量的异常值信息
            cat("\nVariable:", col, "\n")
            cat("Number of outliers:", length(outliers), "\n")
            cat("Outlier values:", paste(round(data[[col]][outliers], 3), collapse=", "), "\n")
            cat("Z-scores:", paste(round(z_scores[outliers], 3), collapse=", "), "\n")
            
            # 打印处理前后的统计信息
            if(replace_with_na) {
              cat("处理后的", col, "统计信息：\n")
              print(summary(data[[col]]))
            }
          }
        }
      }
    }
  }
  
  return(list(data = data, outliers = all_outliers))
}

# 定义wearable数据的列名模式
wearable_patterns <- c("^rhr_", "^bo_", "^steps_", "^deep_sleep$", "^light_sleep$", 
                       "^total_sleep$", "^dream_sleep$", "^awake$", "^daytime_sleep$")

# 处理wearable数据异常值 (使用标准的3个标准差)
wearable_result <- detect_and_handle_outliers(
  features_for_imputation, 
  wearable_patterns, 
  z_threshold = 3, 
  replace_with_na = TRUE
)
features_for_imputation <- wearable_result$data

# 定义OCTA数据的列名模式 
octa_patterns <- c("(SVD|PA|VD).*_0_6_T0$", "Thickness.*_0_6_T0$")

# 处理OCTA数据异常值 (使用更保守的4.5个标准差)
cat("\n\n===== OCTA数据异常值检测 (使用4.5个标准差作为阈值) =====\n")
octa_result <- detect_and_handle_outliers(
  features_for_imputation, 
  octa_patterns, 
  z_threshold = 4.5,  # 使用更保守的阈值
  replace_with_na = TRUE  # 设置为TRUE表示替换为NA，FALSE表示只检测不替换
)
features_for_imputation <- octa_result$data

# 保存所有异常值信息
all_outliers <- c(wearable_result$outliers, octa_result$outliers)

# 打印OCTA数据异常值摘要
cat("\n\n===== OCTA数据异常值摘要 =====\n")
octa_outlier_counts <- sapply(octa_result$outliers, function(x) length(x$indices))
if(length(octa_outlier_counts) > 0) {
  cat("OCTA变量异常值数量:\n")
  print(octa_outlier_counts)
  cat("\n总共检测到", sum(octa_outlier_counts), "个OCTA数据异常值\n")
} else {
  cat("未检测到OCTA数据异常值\n")
}


# Initialize the mice imputation
library(mice)
ini <- mice(features_for_imputation, maxit = 0, printFlag = FALSE)

# # Set all variables to use random forest imputation
# meth <- rep("rf", ncol(features_for_imputation))
# names(meth) <- colnames(features_for_imputation)

# 设置更详细的方法
meth <- make.method(features_for_imputation)
pred <- make.predictorMatrix(features_for_imputation)

# # Calculate the proportion of missing values for each variable
# missing_prop <- colMeans(is.na(features_for_imputation))
# print(missing_prop)
# 
# # 调整预测矩阵，排除高缺失率的变量作为预测变量
# high_missing <- names(which(missing_prop > 0.3))
# pred[, high_missing] <- 0

# Perform random forest imputation using mice
set.seed(123)  # for reproducibility
imp <- mice(features_for_imputation, 
            method = meth,
            m = 5,
            maxit = 10,
            printFlag = FALSE)

# Use first imputed dataset
complete_data <- complete(imp, 1)

complete_data$vision_improved <- factor(ifelse(complete_data$vision_improvement >= 0, 1, 0),
                                        levels = c(0, 1))

# Print missing values summary after imputation
cat("\nMissing values after imputation:\n")
print(colSums(is.na(complete_data)))

# # 使用 PA_Avascular_0_6_T0 的值来填充 PA_OuterRetina_0_6_T0 的缺失值
# complete_data$PA_OuterRetina_0_6_T0[is.na(complete_data$PA_OuterRetina_0_6_T0)] <- 
#   complete_data$PA_Avascular_0_6_T0[is.na(complete_data$PA_OuterRetina_0_6_T0)]





#####OCTA PCA
# Get bloodflow variables
bloodflow_vars <- grep("(SVD|PA|VD).*_0_6_T0$", colnames(complete_data), value = TRUE)
cat("Bloodflow variables:", "\n")
print(bloodflow_vars)

# Get thickness variables
thickness_vars <- grep("Thickness.*_0_6_T0$", colnames(complete_data), value = TRUE)
cat("\nThickness variables:", "\n")
print(thickness_vars)

# Check for missing values in bloodflow data
bloodflow_data <- complete_data[, bloodflow_vars]
cat("\nMissing values in bloodflow data:", "\n")
print(colSums(is.na(bloodflow_data)))

# Check for missing values in thickness data
thickness_data <- complete_data[, thickness_vars]
cat("\nMissing values in thickness data:", "\n")
print(colSums(is.na(thickness_data)))

# Function to handle data preparation and PCA
perform_pca <- function(data, type) {
  # Remove any columns that are all NA
  valid_cols <- colSums(is.na(data)) < nrow(data)
  if(any(!valid_cols)) {
    cat("\nRemoving columns with all NA in", type, "data:", 
        names(data)[!valid_cols], "\n")
    data <- data[, valid_cols, drop = FALSE]
  }
  
  # Handle any remaining NAs by mean imputation
  for(col in colnames(data)) {
    if(any(is.na(data[[col]]))) {
      col_mean <- mean(data[[col]], na.rm = TRUE)
      data[[col]][is.na(data[[col]])] <- col_mean
    }
  }
  
  # Scale the data
  scaled_data <- scale(data)
  
  # Check for infinite values
  if(any(!is.finite(scaled_data))) {
    cat("\nWarning: Infinite values found in scaled", type, "data\n")
    scaled_data[!is.finite(scaled_data)] <- 0
  }
  
  # Perform PCA
  pca_result <- prcomp(scaled_data, scale. = FALSE)
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  n_components <- which(cumsum(var_explained) >= 0.8)[1]
  
  # Extract PC scores
  pc_scores <- pca_result$x[, 1:n_components]
  if(type == "Bloodflow") {
    colnames(pc_scores) <- paste0("OCTA_Flow_PC", 1:n_components)
  } else {
    colnames(pc_scores) <- paste0("OCTA_Thick_PC", 1:n_components)
  }
  
  return(list(
    pca = pca_result,
    var_explained = var_explained,
    n_components = n_components,
    pc_scores = pc_scores
  ))
}

# Perform PCA for bloodflow
bloodflow_pca <- perform_pca(bloodflow_data, "Bloodflow")

# Perform PCA for thickness
thickness_pca <- perform_pca(thickness_data, "Thickness")

# Visualize variance explained for both PCAs
par(mfrow = c(1, 2))

# Bloodflow PCA plot
plot(cumsum(bloodflow_pca$var_explained),
     type = "b",
     xlab = "Number of Components",
     ylab = "Cumulative Proportion of Variance Explained",
     main = "Bloodflow PCA Variance Explanation")
abline(h = 0.8, col = "red", lty = 2)

# Thickness PCA plot
plot(cumsum(thickness_pca$var_explained),
     type = "b",
     xlab = "Number of Components",
     ylab = "Cumulative Proportion of Variance Explained",
     main = "Thickness PCA Variance Explanation")
abline(h = 0.8, col = "red", lty = 2)

# Create final dataset with both sets of PC scores
model_data_final <- cbind(
  dplyr::select(complete_data, -all_of(c(bloodflow_vars, thickness_vars))),
  bloodflow_pca$pc_scores,
  thickness_pca$pc_scores
)

# Print summary statistics
cat("\nNumber of bloodflow components retained:", bloodflow_pca$n_components)
cat("\nCumulative variance explained by bloodflow components:", 
    cumsum(bloodflow_pca$var_explained)[bloodflow_pca$n_components])
cat("\nNumber of thickness components retained:", thickness_pca$n_components)
cat("\nCumulative variance explained by thickness components:", 
    cumsum(thickness_pca$var_explained)[thickness_pca$n_components])

# Print the PCA column names in final dataset
cat("\n\nPCA columns in final dataset:\n")
pca_cols <- grep("^OCTA_", names(model_data_final), value = TRUE)
print(pca_cols)

# Save the results
dir.create("3_data_analysis/2_data_analysis/Prediction/model_data_final", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/Prediction/model_data_final")

save(complete_data, file = "complete_data.rda", compress = "xz")
save(model_data_final, file = "model_data_final.rda", compress = "xz")

