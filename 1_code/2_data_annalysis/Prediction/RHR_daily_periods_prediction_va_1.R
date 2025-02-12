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
octa<- read.csv("2_data/analysis_data/octa_data.csv")

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
    post_vision = case_when(
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
octa_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa, by = c("ID" = "id")) 

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
patient_list <- split(octa_features, octa_features$ID)
processed_data <- purrr::map(patient_list, process_patient)

# 合并结果
octa_features <- bind_rows(processed_data)

# 验证结果
print(ncol(octa_features))
print(head(octa_features))
print(names(octa_features))


#Extract RHR features from day -4 (4 days before surgery), only using RHR_1
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
  )

sleep_features <- daily_sleep_result %>%
  dplyr::select(
    subject_id,
    deep_sleep = matches("day_-1_deep_sleep$"),
    light_sleep = matches("day_-1_light_sleep$"),
    total_sleep = matches("day_-1_total_sleep$"),
    dream_sleep = matches("day_-1_dream_sleep$"),
    awake = matches("day_-1_awake$"),
    daytime_sleep = matches("day_-1_daytime_sleep$"),
    median_sleep_duration=matches("day_-1_median_sleep_duration$"),
    min_sleep_duration=matches("day_-1_min_sleep_duration$"),
    max_sleep_duration=matches("day_-1_max_sleep_duration$"))

steps_features <- daily_steps_result %>%
  dplyr::select(
    subject_id,
    steps_total = matches("day_-1_steps_total$"),
    steps_mean = matches("day_-1_steps_mean $"),
    steps_median = matches("day_-1_steps_median $"),
    steps_max = matches("day_-1_steps_max$"))

#Combine RHR features with vision data
model_data <- rhr_features %>%
  inner_join(disease_data, by = c("subject_id" = "ID")) %>%
  inner_join(vision_data, by = c("subject_id" = "ID")) %>%
  inner_join(bo_features, by = c("subject_id" = "subject_id")) %>%
  inner_join(sleep_features, by = c("subject_id" = "subject_id")) %>%
  inner_join(steps_features, by = c("subject_id" = "subject_id"))%>%
  inner_join(octa_features, by = c("subject_id" = "ID"))
  # # Remove rows with missing values
  # filter(!is.na(vision_improvement),
  #        !is.na(rhr_mean_1),
  #        !is.na(rhr_min_1),
  #        !is.na(rhr_max_1),
  #        !is.na(rhr_sd_1),
  #        !is.na(rhr_cv_1),
  #        !is.na(rhr_iqr_1),
  #        !is.na(age),
  #        !is.na(gender),
  #        !is.na(bo_mean_1),
  #        !is.na(bo_min_1),
  #        !is.na(bo_max_1),
  #        !is.na(bo_sd_1),
  #        !is.na(bo_cv_1),
  #        !is.na(bo_iqr_1),
  #        !is.na(age),
  #        !is.na(gender))



# Convert gender to factor
model_data$gender <- as.factor(model_data$gender)
model_data$cataract_2 <- as.factor(model_data$cataract_2)
model_data$dm_2  <- as.factor(model_data$dm_2 )
model_data$hypertension_2  <- as.factor(model_data$hypertension_2 )

print(summary(model_data))


######OCTA pca
# 准备OCTA变量
retina_vars <- c(
  # SVD related
  "SVD_NerveFiber_0_6_T0",
  "SVD_Superficial_0_6_T0",
  
  # PA related
  "PA_Avascular_0_6_T0",
  "PA_Choriocapillaris_0_6_T0",
  "PA_Choroid_0_6_T0",
  "PA_DCP_0_6_T0",
  "PA_Deep_0_6_T0",
  "PA_ICP_0_6_T0",
  "PA_InnerRetina_0_6_T0",
  "PA_NerveFiber_0_6_T0",
  # "PA_OuterRetina_0_6_T0",
  "PA_PED_0_6_T0",
  "PA_Retina_0_6_T0",
  "PA_Superficial_0_6_T0",
  "PA_SVP_0_6_T0",
  "PA_Vitreous_0_6_T0",
  
  # VD related
  "VD_DCP_0_6_T0",
  "VD_Deep_0_6_T0",
  "VD_ICP_0_6_T0",
  "VD_InnerRetina_0_6_T0",
  "VD_NerveFiber_0_6_T0",
  "VD_Superficial_0_6_T0",
  "VD_SVP_0_6_T0"
)

# 使用 model_data 中的数据进行PCA
retina_data <- model_data[, retina_vars]

# # 标准化数据（保留NA值）
# retina_scaled <- scale(retina_data)

# 检查相关性
cor_matrix <- cor(retina_data)
print("Correlation matrix of retinal features:")
print(round(cor_matrix[1:5, 1:5], 2))

# 最小-最大标准化函数
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# 对每列应用最小-最大标准化
retina_scaled <- as.data.frame(lapply(retina_data, min_max_scale))

# 然后进行插补
imp <- mice(retina_scaled, m=5)
complete_data <- complete(imp, 1)
pca_result <- prcomp(complete_data, scale. = FALSE)  # 因为已经标准化过了
# 使用第一个插补数据集
complete_data <- complete(imp, 1)
# 然后进行PCA
pca_result <- prcomp(complete_data, scale. = TRUE)

# 计算每个主成分解释的方差
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cum_var_explained <- cumsum(var_explained)

# 创建碎石图数据
scree_data <- data.frame(
  PC = 1:length(var_explained),
  Variance = var_explained,
  Cumulative = cum_var_explained
)

# 绘制碎石图
library(ggplot2)
p1 <- ggplot(scree_data, aes(x = PC)) +
  geom_line(aes(y = Variance), color = "blue") +
  geom_point(aes(y = Variance), color = "blue") +
  geom_line(aes(y = Cumulative), color = "red") +
  geom_point(aes(y = Cumulative), color = "red") +
  labs(
    title = "Scree Plot of PCA",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  ) +
  theme_bw() +
  scale_y_continuous(
    name = "Individual Variance Explained",
    sec.axis = sec_axis(~., name = "Cumulative Variance Explained")
  )

print(p1)

# 确定需要保留的主成分数（解释80%的方差）
n_components <- which(cum_var_explained >= 0.8)[1]
print(paste("\nNumber of components needed to explain 80% of variance:", n_components))

# 获取特征对主成分的贡献
loadings <- pca_result$rotation[, 1:min(5, n_components)]

# 提取主成分得分
pc_scores <- pca_result$x[, 1:n_components]
colnames(pc_scores) <- paste0("PC", 1:n_components)

# 创建载荷热图数据
loadings_df <- as.data.frame(loadings)
loadings_df$Feature <- rownames(loadings_df)
loadings_long <- tidyr::pivot_longer(loadings_df, 
                                     cols = -Feature, 
                                     names_to = "PC", 
                                     values_to = "Loading")

# 创建热图
p2 <- ggplot(loadings_long, aes(x = PC, y = Feature, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(title = "Feature Loadings on Principal Components",
       x = "Principal Component",
       y = "Feature")

print(p2)


######lasso
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
    deep_sleep, light_sleep, total_sleep, dream_sleep, awake,
    daytime_sleep, median_sleep_duration, min_sleep_duration,
    max_sleep_duration, 
    # Steps features
    starts_with("steps_"),
    # Demographics and medical history
    age, gender, cataract_2, dm_2, hypertension_2,pre_vision,
    # # OCTA features
    retina_vars,
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

# Initialize the mice imputation
library(mice)
ini <- mice(features_for_imputation, maxit = 0, printFlag = FALSE)

# Set all variables to use random forest imputation
meth <- rep("rf", ncol(features_for_imputation))
names(meth) <- colnames(features_for_imputation)

# Perform random forest imputation using mice
set.seed(123)  # for reproducibility
imp <- mice(features_for_imputation, 
            method = meth,
            m = 5,
            maxit = 10,
            printFlag = FALSE)

# Use first imputed dataset
complete_data <- complete(imp, 1)

# Print missing values summary after imputation
cat("\nMissing values after imputation:\n")
print(colSums(is.na(complete_data)))


# ######impute using rf
# # 加载必要的包
# library(randomForest)
# library(dplyr)
# 
# # 假设你的数据框名为 df
# # 首先将数据分成两部分：有缺失值的行和没有缺失值的行
# df_missing <- complete_data[is.na(complete_data$PA_OuterRetina_0_6_T0), ]
# df_complete <- complete_data[!is.na(complete_data$PA_OuterRetina_0_6_T0), ]
# 
# # 选择预测变量（除了目标变量PA_OuterRetina_0_6_T0之外的所有相关列）
# # 这里需要你根据实际情况调整预测变量的选择
# predictors <- names(complete_data)[!names(complete_data) %in% c("PA_OuterRetina_0_6_T0")]
# 
# # 构建随机森林模型
# rf_model <- randomForest(
#   formula = as.formula(paste("PA_OuterRetina_0_6_T0 ~", 
#                              paste(predictors, collapse = " + "))),
#   data = df_complete,
#   ntree = 500,
#   mtry = floor(sqrt(length(predictors))),
#   importance = TRUE
# )
# 
# # 使用模型预测缺失值
# predicted_values <- predict(rf_model, df_missing)
# 
# # 将预测值填回原始数据框
# complete_data$PA_OuterRetina_0_6_T0[is.na(complete_data$PA_OuterRetina_0_6_T0)] <- predicted_values




###############octa pca
# Get OCTA variables for PCA
retina_vars <- grep("_0_6_T0$", colnames(complete_data), value = TRUE)

# Extract OCTA data from imputed dataset
retina_data <- complete_data[, retina_vars]

# Handle infinite values
retina_data[is.infinite(as.matrix(retina_data))] <- NA

# Check which columns have NAs
na_cols <- colSums(is.na(retina_data))
cat("Columns with NA values:\n")
print(na_cols[na_cols > 0])

# Create predictor matrix for mice
pred_matrix <- matrix(1, ncol = ncol(retina_data), nrow = ncol(retina_data))
diag(pred_matrix) <- 0
colnames(pred_matrix) <- rownames(pred_matrix) <- colnames(retina_data)

# Perform imputation with more iterations and RF method
imp_retina <- mice(retina_data, 
                   m = 5,
                   maxit = 10,  # Increased iterations
                   method = 'rf',
                   predictorMatrix = pred_matrix,
                   printFlag = FALSE)

# Get complete imputed dataset
retina_data <- complete(imp_retina, 1)

# Verify no missing values
cat("\nVerifying no missing values after imputation:\n")
print(colSums(is.na(retina_data)))
cat("\nVerifying no infinite values after imputation:\n")
print(colSums(is.infinite(as.matrix(retina_data))))

# Standardize the data
retina_scaled <- scale(retina_data)

# Double check for any remaining issues
if(any(is.na(retina_scaled)) || any(is.infinite(retina_scaled))) {
  cat("Warning: Still have NA or infinite values after scaling\n")
  # Replace any remaining problematic values with column means
  for(col in 1:ncol(retina_scaled)) {
    bad_vals <- is.na(retina_scaled[,col]) | is.infinite(retina_scaled[,col])
    if(any(bad_vals)) {
      retina_scaled[bad_vals,col] <- mean(retina_scaled[!bad_vals,col])
    }
  }
}

# Perform PCA
pca_result <- prcomp(retina_scaled, scale. = FALSE)  # Already scaled

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
n_components <- which(cumsum(var_explained) >= 0.8)[1]

# Extract PC scores
pc_scores <- pca_result$x[, 1:n_components]
colnames(pc_scores) <- paste0("PC", 1:n_components)

# Create final dataset
model_data_final <- cbind(
  dplyr::select(complete_data, -all_of(retina_vars)),
  pc_scores
)

# Verify dimensions
cat("\nFinal dimensions:\n")
print(dim(model_data_final))

# Verify dimensions
cat("\nDimensions of final dataset:\n")
print(dim(model_data_final))


###########################
# Select features for LASSO
feature_cols <- setdiff(colnames(model_data_final), "vision_improvement")

library(glmnet)

# 准备数据
x <- model_data_final[, feature_cols]
y <- model_data_final$vision_improvement

model_lasso<- glmnet(x, y, nlambda=100, alpha=1)
print(model_lasso)
plot(model_lasso, xvar="lambda", label=TRUE)
plot(model_lasso, xvar="norm", label=TRUE) #或xvar="dev"

#交叉验证，参数解释见下文
# 创建模型矩阵，这会自动处理因子型变量
x <- model.matrix(~., data=model_data_final[, feature_cols])
# 删除截距列（第一列）
x <- x[, -1]

# 现在确保 y 是数值型
y <- as.matrix(model_data_final$vision_improvement)

# 运行交叉验证
set.seed(1234)
alpha1.fit.cv <- cv.glmnet(x, y, type.measure="deviance", alpha=1)
plot(alpha1.fit.cv)
print(alpha1.fit.cv)

coef(alpha1.fit.cv,s=alpha1.fit.cv$lambda.min)

#提取特征
feature_all<-as.data.frame(as.matrix(coef(alpha1.fit.cv,s=alpha1.fit.cv$lambda.min)))
colnames(feature_all)<-"coff"
feature_opt<-feature_all%>%filter(abs(coff)>0)
rownames(feature_opt)

# 加载必要的包
library(caret)
library(randomForest)
library(xgboost)
library(tidyverse)

# 选择LASSO筛选出的特征
selected_features <- c(
  "rhr_sd_1", "rhr_max_1", "rhr_median_1", 
  "deep_sleep", "pre_vision",
  "median_sleep_duration", "steps_total", 
  "dm_2", "PC5", "PC2", "PC4"
)

# 准备建模数据
model_data <- model_data_final[, c(selected_features, "vision_improvement")]

# 设置随机种子
set.seed(123)

# 定义5折交叉验证
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  summaryFunction = defaultSummary
)

# 1. 随机森林模型
rf_model <- train(
  vision_improvement ~ .,
  data = model_data,
  method = "rf",
  trControl = ctrl,
  metric = "RMSE",
  importance = TRUE
)

# 2. XGBoost模型
xgb_model <- train(
  vision_improvement ~ .,
  data = model_data,
  method = "xgbTree",
  trControl = ctrl,
  metric = "RMSE",
  tuneLength = 5
)

# 获取每个模型的预测结果
rf_predictions <- rf_model$pred
xgb_predictions <- xgb_model$pred

# 定义计算性能指标的函数
calculate_metrics <- function(predictions) {
  predictions %>%
    group_by(Resample) %>%
    summarise(
      RMSE = sqrt(mean((obs - pred)^2)),
      MAE = mean(abs(obs - pred)),
      R2 = cor(obs, pred)^2
    )
}

# 计算每个模型的性能指标
rf_performance <- calculate_metrics(rf_predictions)
xgb_performance <- calculate_metrics(xgb_predictions)

# 计算整体性能
calculate_overall_metrics <- function(fold_metrics) {
  summarise(fold_metrics,
            RMSE_mean = mean(RMSE),
            RMSE_sd = sd(RMSE),
            MAE_mean = mean(MAE),
            MAE_sd = sd(MAE),
            R2_mean = mean(R2),
            R2_sd = sd(R2)
  )
}

# 打印随机森林结果
print("\nRandom Forest Results:")
print(rf_model)
print("\nRandom Forest Per-fold Performance:")
print(rf_performance)
print("\nRandom Forest Overall Performance:")
print(calculate_overall_metrics(rf_performance))

# 打印XGBoost结果
print("\nXGBoost Results:")
print(xgb_model)
print("\nXGBoost Per-fold Performance:")
print(xgb_performance)
print("\nXGBoost Overall Performance:")
print(calculate_overall_metrics(xgb_performance))

# 随机森林特征重要性
print("\nRandom Forest Feature Importance:")
print(varImp(rf_model))

# 可视化预测结果比较
predictions_combined <- bind_rows(
  rf_predictions %>% mutate(Model = "Random Forest"),
  xgb_predictions %>% mutate(Model = "XGBoost")
)

ggplot(predictions_combined, aes(x = obs, y = pred, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  facet_wrap(~Model) +
  labs(
    x = "Observed Vision Improvement",
    y = "Predicted Vision Improvement",
    title = "Predicted vs Observed Values: Model Comparison"
  ) +
  theme_minimal()

# 箱线图比较模型性能
bind_rows(
  rf_performance %>% mutate(Model = "Random Forest"),
  xgb_performance %>% mutate(Model = "XGBoost")
) %>%
  ggplot(aes(x = Model, y = R2)) +
  geom_boxplot() +
  labs(
    title = "Model Performance Comparison",
    y = "R-squared"
  ) +
  theme_minimal()






# Verify these columns exist in the dataset
existing_cols <- selected_features[selected_features %in% colnames(model_data_final)]
missing_cols <- selected_features[!selected_features %in% colnames(model_data_final)]

cat("\nExisting columns:\n")
print(existing_cols)
cat("\nMissing columns:\n")
print(missing_cols)

# Prepare modeling data with correct features
model_data <- model_data_final[, c(selected_features, "vision_improved_factor")]

# 如果已经有vision_improved_factor，直接重编码
if("vision_improved_factor" %in% colnames(model_data)) {
  levels(model_data$vision_improved_factor) <- c("0", "1")
} else {
  # 如果不存在，则基于vision_improvement创建
  model_data$vision_improved_factor <- factor(
    ifelse(model_data$vision_improvement >= 0, "1", "0"),
    levels = c("0", "1")
  )
}

# 验证新的因子水平
print(levels(model_data$vision_improved_factor))

# Verify the structure of the modeling data
cat("\nStructure of modeling data:\n")
str(model_data)

# Convert factor levels to valid R variable names
levels(model_data$vision_improved_factor) <- c("Class0", "Class1")

# Verify the new levels
print(levels(model_data$vision_improved_factor))

# Set random seed
set.seed(123)

# Define cross-validation with additional metrics for classification
ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  classProbs = TRUE,  # Required for ROC curves
  summaryFunction = twoClassSummary  # Use classification metrics
)

# 1. Random Forest Classification Model
rf_model <- train(
  vision_improved_factor ~ .,
  data = model_data,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",  # Use AUC-ROC as the metric
  ntree = 500
)

# 2. XGBoost Classification Model
xgb_model <- train(
  vision_improved_factor ~ .,
  data = model_data,
  method = "xgbTree",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5
)

# Calculate ROC curves and plot
library(pROC)

# Function to calculate and plot ROC curves
plot_roc_curves <- function(rf_model, xgb_model) {
  # Get predictions for both models
  rf_probs <- rf_model$pred %>%
    mutate(Model = "Random Forest")
  
  xgb_probs <- xgb_model$pred %>%
    mutate(Model = "XGBoost")
  
  # Combine predictions
  all_probs <- bind_rows(rf_probs, xgb_probs)
  
  # Calculate ROC curves
  roc_curves <- all_probs %>%
    group_by(Model) %>%
    group_modify(~{
      roc_obj <- roc(response = .x$obs, 
                     predictor = .x$Improved,
                     levels = c("NoImprovement", "Improved"))
      data.frame(
        sensitivity = roc_obj$sensitivities,
        specificity = roc_obj$specificities,
        auc = auc(roc_obj)
      )
    })
  
  # Create ROC plot
  p <- ggplot(roc_curves, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
    geom_line(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = "ROC Curves for Vision Improvement Classification",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    annotate("text", x = 0.75, y = 0.25, 
             label = sprintf("RF AUC: %.3f\nXGB AUC: %.3f", 
                             unique(roc_curves$auc[roc_curves$Model == "Random Forest"]),
                             unique(roc_curves$auc[roc_curves$Model == "XGBoost"]))) +
    theme_minimal()
  
  print(p)
  
  # Return AUC values
  return(roc_curves %>% 
           group_by(Model) %>% 
           summarise(AUC = unique(auc)))
}

# Calculate performance metrics for both models
get_model_metrics <- function(model) {
  pred_data <- model$pred
  
  # Calculate metrics per fold
  fold_metrics <- pred_data %>%
    group_by(Resample) %>%
    summarise(
      AUC = as.numeric(roc(obs, Improved)$auc), # Extract numeric AUC value
      Accuracy = mean(pred == obs),
      Sensitivity = sensitivity(factor(pred), factor(obs), positive = "Improved"),
      Specificity = specificity(factor(pred), factor(obs), negative = "NoImprovement")
    )
  
  # Calculate mean and SD across folds
  overall_metrics <- fold_metrics %>%
    summarise(across(everything(), 
                     list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"))
  
  return(list(fold_metrics = fold_metrics, overall_metrics = overall_metrics))
}

# Calculate performance metrics for both models
get_model_metrics <- function(model) {
  pred_data <- model$pred
  
  # Calculate metrics per fold
  fold_metrics <- pred_data %>%
    group_by(Resample) %>%
    summarise(
      AUC = as.numeric(roc(obs, Class1)$auc),  # Changed from Improved to Class1
      Accuracy = mean(pred == obs),
      Sensitivity = sensitivity(factor(pred), factor(obs), positive = "Class1"),  # Updated positive class
      Specificity = specificity(factor(pred), factor(obs), negative = "Class0")  # Updated negative class
    )
  
  # Calculate mean and SD across folds
  overall_metrics <- fold_metrics %>%
    summarise(across(everything(), 
                     list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"))
  
  return(list(fold_metrics = fold_metrics, overall_metrics = overall_metrics))
}

# Function to calculate and plot ROC curves
plot_roc_curves <- function(rf_model, xgb_model) {
  # Get predictions for both models
  rf_probs <- rf_model$pred %>%
    mutate(Model = "Random Forest")
  
  xgb_probs <- xgb_model$pred %>%
    mutate(Model = "XGBoost")
  
  # Combine predictions
  all_probs <- bind_rows(rf_probs, xgb_probs)
  
  # Calculate ROC curves
  roc_curves <- all_probs %>%
    group_by(Model) %>%
    group_modify(~{
      roc_obj <- roc(response = .x$obs, 
                     predictor = .x$Class1,  # Changed from Improved to Class1
                     levels = c("Class0", "Class1"))  # Updated levels
      data.frame(
        sensitivity = roc_obj$sensitivities,
        specificity = roc_obj$specificities,
        auc = as.numeric(auc(roc_obj))
      )
    })
  
  # Create ROC plot
  p <- ggplot(roc_curves, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
    geom_line(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = "ROC Curves for Vision Improvement Classification",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    annotate("text", x = 0.75, y = 0.25, 
             label = sprintf("RF AUC: %.3f\nXGB AUC: %.3f", 
                             unique(roc_curves$auc[roc_curves$Model == "Random Forest"]),
                             unique(roc_curves$auc[roc_curves$Model == "XGBoost"]))) +
    theme_minimal()
  
  print(p)
  
  # Return AUC values
  return(roc_curves %>% 
           group_by(Model) %>% 
           summarise(AUC = unique(auc)))
}

# Calculate and print metrics
rf_metrics <- get_model_metrics(rf_model)
xgb_metrics <- get_model_metrics(xgb_model)

# Print results
cat("\nRandom Forest Results:\n")
print(rf_metrics$overall_metrics)

cat("\nXGBoost Results:\n")
print(xgb_metrics$overall_metrics)

# Plot ROC curves and get AUC values
auc_values <- plot_roc_curves(rf_model, xgb_model)
print(auc_values)

# Plot feature importance for Random Forest
rf_importance <- varImp(rf_model)$importance %>%
  rownames_to_column("Feature") %>%
  arrange(desc(Overall))

importance_plot <- ggplot(rf_importance, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Random Forest Feature Importance",
    x = "Features",
    y = "Importance Score"
  ) +
  theme_minimal()

print(importance_plot)

# Print confusion matrices
cat("\nRandom Forest Confusion Matrix:\n")
confusionMatrix(rf_model)

cat("\nXGBoost Confusion Matrix:\n")
confusionMatrix(xgb_model)
