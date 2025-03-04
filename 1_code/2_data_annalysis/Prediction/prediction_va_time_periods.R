# Load required libraries
library(glmnet)
library(caret)
library(randomForest)
library(xgboost)
library(tidyverse)
library(pROC)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())


# load("3_data_analysis/2_data_analysis/Prediction/model_data_final/model_data_final.rda")
load("3_data_analysis/2_data_analysis/Prediction/model_data_final/complete_data.rda")


#===============================================================================
# 1. LASSO Feature Selection
#===============================================================================
# Select features
feature_cols <- setdiff(colnames(complete_data), 
                        c("vision_improvement", "vision_improved"))

# Prepare data for LASSO
x <- as.matrix(complete_data[, feature_cols])
y <- as.numeric(complete_data$vision_improvement)

# Fit LASSO model
model_lasso <- glmnet(x, y, nlambda=100, alpha=1)
plot(model_lasso, xvar="lambda", label=TRUE)
plot(model_lasso, xvar="norm", label=TRUE) #或xvar="dev"

# Run LASSO with cross-validation
set.seed(1234)
cv_fit <- cv.glmnet(x, y, type.measure="deviance", alpha=1,nfolds = 5)
plot(cv_fit)

lambda.min <- cv_fit$lambda.min
lambda.1se <- cv_fit$lambda.1se
lambda.min
lambda.1se

#指定λ值重新构建模型(通过λ值筛选基因)：
model_lasso_min<- glmnet(x, y, alpha = 1, lambda = lambda.min)
model_lasso_1se<- glmnet(x, y, alpha = 1, lambda = lambda.1se)

#拎出模型使用的基因(存放在beta中)：
head(model_lasso_min$beta)#"."表示这个基因没有被使用

gene_min<- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]#as.numeric后"."会转化为0
gene_1se<- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(gene_min)
length(gene_1se)

lasso.prob <- predict(cv_fit, newx = x,
                      s= c(lambda.min,lambda.1se) )
df<- as.data.frame(cbind(y ,lasso.prob)) #将预测结果和真实生死结果合并，并转换为数据框
colnames(df) <- c("event","pre_min","pre_1se") #列名修改
head(df)

# 计算R²
r2_min <- 1 - (sum((df$event - df$pre_min)^2) / sum((df$event - mean(df$event))^2))
r2_1se <- 1 - (sum((df$event - df$pre_1se)^2) / sum((df$event - mean(df$event))^2))

cat("R² (lambda.min):", r2_min, "\n")
cat("R² (lambda.1se):", r2_1se, "\n")


# Extract selected features
feature_all <- as.data.frame(as.matrix(coef(cv_fit, s=cv_fit$lambda.min)))
colnames(feature_all) <- "coff"
feature_opt <- feature_all %>% filter(abs(coff) > 0)
selected_features <- rownames(feature_opt)

print(selected_features)
print(feature_opt)

#===============================================================================
# 2. Regression Models
#===============================================================================
# Remove intercept from selected features if present
selected_features <- selected_features[selected_features != "(Intercept)"]

# Prepare data for regression
model_data <- complete_data[, c(selected_features, "vision_improvement")]

# Set up 5-fold cross-validation control
ctrl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "final",
  summaryFunction = defaultSummary
)

# # Set up LOOCV control
# ctrl <- trainControl(
#   method = "LOOCV",
#   savePredictions = "final",
#   summaryFunction = defaultSummary
# )

# Train Linear Regression
set.seed(123)
lm_model <- train(
  vision_improvement ~ .,
  data = model_data,
  method = "lm",
  trControl = ctrl,
  metric = "RMSE"
)

# Train Random Forest
set.seed(123)
rf_model <- train(
  vision_improvement ~ .,
  data = model_data,
  method = "rf",
  trControl = ctrl,
  metric = "RMSE",
  importance = TRUE
)

# Train XGBoost
set.seed(123)
xgb_model <- train(
  vision_improvement ~ .,
  data = model_data,
  method = "xgbTree",
  trControl = ctrl,
  metric = "RMSE",
  tuneLength = 5
)


# 尝试弹性网络，但使用较少的参数组合
elastic_net_grid_simple <- expand.grid(
  alpha = c(0.2, 0.5, 0.8),    # 少量alpha值
  lambda = seq(0.01, 0.1, length.out = 5)  # 少量lambda值
)

set.seed(123)
enet_model <- train(
  vision_improvement ~ .,
  data = model_data,
  method = "glmnet",
  trControl = ctrl,
  tuneGrid = elastic_net_grid_simple,
  metric = "RMSE"
)

#===============================================================================
# 3. Evaluate Regression Models
#===============================================================================
# Calculate metrics for Linear Regression
lm_performance <- lm_model$pred %>%
  group_by(Resample) %>%
  summarise(
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    R2 = cor(obs, pred)^2
  )

# Calculate metrics for Random Forest
rf_performance <- rf_model$pred %>%
  group_by(Resample) %>%
  summarise(
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    R2 = cor(obs, pred)^2
  )

# Calculate metrics for XGBoost
xgb_performance <- xgb_model$pred %>%
  group_by(Resample) %>%
  summarise(
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    R2 = cor(obs, pred)^2
  )

# Calculate metrics for enet
enet_performance <- enet_model$pred %>%
  group_by(Resample) %>%
  summarise(
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    R2 = cor(obs, pred)^2
  )


# Calculate overall performance for each model
lm_overall <- summarise(lm_performance,
                        RMSE_mean = mean(RMSE),
                        RMSE_sd = sd(RMSE),
                        MAE_mean = mean(MAE),
                        MAE_sd = sd(MAE),
                        R2_mean = mean(R2),
                        R2_sd = sd(R2))

rf_overall <- summarise(rf_performance,
                        RMSE_mean = mean(RMSE),
                        RMSE_sd = sd(RMSE),
                        MAE_mean = mean(MAE),
                        MAE_sd = sd(MAE),
                        R2_mean = mean(R2),
                        R2_sd = sd(R2))

xgb_overall <- summarise(xgb_performance,
                         RMSE_mean = mean(RMSE),
                         RMSE_sd = sd(RMSE),
                         MAE_mean = mean(MAE),
                         MAE_sd = sd(MAE),
                         R2_mean = mean(R2),
                         R2_sd = sd(R2))

enet_overall <- summarise(enet_performance,
                         RMSE_mean = mean(RMSE),
                         RMSE_sd = sd(RMSE),
                         MAE_mean = mean(MAE),
                         MAE_sd = sd(MAE),
                         R2_mean = mean(R2),
                         R2_sd = sd(R2))

print(lm_overall)
print(rf_overall)
print(xgb_overall)
print(enet_overall)


# # 为每个模型计算整体性能指标
# # 线性回归
# lm_metrics <- data.frame(
#   obs = lm_model$pred$obs,
#   pred = lm_model$pred$pred
# )
# lm_overall <- data.frame(
#   RMSE = sqrt(mean((lm_metrics$obs - lm_metrics$pred)^2)),
#   MAE = mean(abs(lm_metrics$obs - lm_metrics$pred)),
#   R2 = cor(lm_metrics$obs, lm_metrics$pred)^2
# )
# 
# # 随机森林
# rf_metrics <- data.frame(
#   obs = rf_model$pred$obs,
#   pred = rf_model$pred$pred
# )
# rf_overall <- data.frame(
#   RMSE = sqrt(mean((rf_metrics$obs - rf_metrics$pred)^2)),
#   MAE = mean(abs(rf_metrics$obs - rf_metrics$pred)),
#   R2 = cor(rf_metrics$obs, rf_metrics$pred)^2
# )
# 
# # XGBoost
# xgb_metrics <- data.frame(
#   obs = xgb_model$pred$obs,
#   pred = xgb_model$pred$pred
# )
# xgb_overall <- data.frame(
#   RMSE = sqrt(mean((xgb_metrics$obs - xgb_metrics$pred)^2)),
#   MAE = mean(abs(xgb_metrics$obs - xgb_metrics$pred)),
#   R2 = cor(xgb_metrics$obs, xgb_metrics$pred)^2
# )
# 
# # Elastic Net
# enet_metrics <- data.frame(
#   obs = enet_model$pred$obs,
#   pred = enet_model$pred$pred
# )
# enet_overall <- data.frame(
#   RMSE = sqrt(mean((enet_metrics$obs - enet_metrics$pred)^2)),
#   MAE = mean(abs(enet_metrics$obs - enet_metrics$pred)),
#   R2 = cor(enet_metrics$obs, enet_metrics$pred)^2
# )
# 
# # 打印结果
# print("Linear Regression Performance:")
# print(lm_overall)
# print("Random Forest Performance:")
# print(rf_overall)
# print("XGBoost Performance:")
# print(xgb_overall)
# print("Elastic Net Performance:")
# print(enet_overall)


# Separate plots for each model
library(ggplot2)

# Linear Regression Plot
lm_r2 <- paste0("R² = ", round(lm_overall$R2_mean, 3), 
                " (±", round(lm_overall$R2_sd, 3), ")")
lm_plot <- ggplot(lm_model$pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "blue", 
              fill = "lightblue", alpha = 0.2, 
              se = TRUE, level = 0.95) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Observed Vision Improvement",
    y = "Predicted Vision Improvement",
    title = paste0("Linear Regression\n", lm_r2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Random Forest Plot
rf_r2 <- paste0("R² = ", round(rf_overall$R2_mean, 3),
                " (±", round(rf_overall$R2_sd, 3), ")")
rf_plot <- ggplot(rf_model$pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "blue", 
              fill = "lightblue", alpha = 0.2, 
              se = TRUE, level = 0.95) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Observed Vision Improvement",
    y = "Predicted Vision Improvement",
    title = paste0("Random Forest\n", rf_r2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# XGBoost Plot
xgb_r2 <- paste0("R² = ", round(xgb_overall$R2_mean, 3),
                 " (±", round(xgb_overall$R2_sd, 3), ")")
xgb_plot <- ggplot(xgb_model$pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.5, color = "purple") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "blue", 
              fill = "lightblue", alpha = 0.2, 
              se = TRUE, level = 0.95) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Observed Vision Improvement",
    y = "Predicted Vision Improvement",
    title = paste0("XGBoost\n", xgb_r2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Elastic Net Plot
enet_r2 <- paste0("R² = ", round(enet_overall$R2_mean, 3),
                  " (±", round(enet_overall$R2_sd, 3), ")")
enet_plot <- ggplot(enet_model$pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.5, color = "orange") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "blue", 
              fill = "lightblue", alpha = 0.2, 
              se = TRUE, level = 0.95) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Observed Vision Improvement",
    y = "Predicted Vision Improvement",
    title = paste0("Elastic Net\n", enet_r2)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Display each plot
print(lm_plot)
print(rf_plot)
print(xgb_plot)
print(enet_plot)





#===============================================================================
# 4. Classification Models
#===============================================================================
# Prepare data for classification
complete_data_class <- complete_data[, c(selected_features, "vision_improved")]
complete_data_class$vision_improved <- factor(complete_data_class$vision_improved,
                                           levels = c("0", "1"),
                                           labels = c("No", "Yes"))

# Set up cross-validation for classification
ctrl_class <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final",
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Train Logistic Regression
set.seed(123)
logit_model <- train(
  vision_improved ~ .,
  data = complete_data_class,
  method = "glm",
  family = "binomial",
  trControl = ctrl_class,
  metric = "ROC"
)

# Train Random Forest Classification
set.seed(123)
rf_class_model <- train(
  vision_improved ~ .,
  data = complete_data_class,
  method = "rf",
  trControl = ctrl_class,
  metric = "ROC",
  importance = TRUE
)

# Train XGBoost Classification
set.seed(123)
xgb_class_model <- train(
  vision_improved ~ .,
  data = complete_data_class,
  method = "xgbTree",
  trControl = ctrl_class,
  metric = "ROC",
  tuneLength = 5
)

# Train Elastic Net Classification
set.seed(123)
enet_class_model <- train(
  vision_improved ~ .,
  data = complete_data_class,
  method = "glmnet",
  trControl = ctrl_class,
  metric = "ROC",
  tuneLength = 10,
  # Add glmnet specific parameters
  family = "binomial"
)

#===============================================================================
# 5. Evaluate Classification Models
#===============================================================================
# Calculate metrics for Logistic Regression
logit_performance <- logit_model$pred %>%
  group_by(Resample) %>%
  summarise(
    Accuracy = mean(pred == obs),
    Sensitivity = sum(pred == "Yes" & obs == "Yes") / sum(obs == "Yes"),
    Specificity = sum(pred == "No" & obs == "No") / sum(obs == "No"),
    AUC = as.numeric(roc(obs, Yes)$auc)
  )

# Calculate metrics for Random Forest
rf_class_performance <- rf_class_model$pred %>%
  group_by(Resample) %>%
  summarise(
    Accuracy = mean(pred == obs),
    Sensitivity = sum(pred == "Yes" & obs == "Yes") / sum(obs == "Yes"),
    Specificity = sum(pred == "No" & obs == "No") / sum(obs == "No"),
    AUC = as.numeric(roc(obs, Yes)$auc)
  )

# Calculate metrics for XGBoost
xgb_class_performance <- xgb_class_model$pred %>%
  group_by(Resample) %>%
  summarise(
    Accuracy = mean(pred == obs),
    Sensitivity = sum(pred == "Yes" & obs == "Yes") / sum(obs == "Yes"),
    Specificity = sum(pred == "No" & obs == "No") / sum(obs == "No"),
    AUC = as.numeric(roc(obs, Yes)$auc)
  )

# Calculate metrics for Elastic Net
enet_class_performance <- enet_class_model$pred %>%
  group_by(Resample) %>%
  summarise(
    Accuracy = mean(pred == obs),
    Sensitivity = sum(pred == "Yes" & obs == "Yes") / sum(obs == "Yes"),
    Specificity = sum(pred == "No" & obs == "No") / sum(obs == "No"),
    AUC = as.numeric(roc(obs, Yes)$auc)
  )


print("=== Classification Models Performance ===")
print("Logistic Regression:")
print(summary(logit_performance))
print("Random Forest:")
print(summary(rf_class_performance))
print("XGBoost:")
print(summary(xgb_class_performance))
print("Elastic Net:")
print(summary(enet_class_performance))

#===============================================================================
# 6. Visualizations
#===============================================================================


# Classification Models - ROC Curves
# Set up the plot
plot(0, 0, type="n", 
     xlim=c(100,0), 
     ylim=c(0,100),
     xlab="Specificity (%)", 
     ylab="Sensitivity (%)",
     main="ROC Curves for Classification Models",
     axes=FALSE)

# Add axes
axis(1, at=seq(0,100,20))
axis(2, at=seq(0,100,20))

# Plot ROC curve for Logistic Regression
plot.roc(
  logit_model$pred$obs,
  logit_model$pred$Yes,
  col = "#A31621",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "Logistic Regression: %.1f%%",
  print.auc.y = 50,
  add = TRUE
)

# Plot ROC curve for Random Forest
plot.roc(
  rf_class_model$pred$obs,
  rf_class_model$pred$Yes,
  col = "#2694ab",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "Random Forest: %.1f%%",
  print.auc.y = 45,
  add = TRUE
)

# Plot ROC curve for XGBoost
plot.roc(
  xgb_class_model$pred$obs,
  xgb_class_model$pred$Yes,
  col = "#4CAF50",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "XGBoost: %.1f%%",
  print.auc.y = 40,
  add = TRUE
)

# Add Elastic Net ROC curve
plot.roc(
  enet_class_model$pred$obs,
  enet_class_model$pred$Yes,
  col = "#FFA500",
  percent = TRUE,
  lwd = 2,
  print.auc = TRUE,
  print.auc.cex = 1,
  print.auc.pattern = "Elastic Net: %.1f%%",
  print.auc.y = 25,
  add = TRUE
)

# Add legend
legend("bottomright", 
       legend = c("Logistic Regression", "Random Forest", "XGBoost", "Elastic Net"),
       col = c("#A31621", "#2694ab", "#4CAF50","#FFA500"),
       lwd = 2)



all_models_summary <- bind_rows(
  logit_performance %>% mutate(Model = "Logistic Regression"),
  rf_class_performance %>% mutate(Model = "Random Forest"),
  xgb_class_performance %>% mutate(Model = "XGBoost"),
  enet_class_performance %>% mutate(Model = "Elastic Net")
) %>%
  group_by(Model) %>%
  summarise(
    Accuracy = mean(Accuracy),
    Accuracy_SD = sd(Accuracy),
    Sensitivity = mean(Sensitivity),
    Sensitivity_SD = sd(Sensitivity),
    Specificity = mean(Specificity),
    Specificity_SD = sd(Specificity),
    AUC = mean(AUC),
    AUC_SD = sd(AUC)
  )

# Print the summary table
print(all_models_summary)
