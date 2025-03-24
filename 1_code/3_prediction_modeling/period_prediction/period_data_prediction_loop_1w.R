# Comprehensive Analysis for Vision Improvement Prediction
# Modified to match the approach in individual model script

library(glmnet)
library(caret)
library(randomForest)
library(xgboost)
library(tidyverse)
library(pROC)
library(ggplot2)
library(mice)
library(VIM)
setwd(get_project_wd())
rm(list = ls())

# Define the time windows we want to analyze
time_windows <- c(
  "pre_3d",
  "pre_3d_7d", 
  "pre_7d_all", 
  "pre_all", 
  "post_4to6d",
  "post_6d"
)

# Load time window datasets
window_datasets <- list()

for (window in time_windows) {
  file_path <- file.path("3_data_analysis/3_prediction_modeling/1w_prediction/time_period_data", 
                         paste0(window, "_data.csv"))
  
  tryCatch({
    window_datasets[[window]] <- read_csv(file_path, show_col_types = FALSE)
    cat("Loaded dataset for", window, "with", nrow(window_datasets[[window]]), "rows\n")
  }, error = function(e) {
    cat("Could not load dataset for", window, ":", e$message, "\n")
  })
}

# Function to select features based on window prefix
get_feature_prefix <- function(window_label) {
  window_prefix_map <- list(
    "pre_3d" = "pre_surgery_3d",
    "pre_3d_7d" = "pre_surgery_3d_to_7d", 
    "pre_7d_all" = "pre_surgery_7d_all", 
    "pre_all" = "pre_surgery_all", 
    "post_4to6d" = "post_surgery_4to6d",
    "post_6d"= "post_surgery_6d"
  )
  
  prefix <- window_prefix_map[[window_label]]
  if(is.null(prefix)) {
    prefix <- window_label
  }
  
  return(prefix)
}

# Function to analyze data from a specific time window
analyze_time_window_data <- function(window_data, window_label) {
  cat("\n=== Starting analysis for", window_label, "===\n")
  
  # Remove rows where vision_improvement is NA (our target variable)
  window_data <- window_data[!is.na(window_data$vision_improvement), ]
  cat("After removing NA vision_improvement:", nrow(window_data), "rows remain\n")
  
  # Get prefix for feature names
  prefix <- get_feature_prefix(window_label)
  cat("Using prefix for window", window_label, ":", prefix, "\n")
  
  # Convert season to factor if needed
  if ("season" %in% names(window_data) && !is.factor(window_data$season)) {
    window_data$season <- as.factor(window_data$season)
  }
  
  # Select specific features based on time window, similar to individual script
  features_for_imputation <- window_data %>%
    dplyr::select(
      # RHR features
      matches(paste0("^", prefix, "_mean_rhr_steps_1$")),
      matches(paste0("^", prefix, "_min_rhr_steps_1$")),
      matches(paste0("^", prefix, "_max_rhr_steps_1$")),
      matches(paste0("^", prefix, "_median_rhr_steps_1$")),
      matches(paste0("^", prefix, "_sd_rhr_steps_1$")),
      matches(paste0("^", prefix, "_iqr_rhr_steps_1$")),
      # BO features
      matches(paste0("^", prefix, "_bo_mean$")),
      matches(paste0("^", prefix, "_bo_min$")),
      matches(paste0("^", prefix, "_bo_max$")),
      matches(paste0("^", prefix, "_bo_median$")),
      matches(paste0("^", prefix, "_bo_sd$")),
      matches(paste0("^", prefix, "_bo_iqr$")),
      # Sleep features
      matches(paste0("^", prefix, "_deep_sleep_total$")), 
      matches(paste0("^", prefix, "_total_sleep_total$")),
      # Steps features
      matches(paste0("^", prefix, "_steps_total$")),
      matches(paste0("^", prefix, "_steps_mean$")),
      matches(paste0("^", prefix, "_steps_max$")),
      # Demographics and medical history
      age, gender, cataract_2, dm_2, hypertension_2, pre_vision, season,
      # Target variable
      vision_improvement
    )
  
  # Check for missing values
  missing_summary <- sapply(features_for_imputation, function(x) sum(is.na(x)))
  cat("Missing values per feature:", paste(names(missing_summary)[missing_summary > 0], 
                                           "=", missing_summary[missing_summary > 0], 
                                           collapse=", "), "\n")
  
  # Only perform imputation if missing values exist
  if(sum(missing_summary) > 0) {
    # Set random seed for reproducibility
    set.seed(1234)
    
    # Configure imputation methods
    imputation_methods <- make.method(features_for_imputation)
    
    # Create 5 imputed datasets
    imp <- mice(features_for_imputation, m=5, method=imputation_methods, 
                maxit=50, seed=1234, printFlag=FALSE)
    
    # Create a complete dataset using the first imputation
    imputed_data <- complete(imp, 1)
    
    # Replace missing values in the original dataset
    window_data_imputed <- window_data
    for(col in names(features_for_imputation)) {
      if(sum(is.na(window_data[[col]])) > 0) {
        window_data_imputed[[col]] <- imputed_data[[col]]
      }
    }
    cat("Imputation completed\n")
  } else {
    window_data_imputed <- window_data
    cat("No imputation needed\n")
  }
  
  # Extract feature names for model without "vision_improvement"
  feature_names <- setdiff(names(features_for_imputation), "vision_improvement")
  cat("Selected features:", paste(feature_names[1:min(5, length(feature_names))], collapse=", "), "...\n")
  
  # LASSO feature selection
  # Create model matrix for LASSO
  x <- as.matrix(window_data_imputed[, feature_names])
  y <- as.numeric(window_data_imputed$vision_improvement)
  
  # Run LASSO with cross-validation
  set.seed(1234)
  cv_fit <- cv.glmnet(x, y, type.measure="deviance", alpha=1, nfolds=5)
  
  # Extract selected features from LASSO
  feature_all <- as.data.frame(as.matrix(coef(cv_fit, s=cv_fit$lambda.min)))
  colnames(feature_all) <- "coff"
  feature_opt <- feature_all %>% filter(abs(coff) > 0)
  lasso_selected_features <- rownames(feature_opt)
  
  # Remove intercept from selected features if present
  lasso_selected_features <- lasso_selected_features[lasso_selected_features != "(Intercept)"]
  
  cat("LASSO selected", length(lasso_selected_features), "features\n")
  if(length(lasso_selected_features) > 0) {
    cat("Selected features:", paste(lasso_selected_features[1:min(5, length(lasso_selected_features))], collapse=", "), "...\n")
  }
  
  # Prepare data for regression models
  model_data <- window_data_imputed[, c(lasso_selected_features, "vision_improvement")]
  
  # Set up 5-fold cross-validation control
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    savePredictions = "final",
    summaryFunction = defaultSummary
  )
  
  # Train models with consistent random seeds
  
  # Linear Regression
  set.seed(123)
  lm_model <- train(
    vision_improvement ~ .,
    data = model_data,
    method = "lm",
    trControl = ctrl,
    metric = "RMSE"
  )
  
  # Random Forest
  set.seed(123)
  rf_model <- train(
    vision_improvement ~ .,
    data = model_data,
    method = "rf",
    trControl = ctrl,
    metric = "RMSE"
  )
  
  # XGBoost
  set.seed(123)
  xgb_model <- train(
    vision_improvement ~ .,
    data = model_data,
    method = "xgbTree",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5
  )
  
  # LASSO
  lasso_grid <- expand.grid(
    alpha = 1,
    lambda = seq(0.001, 0.1, length.out = 10)
  )
  
  set.seed(123)
  lasso_model <- train(
    vision_improvement ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = lasso_grid,
    metric = "RMSE"
  )
  
  # Elastic Net
  elastic_net_grid <- expand.grid(
    alpha = c(0.2, 0.5, 0.8),
    lambda = seq(0.01, 0.1, length.out = 5)
  )
  
  set.seed(123)
  enet_model <- train(
    vision_improvement ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = elastic_net_grid,
    metric = "RMSE"
  )
  
  # Calculate performance metrics
  calculate_performance <- function(model) {
    model$pred %>%
      group_by(Resample) %>%
      summarise(
        RMSE = sqrt(mean((obs - pred)^2)),
        MAE = mean(abs(obs - pred)),
        R2 = cor(obs, pred)^2
      ) %>%
      summarise(
        RMSE_mean = mean(RMSE),
        RMSE_sd = sd(RMSE),
        MAE_mean = mean(MAE),
        MAE_sd = sd(MAE),
        R2_mean = mean(R2),
        R2_sd = sd(R2)
      )
  }
  
  # Calculate performance for each model
  lm_performance <- calculate_performance(lm_model)
  rf_performance <- calculate_performance(rf_model)
  xgb_performance <- calculate_performance(xgb_model)
  lasso_performance <- calculate_performance(lasso_model)
  enet_performance <- calculate_performance(enet_model)
  
  # Print model performance
  cat("\nModel performances for", window_label, ":\n")
  cat("Linear Regression R²:", lm_performance$R2_mean, "\n")
  cat("Random Forest R²:", rf_performance$R2_mean, "\n")
  cat("XGBoost R²:", xgb_performance$R2_mean, "\n")
  cat("LASSO R²:", lasso_performance$R2_mean, "\n")
  cat("Elastic Net R²:", enet_performance$R2_mean, "\n")
  
  # Combine results
  results <- bind_rows(
    lm_performance %>% mutate(Model = "Linear Regression"),
    rf_performance %>% mutate(Model = "Random Forest"),
    xgb_performance %>% mutate(Model = "XGBoost"),
    lasso_performance %>% mutate(Model = "LASSO"),
    enet_performance %>% mutate(Model = "Elastic Net")
  ) %>%
    mutate(Window = window_label)
  
  return(results)
}

# Main analysis
all_results <- data.frame()

for (window in time_windows) {
  cat("\n\nProcessing", window, "...\n")
  
  # Check if data is available for this window
  if (!window %in% names(window_datasets)) {
    cat("Skipping", window, "because no data is available\n")
    next
  }
  
  window_data <- window_datasets[[window]]
  
  # Handle categorical variables
  categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2", "vision_improved")
  
  # Convert categorical variables to factors
  for (col in intersect(categorical_vars, names(window_data))) {
    window_data[[col]] <- as.factor(window_data[[col]])
  }
  
  # Convert season to factor
  if ("season" %in% names(window_data) && !is.factor(window_data$season)) {
    window_data$season <- as.factor(window_data$season)
  }
  
  # Check if dataset has enough rows
  if (nrow(window_data) < 5) {
    cat("Skipping", window, "because sample size is too small (n =", nrow(window_data), ")\n")
    next
  }
  
  # Count valid samples (non-NA vision_improvement)
  valid_samples <- sum(!is.na(window_data$vision_improvement))
  if (valid_samples < 5) {
    cat("Skipping", window, "because valid sample size is too small (valid samples =", 
        valid_samples, ")\n")
    next
  }
  
  # Analyze data
  tryCatch({
    window_results <- analyze_time_window_data(window_data, window)
    
    # Add to overall results
    if (!is.null(window_results) && nrow(window_results) > 0) {
      all_results <- bind_rows(all_results, window_results)
      cat(window, "processing complete\n")
    } else {
      cat(window, "did not produce valid results\n")
    }
  }, error = function(e) {
    cat("Error analyzing", window, ":", e$message, "\n")
  })
}

# Display results
cat("\n=== FINAL RESULTS ===\n")
if(nrow(all_results) > 0) {
  # Display R² values by Window and Model
  r2_table <- all_results %>%
    dplyr::select(Window, Model, R2_mean) %>%
    pivot_wider(names_from = Model, values_from = R2_mean)
  
  cat("R² by time window and model:\n")
  print(r2_table)
  
  # Find best performing window for each model
  best_by_model <- all_results %>%
    group_by(Model) %>%
    slice_max(order_by = R2_mean, n = 1) %>%
    dplyr::select(Model, Window, R2_mean)
  
  cat("\nBest time window for each model:\n")
  print(best_by_model)
  
  # Find best model for each window
  best_by_window <- all_results %>%
    group_by(Window) %>%
    slice_max(order_by = R2_mean, n = 1) %>%
    dplyr::select(Window, Model, R2_mean)
  
  cat("\nBest model for each time window:\n")
  print(best_by_window)
}

# Save results
# Set output directory
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/time_period_data/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define window order
window_order <- c(
  "pre_3d",
  "pre_3d_7d", 
  "pre_7d_all", 
  "pre_all", 
  "post_4to6d",
  "post_6d"
)

# Check if all windows are actually present in our data
present_windows <- unique(all_results$Window)
valid_window_order <- intersect(window_order, present_windows)

if(length(valid_window_order) > 1) {
  # Ensure Window is an ordered factor with only the windows we actually have
  all_results$Window <- factor(all_results$Window, 
                               levels = valid_window_order, 
                               ordered = TRUE)
} else {
  # If we only have one window, just use it as is
  cat("Only one window present in results. No reordering needed.\n")
}

# Save all results
write.csv(all_results, file.path(output_dir, "all_model_results.csv"), row.names = FALSE)

# Create and save best window results
best_windows <- all_results %>%
  group_by(Model) %>%
  slice_max(order_by = R2_mean, n = 1) %>%
  dplyr::select(Model, Window, R2_mean, R2_sd, RMSE_mean) %>%
  arrange(desc(R2_mean))

write.csv(best_windows, file.path(output_dir, "best_windows_by_model.csv"), row.names = FALSE)

# Save overall model performance
best_models <- all_results %>%
  group_by(Model) %>%
  summarise(
    Avg_R2 = mean(R2_mean),
    Max_R2 = max(R2_mean),
    Min_R2 = min(R2_mean),
    Std_R2 = sd(R2_mean)
  ) %>%
  arrange(desc(Avg_R2))

write.csv(best_models, file.path(output_dir, "overall_model_performance.csv"), row.names = FALSE)

# Check data before plotting
cat("\n=== DATA FOR VISUALIZATION ===\n")
cat("Checking data for line plot...\n")
print(table(all_results$Window, all_results$Model))

# Create and save line plot
line_plot <- ggplot(all_results, aes(x = Window, y = R2_mean, color = Model, group = Model)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd), width = 0.2) +
  labs(
    title = "Vision Improvement Prediction Performance by Time Window",
    subtitle = "Comparing different models across pre- and post-operative time periods",
    x = "Time Window Relative to Surgery",
    y = "R² (coefficient of determination)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
line_plot
# Save line plot
ggsave(file.path(output_dir, "prediction_performance_line.pdf"), line_plot, width = 10, height = 7)

# Create heatmap data
heatmap_data <- all_results %>%
  dplyr::select(Window, Model, R2_mean)

# Check if we have enough data for a heatmap (needs multiple windows)
if(length(unique(heatmap_data$Window)) > 1) {
  # Create heatmap data in wide format
  heatmap_data_wide <- heatmap_data %>%
    pivot_wider(names_from = Model, values_from = R2_mean)
  
  # Print heatmap data to verify
  cat("\nHeatmap data (wide format):\n")
  print(heatmap_data_wide)
  
  # Create heatmap plot
  models <- c("Linear Regression", "Random Forest", "XGBoost", "LASSO", "Elastic Net")
  available_models <- intersect(models, unique(heatmap_data$Model))
  
  # Create base plot
  heatmap_plot <- ggplot(heatmap_data, aes(x = Window, y = Model, fill = R2_mean)) +
    geom_tile() +
    geom_text(aes(label = round(R2_mean, 3)), color = "white") +
    scale_fill_gradient(low = "pink", high = "darkred") +
    labs(
      title = "R² Heatmap by Time Window and Model",
      x = "Time Window Relative to Surgery",
      y = "Model",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # Save heatmap
  ggsave(file.path(output_dir, "r2_heatmap.pdf"), heatmap_plot, width = 12, height = 6)
} else {
  cat("Not enough unique windows for a meaningful heatmap.\n")
}
heatmap_plot
# Create and save summary bar plot
summary_plot <- ggplot(best_windows, aes(x = reorder(Model, -R2_mean), y = R2_mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd), width = 0.2) +
  geom_text(aes(label = paste0(round(R2_mean, 3), " (", Window, ")")), vjust = -0.5) +
  labs(
    title = "Best Prediction Performance by Model",
    subtitle = "Showing highest R² value and corresponding time window",
    x = "Model",
    y = "Maximum R²"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save summary plot
ggsave(file.path(output_dir, "best_model_summary.pdf"), summary_plot, width = 8, height = 6)
