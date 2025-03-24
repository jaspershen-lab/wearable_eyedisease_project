# Comprehensive Analysis for Vision Improvement Prediction
# Modified to match the approach in daily model script

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
  "post_1w",
  "post_7d_30d",
  "post_day23_30", 
  "post_day27_30" 
)

# Load time window datasets
window_datasets <- list()

for (window in time_windows) {
  file_path <- file.path("3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data", 
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
    "pre_3d_to_7d" = "pre_surgery_3d_to_7d", 
    "pre_7d_all" = "pre_surgery_7d_all", 
    "pre_all" = "pre_surgery_all", 
    "post_1w" = "post_surgery_1w",
    "post_7d_30d" = "post_surgery_7d_to_30d",
    "post_day23_30" = "post_surgery_day23_to_30",
    "post_day27_30" = "post_surgery_day27_to_30"
  )
  
  prefix <- window_prefix_map[[window_label]]
  if(is.null(prefix)) {
    prefix <- window_label
  }
  
  return(prefix)
}

# Function to analyze data from a specific time window - modified to match the daily script approach
analyze_time_window_data <- function(window_data, window_label) {
  cat("\n=== Starting analysis for", window_label, "===\n")
  
  # Remove rows where vision_improvement_1m is NA (our target variable)
  window_data <- window_data[!is.na(window_data$vision_improvement_1m), ]
  cat("After removing NA vision_improvement_1m:", nrow(window_data), "rows remain\n")
  
  # Get prefix for feature names
  prefix <- get_feature_prefix(window_label)
  cat("Using prefix for window", window_label, ":", prefix, "\n")
  
  # Convert season to factor if needed
  if ("season" %in% names(window_data) && !is.factor(window_data$season)) {
    window_data$season <- as.factor(window_data$season)
  }
  
  # Select features for imputation - matching the daily script feature set but with prefixes
  features_for_imputation <- window_data %>%
    dplyr::select(
      # RHR features
      matches(paste0("^", prefix, "_mean_rhr_steps_1$")),
      matches(paste0("^", prefix, "_min_rhr_steps_1$")),
      matches(paste0("^", prefix, "_max_rhr_steps_1$")),
      matches(paste0("^", prefix, "_median_rhr_steps_1$")),
      matches(paste0("^", prefix, "_sd_rhr_steps_1$")),
      matches(paste0("^", prefix, "_iqr_rhr_steps_1$")),
      matches(paste0("^", prefix, "_skew_rhr_steps_1$")), # Added to match daily script
      matches(paste0("^", prefix, "_kurt_rhr_steps_1$")), # Added to match daily script
      # BO features
      matches(paste0("^", prefix, "_bo_mean$")),
      matches(paste0("^", prefix, "_bo_min$")),
      matches(paste0("^", prefix, "_bo_max$")),
      matches(paste0("^", prefix, "_bo_median$")),
      matches(paste0("^", prefix, "_bo_sd$")),
      matches(paste0("^", prefix, "_bo_iqr$")),
      matches(paste0("^", prefix, "_bo_skew$")), # Added to match daily script
      matches(paste0("^", prefix, "_bo_kurt$")), # Added to match daily script
      # Sleep features
      matches(paste0("^", prefix, "_deep_sleep_total$")), 
      matches(paste0("^", prefix, "_total_sleep_total$")),
      # Steps features
      matches(paste0("^", prefix, "_steps_total$")),
      matches(paste0("^", prefix, "_steps_mean$")),
      matches(paste0("^", prefix, "_steps_max$")),
      # Demographics and medical history - these don't have prefixes
      age, gender, cataract_2, dm_2, hypertension_2, pre_vision, season, post_vision_1w, vision_improvement_1w,
      # Target variable
      vision_improvement_1m
    )
  
  # First, check for missing data - exactly as in daily script
  missing_summary <- sapply(features_for_imputation, function(x) sum(is.na(x)))
  
  # Only perform imputation if there are missing values - matching daily script approach
  if(sum(missing_summary) > 0) {
    # Set the seed for reproducibility
    set.seed(1234)
    
    # Configure the imputation method
    imputation_methods <- make.method(features_for_imputation)
    
    # Create 5 imputed datasets
    imp <- mice(features_for_imputation, m=5, method=imputation_methods, 
                maxit=50, seed=1234, printFlag=FALSE)
    
    # Create a complete dataset using the first imputation
    imputed_data <- complete(imp, 1)
    
    # Replace the missing values in the original dataset
    window_data_imputed <- window_data
    for(col in names(features_for_imputation)) {
      if(sum(is.na(window_data[[col]])) > 0) {
        window_data_imputed[[col]] <- imputed_data[[col]]
      }
    }
  } else {
    window_data_imputed <- window_data
  }
  
  # Define selected features - all available features in the imputed data
  # We're essentially replicating the original daily script's feature selection
  # but using the available features from our window dataset
  feature_names <- names(features_for_imputation)
  feature_names <- feature_names[feature_names != "vision_improvement_1m"]
  
  # Add debugging code here, just like in daily script
  cat("Window:", window_label, "- Data columns:", paste(colnames(window_data_imputed)[1:min(5, ncol(window_data_imputed))], "...", collapse=", "), "\n")
  cat("Available features:", paste(feature_names[1:min(5, length(feature_names))], "...", collapse=", "), "\n")
  cat("Rows after filtering:", nrow(window_data_imputed), "\n")
  
  # Identify categorical variables
  categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2", "season")
  categorical_in_model <- intersect(categorical_vars, feature_names)
  cat("Categorical variables in model:", paste(categorical_in_model, collapse=", "), "\n")
  
  # Modify the model matrix creation code - exactly matching daily script approach
  if(length(categorical_in_model) > 0) {
    # Build formula string with categorical variables
    formula_str <- paste("~ 0 +", paste(feature_names, collapse = " + "))
    x_formula <- as.formula(formula_str)
    
    # Print debugging information
    cat("Creating model matrix with formula:", formula_str, "\n")
    
    # Use tryCatch to handle potential errors in model.matrix
    x <- tryCatch({
      model.matrix(x_formula, data = window_data_imputed)
    }, error = function(e) {
      cat("Error in model.matrix:", e$message, "\n")
      cat("Trying alternative approach without model.matrix...\n")
      # If model.matrix fails, use direct matrix conversion
      as.matrix(window_data_imputed[, feature_names])
    })
  } else {
    # No categorical variables, use simple matrix conversion
    cat("No categorical variables in model, using simple matrix conversion\n")
    x <- as.matrix(window_data_imputed[, feature_names])
  }
  
  y <- as.numeric(window_data_imputed$vision_improvement_1m)
  
  # Fit LASSO model with cross-validation - matching daily script
  set.seed(1234)
  cv_fit <- cv.glmnet(x, y, type.measure="deviance", alpha=1, nfolds=5)
  
  # Extract LASSO-selected features - matching daily script
  feature_all <- as.data.frame(as.matrix(coef(cv_fit, s=cv_fit$lambda.min)))
  colnames(feature_all) <- "coff"
  feature_opt <- feature_all %>% filter(abs(coff) > 0)
  lasso_selected_features <- rownames(feature_opt)
  
  # Remove intercept from selected features if present
  lasso_selected_features <- lasso_selected_features[lasso_selected_features != "(Intercept)"]
  
  # Handle LASSO-selected features - matching daily script approach
  categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2", "season")
  
  # Check if LASSO didn't select any features
  if(length(lasso_selected_features) == 0) {
    cat("LASSO didn't select any features, using all available features\n")
    original_selected_features <- feature_names
  } else {
    # Extract original variable names from LASSO-selected features
    original_selected_features <- c()
    
    for(feature in lasso_selected_features) {
      # Handle dummy variable names like gender1 or genderFemale
      for(cat_var in categorical_vars) {
        if(grepl(paste0("^", cat_var), feature)) {
          if(!cat_var %in% original_selected_features) {
            original_selected_features <- c(original_selected_features, cat_var)
          }
          next
        }
      }
      
      # If not a dummy variable, it's a continuous variable, add directly
      if(!any(sapply(categorical_vars, function(cv) grepl(paste0("^", cv), feature)))) {
        original_selected_features <- c(original_selected_features, feature)
      }
    }
  }
  
  # Ensure all selected features are in the dataset
  available_original_features <- intersect(original_selected_features, colnames(window_data_imputed))
  
  # Prepare data for training
  model_data <- window_data_imputed[, c(available_original_features, "vision_improvement_1m")]
  
  # Set up cross-validation - matching daily script
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    savePredictions = "final",
    summaryFunction = defaultSummary
  )
  
  # Check dataset size, if too small, reduce folds - matching daily script
  if(nrow(model_data) < 10) {
    warning(paste(window_label, "dataset sample size too small (n =", nrow(model_data), "), switching to LOOCV"))
    ctrl <- trainControl(
      method = "LOOCV",
      savePredictions = "final",
      summaryFunction = defaultSummary
    )
  }
  
  # Train models - using same seed and parameters as daily script
  set.seed(123)
  
  # Linear Regression
  lm_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "lm",
    trControl = ctrl,
    metric = "RMSE"
  )
  
  # Random Forest
  rf_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "rf",
    trControl = ctrl,
    metric = "RMSE"
  )
  
  # XGBoost
  xgb_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "xgbTree",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5
  )
  
  # LASSO - matching grid from daily script
  lasso_grid <- expand.grid(
    alpha = 1,
    lambda = seq(0.001, 0.1, length.out = 10)
  )
  
  lasso_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = lasso_grid,
    metric = "RMSE"
  )
  
  # Elastic Net - matching grid from daily script
  elastic_net_grid <- expand.grid(
    alpha = c(0.2, 0.5, 0.8),
    lambda = seq(0.01, 0.1, length.out = 5)
  )
  
  enet_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = elastic_net_grid,
    metric = "RMSE"
  )
  
  # Calculate performance metrics for each model - matching daily script function
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
  
  # Get performance for each model
  lm_performance <- calculate_performance(lm_model)
  rf_performance <- calculate_performance(rf_model)
  xgb_performance <- calculate_performance(xgb_model)
  lasso_performance <- calculate_performance(lasso_model)
  enet_performance <- calculate_performance(enet_model)
  
  # Combine results - matching daily script
  results <- bind_rows(
    lm_performance %>% mutate(Model = "Linear Regression"),
    rf_performance %>% mutate(Model = "Random Forest"),
    xgb_performance %>% mutate(Model = "XGBoost"),
    lasso_performance %>% mutate(Model = "LASSO"),
    enet_performance %>% mutate(Model = "Elastic Net")
  ) %>%
    mutate(Window = window_label)  # Using Window instead of Day
  
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
  
  # Handle categorical variables - matching daily script
  categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2", "vision_improved_1m")
  
  # Convert categorical variables to factors
  for (col in intersect(categorical_vars, names(window_data))) {
    window_data[[col]] <- as.factor(window_data[[col]])
  }
  
  # Handle season separately
  if ("season" %in% names(window_data)) {
    window_data$season <- as.factor(window_data$season)
  }
  
  # Handle character columns - matching daily script
  char_cols <- sapply(window_data, is.character)
  if(any(char_cols)) {
    for(col in names(window_data)[char_cols]) {
      if(col != "subject_id" && col != "vision_improved_factor_1m" && 
         !(col %in% categorical_vars)) {
        # Try to convert to numeric
        window_data[[col]] <- as.numeric(window_data[[col]])
      }
    }
  }
  
  # Check if dataset has enough rows - matching daily script checks
  if (nrow(window_data) < 5) {
    cat("Skipping", window, "because sample size is too small (n =", nrow(window_data), ")\n")
    next
  }
  
  # Check for target variable - matching daily script
  if(!"vision_improvement_1m" %in% names(window_data)) {
    cat("Skipping", window, "because target variable vision_improvement_1m is missing\n")
    next
  }
  
  # Count valid samples (non-NA vision_improvement_1m) - matching daily script
  valid_samples <- sum(!is.na(window_data$vision_improvement_1m))
  if (valid_samples < 5) {
    cat("Skipping", window, "because valid sample size is too small (valid samples =", 
        valid_samples, ")\n")
    next
  }
  
  # Analyze data with error handling - matching daily script
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

##### Save results - matching the approach in daily script
# Set output directory
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define window order - similar to day_order in daily script
window_order <- c(
  "pre_3d",
  "pre_3d_7d", 
  "pre_7d_all", 
  "pre_all", 
  "post_1w",
  "post_7d_30d",
  "post_day23_30", 
  "post_day27_30" 
)

# Ensure Window is an ordered factor - matching daily script approach
all_results$Window <- factor(all_results$Window, levels = window_order, ordered = TRUE)

# Save all results
write.csv(all_results, file.path(output_dir, "all_model_results.csv"), row.names = FALSE)

# Create and save best window results - matching daily script
best_windows <- all_results %>%
  group_by(Model) %>%
  slice_max(order_by = R2_mean, n = 1) %>%
  dplyr::select(Model, Window, R2_mean, R2_sd, RMSE_mean) %>%
  arrange(desc(R2_mean))

write.csv(best_windows, file.path(output_dir, "best_windows_by_model.csv"), row.names = FALSE)

# Save overall model performance - matching daily script
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

# Create and save line plot - matching daily script
line_plot <- ggplot(all_results, aes(x = Window, y = R2_mean, color = Model, group = Model)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd), width = 0.2) +
  labs(
    title = "1-Month Vision Improvement Prediction Performance by Time Window",
    subtitle = "Comparing different models across pre- and post-operative time periods",
    x = "Time Windows Relative to Surgery",
    y = "R² (coefficient of determination)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
line_plot
# Save line plot
ggsave(file.path(output_dir, "prediction_performance_line.pdf"), line_plot, width = 10, height = 7)

# Create and save heatmap - matching daily script
heatmap_data <- all_results %>%
  dplyr::select(Window, Model, R2_mean) %>%
  pivot_wider(names_from = Model, values_from = R2_mean)

heatmap_plot <- ggplot(heatmap_data, aes(x = Window)) +
  geom_tile(aes(y = "Linear Regression", fill = `Linear Regression`)) +
  geom_tile(aes(y = "Random Forest", fill = `Random Forest`)) +
  geom_tile(aes(y = "XGBoost", fill = `XGBoost`)) +
  geom_tile(aes(y = "LASSO", fill = `LASSO`)) +
  geom_tile(aes(y = "Elastic Net", fill = `Elastic Net`)) +
  geom_text(aes(y = "Linear Regression", label = round(`Linear Regression`, 3)), color = "white") +
  geom_text(aes(y = "Random Forest", label = round(`Random Forest`, 3)), color = "white") +
  geom_text(aes(y = "XGBoost", label = round(`XGBoost`, 3)), color = "white") +
  geom_text(aes(y = "LASSO", label = round(`LASSO`, 3)), color = "white") +
  geom_text(aes(y = "Elastic Net", label = round(`Elastic Net`, 3)), color = "white") +
  scale_fill_gradient(low = "pink", high = "darkred") +
  labs(
    title = "1-Month Prediction R² Heatmap by Time Window and Model",
    x = "Time Windows Relative to Surgery",
    y = "Model",
    fill = "R²"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
heatmap_plot
# Save heatmap
ggsave(file.path(output_dir, "r2_heatmap.pdf"), heatmap_plot, width = 12, height = 6)

# Create and save summary bar plot - matching daily script
summary_plot <- ggplot(best_windows, aes(x = reorder(Model, -R2_mean), y = R2_mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd), width = 0.2) +
  geom_text(aes(label = paste0(round(R2_mean, 3), " (", Window, ")")), vjust = -0.5) +
  labs(
    title = "Best 1-Month Prediction Performance by Model",
    subtitle = "Showing highest R² value and corresponding time window",
    x = "Model",
    y = "Maximum R²"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save summary plot
ggsave(file.path(output_dir, "best_model_summary.pdf"), summary_plot, width = 8, height = 6)

# Create time window contribution analysis - matching one-week analysis in daily script
windows_contribution <- all_results %>%
  filter(R2_mean > 0) %>%  # Only select valid windows
  dplyr::select(Window, Model, R2_mean) %>%
  group_by(Window) %>%
  summarise(max_R2 = max(R2_mean)) %>%
  arrange(desc(max_R2))

# Export time window contribution analysis
write.csv(windows_contribution, file.path(output_dir, "time_window_contribution.csv"), row.names = FALSE)
