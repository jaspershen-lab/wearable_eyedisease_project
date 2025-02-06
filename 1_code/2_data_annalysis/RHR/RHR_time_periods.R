library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


########read data
load(
  "3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/RHR/", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/RHR/")


###############rhr by time window
# Process baseline info
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# Function to calculate RHR for a specific label
calculate_rhr_by_label <- function(data, label_value) {
  # Get sample info with the specific label
  filtered_sample_info <- data@sample_info %>%
    dplyr::filter(label == label_value)
  
  # Get corresponding heart rate data
  heart_rate_values <- data@expression_data[1, filtered_sample_info$sample_id] %>%
    as.numeric()
  
  # Create data frame with sample info and heart rate values
  data.frame(
    sample_id = filtered_sample_info$sample_id,
    heart_rate = heart_rate_values,
    timestamp = filtered_sample_info$measure_time,
    subject_id = filtered_sample_info$subject_id
  ) %>%
    arrange(subject_id, timestamp)
}

# Calculate RHR for both labels
rhr_label_1 <- calculate_rhr_by_label(heart_rate_data, "<1")
rhr_label_50 <- calculate_rhr_by_label(heart_rate_data, "<50")

# Function to process RHR data with time periods
process_rhr_data <- function(rhr_data, label_suffix) {
  rhr_data %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      days_to_surgery = as.numeric(difftime(surgery_time_1, timestamp, units = "days"))
    ) %>%
    group_by(subject_id) %>%
    mutate(
      has_future_surgery = all(days_to_surgery > 0),
      min_days = min(days_to_surgery),
      max_days = max(days_to_surgery)
    ) %>%
    mutate(
      time_period = case_when(
        is.na(surgery_time_1) ~ "no_surgery_date",
        days_to_surgery >= 0 & days_to_surgery < 3 ~ "pre_surgery_3d",
        days_to_surgery >= 3 & days_to_surgery < 7 ~ "pre_surgery_3d_to_7d",
        days_to_surgery >= 7 ~ "pre_surgery_over_7d",
        days_to_surgery < 0 & days_to_surgery >= -7 ~ "post_surgery_1w",
        days_to_surgery < -7 & days_to_surgery >= -30 ~ "post_surgery_1m",
        days_to_surgery < -30 ~ "post_surgery_over_1m",
        TRUE ~ "other"
      ),
      time_period = case_when(
        has_future_surgery & max_days >= 7 ~ case_when(
          days_to_surgery <= min_days + 3 ~ "pre_surgery_3d",
          days_to_surgery <= min_days + 7 ~ "pre_surgery_3d_to_7d",
          TRUE ~ "pre_surgery_over_7d"
        ),
        TRUE ~ time_period
      )
    ) %>%
    ungroup() %>%
    group_by(subject_id, time_period) %>%
    summarise(
      n_measurements = n(),
      mean_hr = mean(heart_rate),
      sd_hr = sd(heart_rate),
      .groups = "drop"
    ) %>%
    mutate(label = label_suffix)
}

# Process both datasets
rhr_summary_1 <- process_rhr_data(rhr_label_1, "steps_1")
rhr_summary_50 <- process_rhr_data(rhr_label_50, "steps_50")

# Combine the results
combined_rhr_summary <- bind_rows(rhr_summary_1, rhr_summary_50) %>%
  mutate(
    label_time = paste(label, time_period, sep = "_")
  ) %>%
  dplyr::select(subject_id, label_time, mean_hr, sd_hr, n_measurements)


# Create wide format table with surgery date
rhr_summary_wide <- combined_rhr_summary %>%
  dplyr::select(subject_id, label_time, mean_hr) %>%
  pivot_wider(
    names_from = label_time,
    values_from = mean_hr,
    names_prefix = "rhr_"
  ) %>%
  # Add surgery date
  left_join(
    baseline_info_processed %>%
      mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
      dplyr::select(ID, surgery_date),
    by = c("subject_id" = "ID")
  ) %>%
  # Reorder columns to put surgery_date at the beginning
  dplyr::select(subject_id, surgery_date, everything())

# Add summary statistics
rhr_summary_stats <- combined_rhr_summary %>%
  group_by(label_time) %>%
  summarise(
    mean_rhr = mean(mean_hr, na.rm = TRUE),
    sd_rhr = sd(mean_hr, na.rm = TRUE),
    median_rhr = median(mean_hr, na.rm = TRUE),
    n_subjects = n(),
    .groups = "drop"
  )
