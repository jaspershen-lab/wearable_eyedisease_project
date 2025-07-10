library(tidyverse)
library(lubridate)
library(ggplot2)
library(tidymass)
library(r4projects)
library(ggside)
library(patchwork)

setwd(get_project_wd())
rm(list = ls())

load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/7_figures/figure1/heart_rate_dot_plot_ppv_diabetes", recursive = TRUE)
setwd("3_data_analysis/7_figures/figure1//heart_rate_dot_plot_ppv_diabetes")


# Filter for PPV diabetes patients first - ensure consistency
ppv_diabetes_baseline <- baseline_info %>%
  mutate(
    surgery_type = case_when(
      surgery_1..0.PI.1.other. == 0 ~ "Anterior (Cataract)",
      surgery_1..0.PI.1.other. == 1 ~ "Posterior (PPV)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(surgery_type == "Posterior (PPV)" & diabetes_history == 1)

# Get heart rate data subject IDs
heart_rate_ids <- unique(heart_rate_data@sample_info$subject_id)

# Find intersection - PPV diabetes patients who also have heart rate data
ppv_diabetes_ids <- ppv_diabetes_baseline %>%
  filter(toupper(ID) %in% heart_rate_ids) %>%
  pull(ID) %>%
  toupper()

cat("Total PPV diabetes patients in baseline:", nrow(ppv_diabetes_baseline), "\n")
cat("PPV diabetes patients with heart rate data:", length(ppv_diabetes_ids), "\n")
print("PPV diabetes patient IDs with heart rate data:")
print(ppv_diabetes_ids)

# Function to calculate daily hour coverage for PPV diabetes patients only
calculate_hourly_coverage_ppv_diabetes <- function(heart_rate_data, ppv_diabetes_ids) {
  # Convert heart_rate_data from mass_dataset to data frame for processing
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame() %>%
    # Filter for PPV diabetes patients only
    filter(subject_id %in% ppv_diabetes_ids)
  
  cat("Heart rate data subjects after filtering:", length(unique(hr_df$subject_id)), "\n")
  
  # Create ID mapping for anonymization (1-20)
  unique_ids <- sort(unique(hr_df$subject_id))
  id_mapping <- data.frame(
    original_id = unique_ids,
    anonymous_id = paste0("P", sprintf("%02d", 1:length(unique_ids)))
  )
  
  cat("ID mapping created for", nrow(id_mapping), "patients\n")
  print(id_mapping)
  
  # Get surgery dates for PPV diabetes patients
  baseline_info_filtered <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1),
      ID = toupper(ID)
    ) %>%
    filter(ID %in% ppv_diabetes_ids) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Calculate hours coverage
  hourly_coverage <- hr_df %>%
    # Join with surgery dates
    left_join(baseline_info_filtered, by = c("subject_id" = "ID")) %>%
    # Add anonymous ID mapping
    left_join(id_mapping, by = c("subject_id" = "original_id")) %>%
    # Calculate days relative to surgery
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = hour(measure_time)
    ) %>%
    # Filter to desired range (-4 to 30 days)
    filter(
      day_point >= -4,
      day_point <= 30
    ) %>%
    # Count unique hours per subject per day using anonymous ID
    group_by(anonymous_id, day_point) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # Create complete grid with all anonymous subject-day combinations
  all_anonymous_ids <- id_mapping$anonymous_id
  all_days <- seq(-4, 30)
  complete_grid <- expand.grid(
    anonymous_id = all_anonymous_ids,
    day_point = all_days
  )
  
  # Join with actual coverage and fill missing with 0
  final_coverage <- complete_grid %>%
    left_join(hourly_coverage, by = c("anonymous_id", "day_point")) %>%
    mutate(hours_covered = coalesce(hours_covered, 0))
  
  # Return both coverage data and ID mapping
  return(list(
    coverage = final_coverage,
    id_mapping = id_mapping
  ))
}

# Create heatmap function for PPV diabetes patients
create_hr_coverage_heatmap_ppv <- function(hourly_coverage) {
  # Sort anonymous IDs numerically
  hourly_coverage$anonymous_id <- factor(
    hourly_coverage$anonymous_id,
    levels = rev(sort(unique(hourly_coverage$anonymous_id)))
  )
  
  # Create the plot
  p <- ggplot(hourly_coverage, aes(x = day_point, y = anonymous_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,   # 实心圆点，没有边框
               alpha = 0.7) +
    # 分别设置color和fill的渐变
    scale_color_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_fill_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    scale_alpha(
      range = c(0.3, 1),
      guide = "none"
    ) +
    # 设置主题
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = "Heart Rate Data Coverage Pattern - PPV Diabetes Patients",
      subtitle = paste("N =", length(unique(hourly_coverage$anonymous_id)), "diabetic PPV patients"),
      x = "Days Relative to Surgery",
      y = "Patient ID"
    )
  
  return(p)
}

# Run analysis for PPV diabetes patients only
coverage_results <- calculate_hourly_coverage_ppv_diabetes(heart_rate_data, ppv_diabetes_ids)
hourly_coverage_ppv <- coverage_results$coverage
id_mapping <- coverage_results$id_mapping

# Save ID mapping for reference
write.csv(id_mapping, "id_mapping_ppv_diabetes.csv", row.names = FALSE)

p_main <- create_hr_coverage_heatmap_ppv(hourly_coverage_ppv)
print(p_main)

ggsave(filename = "heart_rate_data_dot_plot_ppv_diabetes.pdf", plot = p_main, width = 8, height = 7)

# Create perfect alignment plot for PPV diabetes patients
# Create perfect alignment plot for PPV diabetes patients
create_perfect_alignment_plot_ppv <- function(hourly_coverage) {
  # Sort anonymous IDs numerically 
  hourly_coverage$anonymous_id <- factor(
    hourly_coverage$anonymous_id,
    levels = rev(sort(unique(hourly_coverage$anonymous_id)))
  )
  
  # Calculate summary data for the additional plots
  # Daily totals for each day (for top bar chart)
  daily_totals <- hourly_coverage %>%
    group_by(day_point) %>%
    summarise(daily_hours = sum(hours_covered), .groups = "drop")
  
  # Total hours per subject across all days (for right bar chart)
  subject_totals <- hourly_coverage %>%
    group_by(anonymous_id) %>%
    summarise(total_hours = sum(hours_covered), .groups = "drop")
  
  # Get the range of days and unique subject IDs for consistent axis limits
  day_range <- range(hourly_coverage$day_point)
  unique_subjects <- levels(hourly_coverage$anonymous_id)
  
  # 设置一致的横坐标断点
  x_breaks <- seq(from = ceiling(day_range[1]/5)*5, 
                  to = floor(day_range[2]/5)*5, 
                  by = 5)
  
  # 如果数据范围不包含标准断点，则添加边界值
  if(day_range[1] < min(x_breaks)) {
    x_breaks <- c(day_range[1], x_breaks)
  }
  if(day_range[2] > max(x_breaks)) {
    x_breaks <- c(x_breaks, day_range[2])
  }
  x_breaks <- unique(sort(x_breaks))
  
  # Define shared theme elements for consistency
  base_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 1. Create main heatmap
  p_main <- ggplot(hourly_coverage, aes(x = day_point, y = anonymous_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,
               alpha = 0.7) +
    scale_color_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_fill_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
                       breaks = x_breaks,
                       expand = c(0.01, 0.01)) +
    base_theme +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      legend.position = "right"
    ) +
    labs(
      title = "Heart Rate Data Coverage - PPV Diabetes Patients",
      subtitle = paste("N =", length(unique_subjects), "participants"),
      x = "Days Relative to Surgery",
      y = "Patient ID"
    )
  
  # 2. Create top bar chart
  max_daily_hours <- max(daily_totals$daily_hours) * 1.05
  p_top <- ggplot(daily_totals, aes(x = day_point, y = daily_hours)) +
    geom_col(fill = "#9999cc", width = 0.8) +
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
                       breaks = x_breaks,
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, max_daily_hours),
                       expand = c(0, 0)) +
    base_theme +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)
    ) +
    labs(y = "Daily\nHours")
  
  # 3. Create right bar chart
  max_total_hours <- max(subject_totals$total_hours) * 1.05
  p_right <- ggplot(subject_totals, aes(y = anonymous_id, x = total_hours)) +
    geom_col(fill = "#9999cc") +
    scale_y_discrete(limits = unique_subjects,
                     expand = c(0.01, 0.01)) +
    scale_x_continuous(limits = c(0, max_total_hours),
                       expand = c(0, 0)) +
    base_theme +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)
    ) +
    labs(x = "Total Hours")
  
  # 4. Create empty plot for top-right corner
  p_empty <- ggplot() + 
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Combine plots
  combined_plot <- (p_top | p_empty) / (p_main | p_right) +
    plot_layout(
      widths = c(4, 1),
      heights = c(1, 5),
      guides = "collect"
    ) &
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
  
  return(combined_plot)
}

# Create the combined plot
combined_plot_ppv <- create_perfect_alignment_plot_ppv(hourly_coverage_ppv)
print(combined_plot_ppv)

ggsave(filename = "heart_rate_data_perfectly_aligned_ppv_diabetes.pdf", 
       plot = combined_plot_ppv, 
       width = 14,
       height = 8,  # Reduced height since fewer subjects
       dpi = 300)

# Create histogram of participants with data per day for PPV diabetes patients
participants_per_day_ppv <- hourly_coverage_ppv %>%
  filter(hours_covered > 0) %>%
  group_by(day_point) %>%
  summarise(
    participant_count = n_distinct(anonymous_id)
  )

p_histogram_ppv <- ggplot(participants_per_day_ppv, aes(x = day_point, y = participant_count)) +
  geom_bar(stat = "identity", fill = "#9999cc", color = "black", alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Number of PPV Diabetes Participants with Wearable Data per Day",
    subtitle = paste("N =", length(unique(hourly_coverage_ppv$anonymous_id)), "diabetic PPV patients"),
    x = "Days Relative to Surgery",
    y = "Number of Participants"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_continuous(breaks = seq(min(participants_per_day_ppv$day_point), 
                                  max(participants_per_day_ppv$day_point), 
                                  by = 5)) +
  scale_y_continuous(limits = c(0, length(unique(hourly_coverage_ppv$anonymous_id))))

print(p_histogram_ppv)
ggsave(filename = "heart_rate_data_histogram_ppv_diabetes.pdf", plot = p_histogram_ppv, width = 8, height = 7)

# Create combined plot
combined_plot_main <- p_main + p_histogram_ppv + 
  plot_layout(
    design = "
    AA
    BB
    ",
    heights = c(3, 1)
  )

print(combined_plot_main)
ggsave(filename = "heart_rate_data_combined_ppv_diabetes.pdf", plot = combined_plot_main, width = 8, height = 7)

# Calculate presurgery days for PPV diabetes patients
calculate_presurgery_days_ppv <- function(hourly_coverage) {
  presurgery_data <- hourly_coverage %>%
    filter(day_point < 0) %>%
    filter(hours_covered > 0) %>%
    group_by(anonymous_id) %>%
    summarise(
      presurgery_days_worn = n_distinct(day_point),
      .groups = "drop"
    ) %>%
    arrange(anonymous_id)
  
  return(presurgery_data)
}

presurgery_days_ppv <- calculate_presurgery_days_ppv(hourly_coverage_ppv)

print(presurgery_days_ppv)
write.csv(presurgery_days_ppv, "presurgery_wearable_days_ppv_diabetes.csv", row.names = FALSE)

# Create presurgery days bar plot
p_presurgery_ppv <- ggplot(presurgery_days_ppv, aes(x = anonymous_id, y = presurgery_days_worn)) +
  geom_bar(stat = "identity", fill = "#9999cc", color = "black", alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Presurgery Wearable Days - PPV Diabetes Patients",
    subtitle = paste("N =", nrow(presurgery_days_ppv), "diabetic PPV patients"),
    x = "Patient ID",
    y = "Presurgery Days"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

print(p_presurgery_ppv)
ggsave(filename = "presurgery_wearable_days_ppv_diabetes.pdf", plot = p_presurgery_ppv, width = 8, height = 6)

# Calculate summary statistics for PPV diabetes patients
presurgery_summary_ppv <- presurgery_days_ppv %>%
  summarise(
    mean_days = mean(presurgery_days_worn),
    median_days = median(presurgery_days_worn),
    min_days = min(presurgery_days_worn),
    max_days = max(presurgery_days_worn),
    sd_days = sd(presurgery_days_worn)
  )

print("Presurgery days summary for PPV diabetes patients:")
print(presurgery_summary_ppv)
write.csv(presurgery_summary_ppv, "presurgery_days_summary_ppv_diabetes.csv", row.names = FALSE)

cat("\n============ ANALYSIS COMPLETE ============\n")
cat("Generated files for PPV diabetes patients:\n")
cat("1. heart_rate_data_dot_plot_ppv_diabetes.pdf\n")
cat("2. heart_rate_data_perfectly_aligned_ppv_diabetes.pdf\n")
cat("3. heart_rate_data_histogram_ppv_diabetes.pdf\n")
cat("4. heart_rate_data_combined_ppv_diabetes.pdf\n")
cat("5. presurgery_wearable_days_ppv_diabetes.pdf\n")
cat("6. presurgery_wearable_days_ppv_diabetes.csv\n")
cat("7. presurgery_days_summary_ppv_diabetes.csv\n")
cat("8. id_mapping_ppv_diabetes.csv - ID mapping table\n").x = element_blank(),
axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
).x = element_blank(),
axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
)

print(p_presurgery_ppv)
ggsave(filename = "presurgery_wearable_days_ppv_diabetes.pdf", plot = p_presurgery_ppv, width = 8, height = 6)

# Calculate summary statistics for PPV diabetes patients
presurgery_summary_ppv <- presurgery_days_ppv %>%
  summarise(
    mean_days = mean(presurgery_days_worn),
    median_days = median(presurgery_days_worn),
    min_days = min(presurgery_days_worn),
    max_days = max(presurgery_days_worn),
    sd_days = sd(presurgery_days_worn)
  )

print("Presurgery days summary for PPV diabetes patients:")
print(presurgery_summary_ppv)
write.csv(presurgery_summary_ppv, "presurgery_days_summary_ppv_diabetes.csv", row.names = FALSE)