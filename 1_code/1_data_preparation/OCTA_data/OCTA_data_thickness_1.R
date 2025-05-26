library(tidyverse)
library(r4projects)
library(lubridate) # 添加处理日期的包
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Function to extract measurement info from filename
extract_measure_info <- function(filename) {
  # Extract layer name using all possible layers we found
  layer <- str_extract(
    filename,
    paste0(
      "(",
      paste(
        "Choroid",
        "GCL\\+IPL",
        "ILM\\+BM",
        "INL",
        "OuterRetina",
        "PED",
        "Retina",
        "RNFL",
        "RNFL\\+GCL\\+IPL",
        sep = "|"
      ),
      ")"
    )
  )
  
  # Extract measurement type (always Thickness in this case)
  measure_type <- "Thickness"
  
  list(measure_type = measure_type, layer = layer)
}

# 加载baseline_info数据，获取手术日期
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 确保手术日期是日期格式
baseline_info$surgery_time_1 <- as.Date(baseline_info$surgery_time_1)

# 创建ID和手术日期的查找表
surgery_dates <- baseline_info %>%
  select(ID, surgery_time_1) %>%
  filter(!is.na(surgery_time_1))

# 创建一个数据框来存储每次检查与手术的相对时间
exam_surgery_info <- data.frame(
  id = character(),
  exam_date = as.Date(character()),
  time_point = character(),
  days_from_surgery = numeric(),
  stringsAsFactors = FALSE
)

# Initialize empty list for all results
result_list <- list()

# Main processing
all_subjects <- dir("2_data/OCTA/OCTA_thickness_3/")

for(i in 1:length(all_subjects)) {
  cat("\nProcessing subject", i, "(", all_subjects[i], ")\n")
  current_id <- paste0("SH", stringr::str_extract(all_subjects[i], "[0-9]{3}"))
  
  # 获取当前受试者的手术日期
  subject_surgery_date <- surgery_dates %>%
    filter(ID == current_id) %>%
    pull(surgery_time_1)
  
  # 如果没有手术日期数据，则跳过此患者
  if (length(subject_surgery_date) == 0 || all(is.na(subject_surgery_date))) {
    cat("No surgery date found for subject", current_id, "- skipping\n")
    next
  }
  
  # 确保获取单个有效的手术日期
  subject_surgery_date <- subject_surgery_date[1]
  
  # Get all 6x6 files for this subject
  all_files <- dir(file.path("2_data/OCTA/OCTA_thickness_3/", all_subjects[i]))
  all_files <- all_files[stringr::str_detect(all_files, "26x21")]
  
  if(length(all_files) == 0) next
  
  # 按原来的方式提取日期并分配时间点
  dates <- str_extract(all_files, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
  unique_dates <- sort(unique(dates))
  time_points <- paste0("T", seq_along(unique_dates) - 1)
  names(time_points) <- unique_dates
  
  subject_results <- list()
  
  for(temp_folder in all_files) {
    current_date_str <- str_extract(temp_folder, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    current_date <- as.Date(current_date_str)
    which_eye <- str_extract(temp_folder, "O[DS]{1}")
    current_time <- time_points[current_date_str]
    
    # 计算与手术日期的差异（天数）
    days_from_surgery <- NA
    if (!is.null(subject_surgery_date) && !is.na(subject_surgery_date)) {
      days_from_surgery <- as.numeric(difftime(current_date, subject_surgery_date, units = "days"))
    }
    
    # 添加到相对时间信息表
    exam_surgery_info <- rbind(exam_surgery_info, data.frame(
      id = current_id,
      exam_date = current_date,
      time_point = current_time,
      days_from_surgery = days_from_surgery,
      stringsAsFactors = FALSE
    ))
    
    # Get all thickness files
    quant_path <- file.path("2_data/OCTA/OCTA_thickness_3/", 
                            all_subjects[i], temp_folder, "Quantization")
    all_data <- dir(quant_path)
    all_data <- all_data[str_detect(all_data, "Thickness") & str_detect(all_data, "ETDRS")]
    
    if(length(all_data) == 0) next
    
    # Process each data file
    for(temp_data in all_data) {
      tryCatch({
        # Extract measurement info
        measure_info <- extract_measure_info(temp_data)
        if(is.na(measure_info$layer)) next
        
        # Read and process data
        data <- read_csv(
          file.path(quant_path, temp_data),
          locale = locale(encoding = "UTF-8"),
          col_names = FALSE,
          show_col_types = FALSE
        )
        
        # Process data
        result <- data[-1,] %>%
          dplyr::select(1:2) %>%
          setNames(c("region", "value")) %>%
          mutate(
            region = case_when(
              grepl("0-21", region) ~ "0_21",
              grepl("0-3", region) ~ "0_3",
              grepl("0-6", region) ~ "0_6",
              grepl("0-9", region) ~ "0_9",
              grepl("0-12", region) ~ "0_12",
              grepl("0-15", region) ~ "0-15",
              TRUE ~ NA_character_
            ),
            value = as.numeric(as.character(value))
          ) %>%
          filter(!is.na(region), !is.na(value))
        
        if (nrow(result) == 0) next
        
        # Create wide format column names with time point
        for(j in 1:nrow(result)) {
          col_name <- paste(
            measure_info$measure_type,
            measure_info$layer,
            which_eye,
            result$region[j],
            current_time,
            sep = "_"
          )
          
          # Store result with subject ID
          if(length(subject_results) == 0) {
            subject_results[[1]] <- tibble(id = current_id)
          }
          subject_results[[1]][[col_name]] <- result$value[j]
        }
        
      }, error = function(e) {
        cat("Error processing file:", e$message, "\n")
      })
    }
  }
  
  # Add subject results to main list
  if(length(subject_results) > 0) {
    result_list[[i]] <- subject_results[[1]]
  }
}

# Combine all results
final_data <- bind_rows(result_list) %>%
  arrange(.data$id)

# 分类时间点与手术的相对关系
exam_surgery_summary <- exam_surgery_info %>%
  mutate(surgery_relation = case_when(
    is.na(days_from_surgery) ~ "Unknown",
    days_from_surgery < 0 ~ "Pre-surgery",
    days_from_surgery <= 14 ~ "0-2 weeks post-surgery",
    days_from_surgery <= 45 ~ "2-6 weeks post-surgery",
    days_from_surgery <= 90 ~ "6-12 weeks post-surgery",
    TRUE ~ ">12 weeks post-surgery"
  ))

# 创建每位患者的时间点概要
time_point_summary <- exam_surgery_summary %>%
  group_by(id, time_point) %>%
  summarize(
    exam_date = dplyr::first(exam_date),
    days_from_surgery = dplyr::first(days_from_surgery),
    surgery_relation = dplyr::first(surgery_relation),
    .groups = "drop"
  ) %>%
  arrange(id, time_point)

###create working directory
if (!dir.exists(
  "2_data/analysis_data/"
)) {
  dir.create(
    "2_data/analysis_data",
    showWarnings = FALSE,
    recursive = TRUE
  )
}

setwd(
  "2_data/analysis_data/"
)

# Save results
write_csv(final_data, "octa_data_thickness_1.csv")
write_csv(time_point_summary, "octa_thickness_timepoints.csv")

# Display results
cat("\nFinal data dimensions:", nrow(final_data), "x", ncol(final_data), "\n")
cat("\nFirst few rows:\n")
print(head(final_data))

# 显示时间点与手术关系的汇总
cat("\nTime point relation to surgery:\n")
print(table(exam_surgery_summary$time_point, exam_surgery_summary$surgery_relation))

colnames(final_data)
