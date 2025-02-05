library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_bloodoxygen_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  colnames(data) <- c("数据时间", "平均血氧值(%)", "平均血氧值.测量时间", "外部ID")
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)

blood_oxygen_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat(subject_id, " ")
    
    hr_files <- list.files(
      folder,
      pattern = "t_apvddepp_continuousbloodoxygensaturation.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hr_files) > 0) {
      subject_data <- map_df(hr_files, read_bloodoxygen_file)
      colnames(subject_data) <-
        c("data_time", "blood_oxygen", "measure_time", "subject_id")
      subject_data <-
        subject_data %>%
        dplyr::select(subject_id, measure_time, blood_oxygen)
      
      subject_data <-
        subject_data %>%
        dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
        dplyr::mutate(measure_time = as.POSIXct(measure_time)) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else{
      return(NULL)
    }
  })


blood_oxygen_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(measure_time, blood_oxygen)) +
  geom_line()


blood_oxygen_data <-
  blood_oxygen_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  blood_oxygen_data %>%
  dplyr::select(sample_id, subject_id, measure_time)

expression_data <-
  blood_oxygen_data %>%
  dplyr::select(blood_oxygen) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

variable_info <-
  data.frame(variable_id = "blood_oxygen")

sample_info$class <- "Subject"

library(massdataset)

blood_oxygen_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

dir.create("3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen", recursive = TRUE)
save(blood_oxygen_data, file = "3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen/blood_oxygen_data.rda", compress = "xz")
