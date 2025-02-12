library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)
library(readxl)

###mental health combine
data_mental <- read_excel("2_data/mental questionnaire/mental_health_2.xlsx")
data_baseline <- read_excel("2_data/baseline questionnaire/baseline_questionnaire_3.xlsx")
data_va <- read_excel("2_data/va iop/va_iop_4.xlsx")

###create working directory
if (!dir.exists(
  "2_data/analysis_data"
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

# 1. 按ID和时间排序
data_mental <- data_mental %>%
  arrange(ID, Time_to_submit_answers)

# 2. 为每个ID的回答添加编号
data_mental <- data_mental %>%
  group_by(ID) %>%
  mutate(visit_number = row_number()) %>%
  ungroup()

# 3. 转换为长格式并添加访问编号
data_mental_long <- data_mental %>%
  pivot_longer(
    cols = starts_with("phq9_") | starts_with("gad7_") |
      starts_with("depression_diagnosed") | starts_with("anxiety_diagnosed"),
    names_to = "question",
    values_to = "response"
  ) %>%
  mutate(question = paste0(question, "_visit", visit_number))

# 4. 转换回宽格式，只用ID作为标识符
data_mental_wide <- data_mental_long %>%
  dplyr::select(-visit_number) %>%
  pivot_wider(
    id_cols = ID,  # 只使用ID作为标识符
    names_from = question,
    values_from = response
  )
  

write.csv(data_mental_wide, "mental_health_transformed.csv", row.names = FALSE)


### va
data_va <- data_va %>%
  mutate(across(where(is.character), ~ na_if(., ".")))


###combine with baseline
data_merged <- data_baseline %>%
  left_join(data_mental_wide, by = "ID") %>%
  left_join(data_va, by = "ID")

head(data_merged)

names(data_merged)[names(data_merged) == "surger_eye_1"] <- "surgery_eye_1"
names(data_merged)[names(data_merged) == "surger_eye_2"] <- "surgery_eye_2"
names(data_merged)[names(data_merged) == "gender.x"] <- "gender"
names(data_merged)[names(data_merged) == "age.x"] <- "age"


write.csv(data_merged, "baseline_info.csv", row.names = FALSE)

