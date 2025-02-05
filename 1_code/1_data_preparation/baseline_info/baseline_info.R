library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)
library(readxl)

###mental health combine
data_mental <- read_excel("2_data/mental questionnaire/mental_health_2.xlsx")
data_baseline <- read_excel("2_data/baseline questionnaire/baseline_questionnaire_2.xlsx")
data_va <- read_excel("2_data/va iop/va_iop_3.xlsx")

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

# 按 id 和 Time_to_submit_answers 排序
data_mental <- data_mental %>%
  arrange(ID, Time_to_submit_answers)

# 为每个 id 的回答添加编号
data_mental <- data_mental %>%
  group_by(ID) %>%
  mutate(answer_instance = row_number()) %>%
  ungroup()

# 修改 PHQ-9 和 GAD-7 变量名以包含回答编号
data_mental_long <- data_mental %>%
  pivot_longer(
    cols = starts_with("phq9_") | starts_with("gad7_")|starts_with("depression_diagnosed") | starts_with("anxiety_diagnosed"), 
    names_to = "question", 
    values_to = "response"
  ) %>%
  mutate(question = paste0(question, "_", answer_instance))

# 转换为宽格式
data_mental_wide <- data_mental_long %>%
  dplyr::select(-answer_instance) %>%
  pivot_wider(names_from = question, values_from = response)

write.csv(data_mental_wide, "mental_health_transformed.csv", row.names = FALSE)


### va
data_va <- data_va %>%
  mutate(across(where(is.character), ~ na_if(., ".")))


###combine with baseline
data_merged <- data_baseline %>%
  left_join(data_mental_wide, by = "ID") %>%
  left_join(data_va, by = "ID")

head(data_merged)

write.csv(data_merged, "baseline_info.csv", row.names = FALSE)
