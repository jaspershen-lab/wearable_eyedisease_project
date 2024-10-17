library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

###read data
data2 <-
  readxl::read_xlsx(
    "2_data/wearable data/",
    sheet = 1
  )


###
dir.create()
setwd()
 ##test