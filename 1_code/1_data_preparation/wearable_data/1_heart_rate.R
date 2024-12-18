library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

###read data
all_subjects <- dir("2_data/OCTA/OCTA_bloodflow/")

for(i in 1:all_subjects){
  cat(i, "")
  subject_id <-
    paste0("SH",stringr::str_extract(i, "[0-9]{3}"))
    
  all_files <-
    dir(file.path("2_data/OCTA/OCTA_bloodflow/",i))

  all_files <-
    all_files[stringr::str_detect(all_files, "6x6")]
  
  all_files %>% 
  purrr::map(function(temp_folder){
    cat(temp_folder)
    date <- 
      stringr::str_extract(temp_folder, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
    which_eye <-
      stringr::str_extract(temp_folder, "O[DS]{1}")
    
    all_data <-
      dir(file.path("2_data/OCTA/OCTA_bloodflow/", i, temp_folder, "Quantization"))
    
    all_data %>% 
    purrr::map(function(temp_data){
      
      data <- read_csv(
        file.path("2_data/OCTA/OCTA_bloodflow/", i, temp_folder, "Quantization", temp_data),
        locale = locale(encoding = "UTF-8"),
        col_names = FALSE
      )
      
      data <-
        data[-1,]
      
      
      
      
    })
    
    
    
    
  })
  
    
}

