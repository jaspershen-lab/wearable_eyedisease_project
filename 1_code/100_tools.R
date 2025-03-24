library(tidyverse)
library(ggplot2)



# 糖尿病分组颜色
diabetes_colors <- c("#3182bd", "#e6550d")
names(diabetes_colors) <- c("Diabetes", "No Diabetes")

# 非糖尿病白内障vs糖尿病PPV分组颜色
surgery_colors <- c("#e6550d","#3182bd")
names(surgery_colors) <- c("No_Diabetes_Cataract", "Diabetes_PPV")

# 三组分组颜色
original_colors <- c("#e6550d", "#31a354", "#3182bd")
names(original_colors) <- c("No_Diabetes_Surgery0", "Diabetes_Surgery0", "Diabetes_Surgery1")

# 一致性主题设置
theme_custom <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )