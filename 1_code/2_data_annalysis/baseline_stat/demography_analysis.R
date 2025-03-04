library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)

library(ComplexHeatmap)
library(circlize)
library(grid)
library(tidyverse)

# 读取数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")



#############
dir.create("3_data_analysis/2_data_analysis/baseline_stat/demography_analysis", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/baseline_stat/demography_analysis")


# 定义白内障和糖尿病状态
baseline_info <- baseline_info %>%
  mutate(
    cataract_2 = case_when(
      cataract == 1 ~ 0,  # 无白内障
      cataract %in% c(2, 3, 4) ~ 1, # 有白内障
      TRUE ~ NA_real_
    ),
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  # 有糖尿病史
      diabetes_history == 2 ~ 0,  # 无糖尿病史
      TRUE ~ NA_real_
    )
  )


# 设置行名
rownames(baseline_info) <- paste0("Sample_", 1:nrow(baseline_info))

# 定义颜色方案
gender_colors <- c("0" = "#eac4d5", "1" = "#95b8d1")  # 0=女性(粉色), 1=男性(蓝色)
disease_colors <- c("0" = "#fbf2c4", "1" = "#b8e0d4")   # 0=无(白色), 1=有(橙色)

# 设置透明度
alpha_value = 0.7

# 创建注释对象
ha = columnAnnotation(
  # 年龄条形图
  Age = anno_barplot(baseline_info$age,
                     gp = gpar(fill = scales::alpha("#809bce", alpha_value)),
                     border = TRUE,
                     width = unit(2, "cm")),
  
  # 性别、白内障和糖网的分类注释
  Gender = factor(baseline_info$gender),  # 转换为因子
  Cataract = factor(baseline_info$cataract_2),  # 转换为因子
  Diabetes = factor(baseline_info$dm_2),  # 转换为因子
  
  # 设置分类变量的颜色
  col = list(
    Gender = gender_colors,
    Cataract = disease_colors,
    Diabetes = disease_colors
  ),
  
  # 设置注释的样式
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm"),
  
  # 添加图例
  show_legend = TRUE,
  annotation_legend_param = list(
    Gender = list(
      title = "Gender",
      at = c("0", "1"),
      labels = c("Female", "Male")
    ),
    Cataract = list(
      title = "Cataract",
      at = c("0", "1"),
      labels = c("No", "Yes")
    ),
    Diabetes = list(
      title = "Diabetes",
      at = c("0", "1"),
      labels = c("No", "Yes")
    )
  ),
  
  # 设置列间距
  gap = unit(c(2, 1, 1, 1), "mm")
)

# 创建一个空矩阵用于绘制热图
mat = matrix(0, nrow = 1, ncol = nrow(baseline_info))

# 创建并绘制热图
ht = Heatmap(mat,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             show_heatmap_legend = FALSE,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             height = unit(0.2, "mm"),  # 减小高度
             col = "white")  # 将矩阵颜色设为白色

# 绘制并设置图例位置
pdf("demography_analysis.pdf", width = 9, height = 7)
draw(ht, annotation_legend_side = "bottom")
dev.off()
