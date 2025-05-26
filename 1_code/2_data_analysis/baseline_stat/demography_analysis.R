library(tidyverse)
library(tidymass)
library(r4projects)

library(lubridate)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
library(ggsci)

library(tableone)
library(knitr)

setwd(get_project_wd())
rm(list = ls())

# Load data
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Get unique subject IDs from heart_rate_data
heart_rate_ids <- unique(heart_rate_data@sample_info$subject_id)



# Create output directory
dir.create("3_data_analysis/2_data_analysis/baseline_stat/demography_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/2_data_analysis/baseline_stat/demography_analysis")



#####comparison table
# Add surgery type labels
baseline_info <- baseline_info %>%
  mutate(
    surgery_type = case_when(
      surgery_1..0.PI.1.other. == 0 ~ "Anterior (Cataract)",
      surgery_1..0.PI.1.other. == 1 ~ "Posterior (PPV)",
      TRUE ~ NA_character_
    )
  )

# Filter out data with NA surgery type
baseline_info_filtered <- baseline_info %>%
  filter(!is.na(surgery_type))

# Step 3: Add surgery type labels
baseline_info_filtered <- baseline_info_filtered %>%
  mutate(
    surgery_type = case_when(
      surgery_1..0.PI.1.other. == 0 ~ "Anterior (Cataract)",
      surgery_1..0.PI.1.other. == 1 ~ "Posterior (PPV)",
      TRUE ~ NA_character_
    )
  )

# Filter baseline_info to only include participants with heart rate data
baseline_info_filtered <- baseline_info %>% 
  filter(ID %in% heart_rate_ids)


# Step 5: Ensure variable types are correct and create derived variables
baseline_info_filtered <- baseline_info_filtered %>%
  mutate(
    # Convert gender to factor
    gender_factor = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
    
    # Diabetes status
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  # Has diabetes history
      diabetes_history == 2 ~ 0,  # No diabetes history
      TRUE ~ NA_real_
    ),
    dm_2 = factor(dm_2, levels = c(0, 1), labels = c("No", "Yes")),
    
    # Cataract status
    cataract_2 = case_when(
      cataract == 1 ~ 0,  # No cataract
      cataract %in% c(2, 3, 4) ~ 1, # Has cataract
      TRUE ~ NA_real_
    ),
    cataract_2 = factor(cataract_2, levels = c(0, 1), labels = c("No", "Yes")),
    
    # Ensure BMI is numeric
    bmi = as.numeric(bmi),
    
    # Ensure age is numeric
    age = as.numeric(age)
  )

# Count samples in each group
n_anterior <- sum(baseline_info_filtered$surgery_type == "Anterior (Cataract)", na.rm=TRUE)
n_posterior <- sum(baseline_info_filtered$surgery_type == "Posterior (PPV)", na.rm=TRUE)

# Print group counts
cat("\n============ SAMPLE SIZES ============\n")
cat("Total participants with heart rate data:", nrow(baseline_info_filtered), "\n")
cat("Anterior (Cataract) group:", n_anterior, "\n")
cat("Posterior (PPV) group:", n_posterior, "\n")

# Function to calculate and print statistics for continuous variables
print_continuous_stats <- function(data, var_name, var) {
  anterior_data <- data[[var]][data$surgery_type == "Anterior (Cataract)"]
  posterior_data <- data[[var]][data$surgery_type == "Posterior (PPV)"]
  
  # Calculate statistics
  anterior_mean <- mean(anterior_data, na.rm=TRUE)
  anterior_sd <- sd(anterior_data, na.rm=TRUE)
  posterior_mean <- mean(posterior_data, na.rm=TRUE)
  posterior_sd <- sd(posterior_data, na.rm=TRUE)
  
  # Calculate p-value using Wilcoxon test
  p_value <- wilcox.test(anterior_data, posterior_data)$p.value
  
  # Print results
  cat("\n============", var_name, "============\n")
  cat("Valid N:", sum(!is.na(data[[var]])), "\n")
  cat("Anterior (Cataract):", sprintf("%.2f (%.2f)", anterior_mean, anterior_sd), "\n")
  cat("Posterior (PPV):", sprintf("%.2f (%.2f)", posterior_mean, posterior_sd), "\n")
  cat("p-value:", ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), "\n")
}

# Function to calculate and print statistics for categorical variables
print_categorical_stats <- function(data, var_name, var) {
  tab <- table(data$surgery_type, data[[var]])
  props <- prop.table(tab, 1) * 100  # Row percentages
  
  # Calculate p-value using Chi-square or Fisher's exact test
  test_result <- if (min(tab) < 5 || ncol(tab) * nrow(tab) < 20) {
    fisher.test(tab)
  } else {
    chisq.test(tab)
  }
  p_value <- test_result$p.value
  
  # Print results
  cat("\n============", var_name, "============\n")
  cat("Valid N:", sum(!is.na(data[[var]])), "\n")
  cat("p-value:", ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), "\n\n")
  
  # Print counts and percentages for each level
  for (level in unique(data[[var]])) {
    cat(level, ":\n")
    anterior_count <- tab["Anterior (Cataract)", level]
    anterior_pct <- props["Anterior (Cataract)", level]
    posterior_count <- tab["Posterior (PPV)", level]
    posterior_pct <- props["Posterior (PPV)", level]
    
    cat("  Anterior (Cataract):", sprintf("%d (%.1f%%)", anterior_count, anterior_pct), "\n")
    cat("  Posterior (PPV):", sprintf("%d (%.1f%%)", posterior_count, posterior_pct), "\n\n")
  }
}

# Print statistics for each variable
print_continuous_stats(baseline_info_filtered, "AGE", "age")
print_continuous_stats(baseline_info_filtered, "BMI", "bmi")
print_categorical_stats(baseline_info_filtered, "GENDER", "gender_factor")
print_categorical_stats(baseline_info_filtered, "DIABETES", "dm_2")
print_categorical_stats(baseline_info_filtered, "CATARACT", "cataract_2")




########绘图
# # Define cataract and diabetes status
baseline_info_filtered <- baseline_info_filtered %>%
  mutate(
    cataract_2 = case_when(
      cataract == 1 ~ 0,  # No cataract
      cataract %in% c(2, 3, 4) ~ 1, # Has cataract
      TRUE ~ NA_real_
    ),
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  # Has diabetes history
      diabetes_history == 2 ~ 0,  # No diabetes history
      TRUE ~ NA_real_
    ),
    # Make sure BMI is numeric
    bmi = as.numeric(bmi)
  )

# 准备数据框
# 对 baseline_info_filtered 作相应处理，创建与 circos 兼容的数据框
df <- baseline_info_filtered %>%
  # 为每个样本分配因子水平，保持顺序
  # 使用ID列而不是sample_id
  mutate(
    factors = ID,
    x = 1,
    y = 1
  ) %>%
  # 确保因子水平的顺序与数据框一致
  mutate(factors = factor(factors, levels = factors))

# 定义颜色 - 修改为分开的疾病颜色
gender_colors <- c(Female = "#eac4d5", Male = "#95b8d1")  # 0=female(pink), 1=male(blue)
# 分别为白内障和糖尿病定义不同的颜色方案
cataract_colors <- c(No = "#fbf2c4", Yes = "#f4a582")      # 0=no(light yellow), 1=yes(orange-red)
diabetes_colors <- c(No = "#e5f5e0", Yes = "#74c476")      # 0=no(light green), 1=yes(darker green)

# 将数值型转换为颜色标签
df <- df %>%
  mutate(
    gender_label = if_else(gender == 0, "Female", "Male"),
    cataract_label = if_else(cataract_2 == 0, "No", "Yes"),
    diabetes_label = if_else(dm_2 == 0, "No", "Yes")
  )

# 设置 circos 参数
circos.par(
  "track.height" = 0.2,
  start.degree = 20,  # 从左上角开始，对应第二段代码中的 pi/9
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 20),  # 最后一个样本后添加小间隙
  cell.padding = c(0, 0, 0, 0)
)

# 初始化 circos
circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.5, 1.5))

# 绘制第一轨道：年龄
# 确保age是数值型
df$age <- as.numeric(df$age)
age_range <- range(df$age, na.rm = TRUE)
temp_value <- df$age

circos.track(
  factors = df$factors,
  y = temp_value,
  ylim = c(0.8 * min(temp_value, na.rm = TRUE), 1.1 * max(temp_value, na.rm = TRUE)),
  bg.border = "black",
  bg.col = "#FFFFFF",  # 白色背景，使线条更明显
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 仅在第一个扇区添加y轴
    if(i == 1) {
      circos.yaxis(
        side = "left",
        at = c(round(min(temp_value, na.rm = TRUE), 0), 
               round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)) / 2, 0), 
               round(max(temp_value, na.rm = TRUE), 0)),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.6,
        labels.niceFacing = FALSE
      )
    }
    
    # 使用颜色渐变表示年龄
    current_age <- temp_value[i]
    if(!is.na(current_age)) {
      # 计算归一化值用于颜色
      age_normalize <- (current_age - age_range[1]) / (age_range[2] - age_range[1])
      # 确保值在0-1之间
      age_normalize <- max(0, min(1, age_normalize))
      age_color <- colorRampPalette(c("#e0f3db", "#31a354"))(100)[ceiling(age_normalize * 99) + 1]
      
      # 绘制线条
      circos.lines(
        x = c(mean(xlim, na.rm = TRUE), mean(xlim, na.rm = TRUE)),
        y = c(ylim[1], current_age),
        type = "l",
        col = age_color,
        lwd = 3
      )
    }
  }
)

# 添加年龄图例
color_breaks <- seq(age_range[1], age_range[2], length.out = 5)
color_labels <- round(color_breaks)
pushViewport(viewport(x = 0.87, y = 0.85, width = 0.2, height = 0.1))
color_bar <- colorRampPalette(c("#e0f3db", "#31a354"))(100)
for(i in 1:100) {
  grid.rect(x = 0.1 + (i-1)*0.8/100, y = 0.5, width = 0.8/100, height = 0.3,
            gp = gpar(fill = color_bar[i], col = NA))
}
for(i in 1:5) {
  grid.text(color_labels[i], x = 0.1 + (i-1)*0.8/4, y = 0.2, just = "center", gp = gpar(fontsize = 8))
}
grid.text("Age", x = 0.5, y = 0.85, just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
popViewport()

# 绘制第二轨道：BMI
# 确保BMI是数值型
df$bmi <- as.numeric(df$bmi)
bmi_range <- range(df$bmi, na.rm = TRUE)
temp_value <- df$bmi

circos.track(
  factors = df$factors,
  y = temp_value,
  ylim = c(0.8 * min(temp_value, na.rm = TRUE), 1.1 * max(temp_value, na.rm = TRUE)),
  bg.border = "black",
  bg.col = "#FFFFFF",  # 白色背景，使线条更明显
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 仅在第一个扇区添加y轴
    if(i == 1) {
      circos.yaxis(
        side = "left",
        at = c(round(min(temp_value, na.rm = TRUE), 1), 
               round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)) / 2, 1), 
               round(max(temp_value, na.rm = TRUE), 1)),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.6,
        labels.niceFacing = FALSE
      )
    }
    
    # 使用颜色渐变表示BMI - 使用不同于年龄的配色方案
    current_bmi <- temp_value[i]
    if(!is.na(current_bmi)) {
      # 计算归一化值用于颜色
      bmi_normalize <- (current_bmi - bmi_range[1]) / (bmi_range[2] - bmi_range[1])
      # 确保值在0-1之间
      bmi_normalize <- max(0, min(1, bmi_normalize))
      # 使用蓝-紫色系表示BMI
      bmi_color <- colorRampPalette(c("#c6dbef", "#6baed6", "#2171b5"))(100)[ceiling(bmi_normalize * 99) + 1]
      
      # 绘制线条
      circos.lines(
        x = c(mean(xlim, na.rm = TRUE), mean(xlim, na.rm = TRUE)),
        y = c(ylim[1], current_bmi),
        type = "l",
        col = bmi_color,
        lwd = 3
      )
    }
  }
)

# 添加BMI图例
color_breaks <- seq(bmi_range[1], bmi_range[2], length.out = 5)
color_labels <- round(color_breaks, 1)
pushViewport(viewport(x = 0.87, y = 0.70, width = 0.2, height = 0.1))
color_bar <- colorRampPalette(c("#c6dbef", "#6baed6", "#2171b5"))(100)
for(i in 1:100) {
  grid.rect(x = 0.1 + (i-1)*0.8/100, y = 0.5, width = 0.8/100, height = 0.3,
            gp = gpar(fill = color_bar[i], col = NA))
}
for(i in 1:5) {
  grid.text(color_labels[i], x = 0.1 + (i-1)*0.8/4, y = 0.2, just = "center", gp = gpar(fontsize = 8))
}
grid.text("BMI", x = 0.5, y = 0.85, just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
popViewport()

# 绘制第三轨道：性别
circos.track(
  factors = df$factors,
  y = df$y,
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 获取性别颜色
    gender_col <- gender_colors[df$gender_label[i]]
    if(is.na(gender_col)) gender_col <- "grey"
    
    # 绘制矩形
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = gender_col,
      border = NA
    )
  }
)

# 绘制第四轨道：白内障状态 (使用新的颜色方案)
circos.track(
  factors = df$factors,
  y = df$y,
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 获取白内障颜色 - 使用专门的白内障颜色方案
    cataract_col <- cataract_colors[df$cataract_label[i]]
    if(is.na(cataract_col)) cataract_col <- "grey"
    
    # 绘制矩形
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = cataract_col,
      border = NA
    )
  }
)

# 绘制第五轨道：糖尿病状态 (使用新的颜色方案)
circos.track(
  factors = df$factors,
  y = df$y,
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 获取糖尿病颜色 - 使用专门的糖尿病颜色方案
    diabetes_col <- diabetes_colors[df$diabetes_label[i]]
    if(is.na(diabetes_col)) diabetes_col <- "grey"
    
    # 绘制矩形
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = diabetes_col,
      border = NA
    )
  }
)

# 添加中心文字 - 显示参与者数量
# 在circlize中添加中心文字需要使用grid系统
grid.text(label = nrow(df), x = 0.5, y = 0.5, 
          gp = gpar(fontsize = 24, fontface = "bold"))
grid.text(label = "Participants", x = 0.5, y = 0.43, 
          gp = gpar(fontsize = 12))

# 添加图例 - 修改为分开的疾病图例
# 创建图例的函数
draw_legend <- function() {
  # 创建图例区域
  pushViewport(viewport(x = 0.87, y = 0.4, width = 0.2, height = 0.5))
  
  # 图例标题
  grid.text("Gender:", x = 0.1, y = 0.9, just = "left", gp = gpar(fontsize = 12))
  
  # 性别图例
  grid.rect(x = 0.15, y = 0.82, width = 0.1, height = 0.05, 
            gp = gpar(fill = gender_colors["Female"], col = NA))
  grid.text("Female", x = 0.3, y = 0.82, just = "left", gp = gpar(fontsize = 10))
  
  grid.rect(x = 0.15, y = 0.74, width = 0.1, height = 0.05, 
            gp = gpar(fill = gender_colors["Male"], col = NA))
  grid.text("Male", x = 0.3, y = 0.74, just = "left", gp = gpar(fontsize = 10))
  
  # 白内障状态图例 - 独立图例
  grid.text("Cataract:", x = 0.1, y = 0.62, just = "left", gp = gpar(fontsize = 12))
  
  grid.rect(x = 0.15, y = 0.54, width = 0.1, height = 0.05, 
            gp = gpar(fill = cataract_colors["No"], col = NA))
  grid.text("No", x = 0.3, y = 0.54, just = "left", gp = gpar(fontsize = 10))
  
  grid.rect(x = 0.15, y = 0.46, width = 0.1, height = 0.05, 
            gp = gpar(fill = cataract_colors["Yes"], col = NA))
  grid.text("Yes", x = 0.3, y = 0.46, just = "left", gp = gpar(fontsize = 10))
  
  # 糖尿病状态图例 - 独立图例
  grid.text("Diabetes:", x = 0.1, y = 0.34, just = "left", gp = gpar(fontsize = 12))
  
  grid.rect(x = 0.15, y = 0.26, width = 0.1, height = 0.05, 
            gp = gpar(fill = diabetes_colors["No"], col = NA))
  grid.text("No", x = 0.3, y = 0.26, just = "left", gp = gpar(fontsize = 10))
  
  grid.rect(x = 0.15, y = 0.18, width = 0.1, height = 0.05, 
            gp = gpar(fill = diabetes_colors["Yes"], col = NA))
  grid.text("Yes", x = 0.3, y = 0.18, just = "left", gp = gpar(fontsize = 10))
  
  popViewport()
}

# 绘制图例
draw_legend()

# 重置circos参数
circos.clear()

# 保存图表
dev.copy(pdf, "circular_demography_plot_with_separate_diseases.pdf", width = 14, height = 8)
dev.off()




# 
# # ===== HEATMAP VISUALIZATION =====
# 
# # Set row names
# rownames(baseline_info_filtered) <- paste0("Sample_", 1:nrow(baseline_info_filtered))
# 
# # Define color schemes
# gender_colors <- c("0" = "#eac4d5", "1" = "#95b8d1")  # 0=female(pink), 1=male(blue)
# disease_colors <- c("0" = "#fbf2c4", "1" = "#b8e0d4")   # 0=no(light yellow), 1=yes(light green)
# 
# # Set transparency
# alpha_value = 0.7
# 
# # Create annotation object
# ha = columnAnnotation(
#   # Age barplot
#   Age = anno_barplot(baseline_info_filtered$age,
#                      gp = gpar(fill = scales::alpha("#809bce", alpha_value)),
#                      border = TRUE,
#                      width = unit(2, "cm")),
#   
#   # Gender, cataract and diabetes category annotations
#   Gender = factor(baseline_info_filtered$gender),  # Convert to factor
#   Cataract = factor(baseline_info_filtered$cataract_2),  # Convert to factor
#   Diabetes = factor(baseline_info_filtered$dm_2),  # Convert to factor
#   
#   # Set colors for categorical variables
#   col = list(
#     Gender = gender_colors,
#     Cataract = disease_colors,
#     Diabetes = disease_colors
#   ),
#   
#   # Set annotation style
#   annotation_name_gp = gpar(fontsize = 10),
#   annotation_name_side = "left",
#   simple_anno_size = unit(0.5, "cm"),
#   
#   # Add legend
#   show_legend = TRUE,
#   annotation_legend_param = list(
#     Gender = list(
#       title = "Gender",
#       at = c("0", "1"),
#       labels = c("Female", "Male")
#     ),
#     Cataract = list(
#       title = "Cataract",
#       at = c("0", "1"),
#       labels = c("No", "Yes")
#     ),
#     Diabetes = list(
#       title = "Diabetes",
#       at = c("0", "1"),
#       labels = c("No", "Yes")
#     )
#   ),
#   
#   # Set column spacing
#   gap = unit(c(2, 1, 1, 1), "mm")
# )
# 
# # Create an empty matrix for plotting heatmap
# mat = matrix(0, nrow = 1, ncol = nrow(baseline_info_filtered))
# 
# # Create and plot heatmap
# ht = Heatmap(mat,
#              top_annotation = ha,
#              show_row_names = FALSE,
#              show_column_names = FALSE,
#              show_heatmap_legend = FALSE,
#              cluster_rows = FALSE,
#              cluster_columns = FALSE,
#              height = unit(0.2, "mm"),  # Reduce height
#              col = "white")  # Set matrix color to white
# 
# # Plot and set legend position
# pdf("demography_analysis.pdf", width = 9, height = 7)
# draw(ht, annotation_legend_side = "bottom")
# dev.off()
# 
# 
# # ===== CIRCULAR VISUALIZATION =====
# 
# # Create sample IDs
# n_samples <- nrow(baseline_info_filtered)
# baseline_info_filtered$sample_id <- paste0("SF", 1500 + 1:n_samples)
# 
# # Convert data to long format for ggplot
# data_long <- baseline_info_filtered %>%
#   dplyr::select(sample_id, age, gender, cataract_2, dm_2) %>%
#   pivot_longer(cols = c(age, gender, cataract_2, dm_2),
#                names_to = "variable",
#                values_to = "value")
# 
# # Create a factor for variables to control their order (inner to outer)
# data_long$variable <- factor(data_long$variable, 
#                              levels = c("age", "gender", "cataract_2", "dm_2"),
#                              labels = c("Age", "Gender", "Cataract", "Diabetes"))
# 
# # Set angle range for the notched circle
# # Create a notch at approx 45 degrees in the top-left
# start_angle <- pi/9       # About 20 degrees, start from top-left
# end_angle <- 2*pi - pi/36 # End slightly less than 360 degrees, ensuring a notch
# 
# # Calculate angle position for each sample
# sample_ids <- unique(data_long$sample_id)
# n_samples <- length(sample_ids)
# 
# # Calculate uniform angle distribution within the specified range
# angles <- seq(start_angle, end_angle, length.out = n_samples)
# data_long$angle <- rep(angles, 4)  # Four variables
# 
# # Create a mapping for radial position
# data_long <- data_long %>%
#   mutate(
#     # Inner radius: between 1.5 and 5, increasing by variable
#     inner_radius = case_when(
#       variable == "Age" ~ 1.5,
#       variable == "Gender" ~ 2.5,
#       variable == "Cataract" ~ 3.5,
#       variable == "Diabetes" ~ 4.5
#     ),
#     # Outer radius: inner radius plus 0.8
#     outer_radius = inner_radius + 0.8
#   )
# 
# # Create color mapping function
# get_color <- function(variable, value) {
#   if (is.na(value)) return("#CCCCCC")  # Gray for NA values
#   
#   if (variable == "Age") {
#     # Continuous color gradient: age - using green gradient
#     age_range <- range(baseline_info_filtered$age, na.rm = TRUE)
#     age_normalize <- (value - age_range[1]) / (age_range[2] - age_range[1])
#     # Gradient from light green to dark green
#     return(colorRampPalette(c("#e0f3db", "#31a354"))(100)[ceiling(age_normalize * 99) + 1])
#   } else if (variable == "Gender") {
#     # Categorical color: gender (0=female pink, 1=male blue)
#     return(ifelse(value == 0, "#eac4d5", "#95b8d1"))
#   } else if (variable == "Cataract" || variable == "Diabetes") {
#     # Categorical color: disease status (0=no yellow, 1=yes green)
#     return(ifelse(value == 0, "#fbf2c4", "#b8e0d4"))
#   }
# }
# 
# # Apply color mapping
# data_long$color <- mapply(get_color, data_long$variable, data_long$value)
# 
# # Convert to polar coordinate data: calculate arc start and end points
# arc_data <- data_long %>%
#   mutate(
#     # Calculate angle width for each sample arc
#     angle_width = (end_angle - start_angle) / n_samples,
#     # Calculate start and end angles
#     start_angle = angle - angle_width/2,
#     end_angle = angle + angle_width/2
#   )
# 
# # Create point data for arcs (need multiple points to draw curves)
# points_per_arc <- 10
# polygons <- list()
# counter <- 1
# 
# for (i in 1:nrow(arc_data)) {
#   arc <- arc_data[i, ]
#   angles <- seq(arc$start_angle, arc$end_angle, length.out = points_per_arc)
#   
#   # Create polygon points
#   inner_x <- arc$inner_radius * cos(angles)
#   inner_y <- arc$inner_radius * sin(angles)
#   
#   outer_x <- arc$outer_radius * cos(rev(angles))
#   outer_y <- arc$outer_radius * sin(rev(angles))
#   
#   x <- c(inner_x, outer_x)
#   y <- c(inner_y, outer_y)
#   
#   polygons[[counter]] <- data.frame(
#     x = x,
#     y = y,
#     id = counter,
#     variable = arc$variable,
#     sample_id = arc$sample_id,
#     color = arc$color
#   )
#   
#   counter <- counter + 1
# }
# 
# # Merge all polygons
# all_polygons <- do.call(rbind, polygons)
# 
# # Add variable label positions
# var_labels <- data.frame(
#   variable = c("Age", "Gender", "Cataract", "Diabetes"),
#   radius = c(1.5, 2.5, 3.5, 4.5),
#   angle = rep(5.8, 4)  # Position in bottom right
# )
# 
# var_labels <- var_labels %>%
#   mutate(
#     x = radius * cos(angle),
#     y = radius * sin(angle)
#   )
# 
# # Create circular heatmap
# p <- ggplot() +
#   # Draw arc polygons
#   geom_polygon(data = all_polygons, 
#                aes(x = x, y = y, group = id, fill = color),
#                color = "white", size = 0.1) + # Add white border
#   scale_fill_identity() +
#   
#   # Add variable labels
#   geom_text(data = var_labels, 
#             aes(x = x, y = y, label = variable),
#             hjust = 0, size = 3.5) +
#   
#   # Add center text - number of participants with heart rate data
#   annotate("text", x = 0, y = 0.7, label = n_samples, size = 8, fontface = "bold") +
#   annotate("text", x = 0, y = -0.3, label = "Participants", size = 5) +
#   
#   # Remove default theme elements
#   theme_void() +
#   
#   # Ensure the plot is a perfect circle
#   coord_fixed() +
#   
#   # Extend plot range to fit all content
#   xlim(-6, 6) + 
#   ylim(-6, 6)
# 
# # Add legend on right side
# # Gender legend
# gender_legend <- data.frame(
#   x = rep(4.0, 2),
#   y = c(4.5, 4.0),
#   color = c("#eac4d5", "#95b8d1"),
#   label = c("Female", "Male")
# )
# 
# # Disease status legend
# disease_legend <- data.frame(
#   x = rep(4.0, 2),
#   y = c(3.0, 2.5),
#   color = c("#fbf2c4", "#b8e0d4"),
#   label = c("No", "Yes")
# )
# 
# p <- p +
#   # Add legend title
#   annotate("text", x = 4.0, y = 5.0, label = "Gender:", hjust = 0, size = 4) +
#   annotate("text", x = 4.0, y = 3.5, label = "Disease Status:", hjust = 0, size = 4) +
#   
#   # Add legend points
#   geom_point(data = gender_legend, aes(x = x, y = y, color = color), size = 4) +
#   geom_point(data = disease_legend, aes(x = x, y = y, color = color), size = 4) +
#   
#   # Add legend text
#   geom_text(data = gender_legend, aes(x = x + 0.5, y = y, label = label), hjust = 0, size = 3.5) +
#   geom_text(data = disease_legend, aes(x = x + 0.5, y = y, label = label), hjust = 0, size = 3.5) +
#   
#   scale_color_identity()
# 
# p
# 
# # Save the plot
# ggsave("circular_demography_plot.pdf", plot = p, width = 8, height = 8)






