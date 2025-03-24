library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
library(lubridate)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
rm(list = ls())

# Load data
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Get unique subject IDs from heart_rate_data
heart_rate_ids <- unique(heart_rate_data@sample_info$subject_id)

# Filter baseline_info to only include participants with heart rate data
baseline_info_filtered <- baseline_info %>% 
  filter(ID %in% heart_rate_ids)

# Create output directory
dir.create("3_data_analysis/2_data_analysis/baseline_stat/demography_analysis", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/2_data_analysis/baseline_stat/demography_analysis")


library(tidyverse)
library(circlize)
library(ggsci)

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
    )
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

# 定义颜色
gender_colors <- c(Female = "#eac4d5", Male = "#95b8d1")  # 0=female(pink), 1=male(blue)
disease_colors <- c(No = "#fbf2c4", Yes = "#b8e0d4")      # 0=no(light yellow), 1=yes(light green)

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
      
      # # 在线条顶端添加小圆点
      # circos.points(
      #   x = mean(xlim, na.rm = TRUE),
      #   y = current_age,
      #   pch = 16,
      #   cex = 0.8,
      #   col = age_color
      # )
    }
  }
)

# 添加年龄图例
color_breaks <- seq(age_range[1], age_range[2], length.out = 5)
color_labels <- round(color_breaks)
pushViewport(viewport(x = 0.87, y = 0.8, width = 0.2, height = 0.1))
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



# 绘制第二轨道：性别
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

# 绘制第三轨道：白内障状态
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
    
    # 获取白内障颜色
    cataract_col <- disease_colors[df$cataract_label[i]]
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

# 绘制第四轨道：糖尿病状态
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
    
    # 获取糖尿病颜色
    diabetes_col <- disease_colors[df$diabetes_label[i]]
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
          gp = gpar(fontsize = 16))

# 添加图例
# 创建图例的函数
draw_legend <- function() {
  # 创建图例区域
  pushViewport(viewport(x = 0.87, y = 0.5, width = 0.2, height = 0.5))
  
  # 图例标题
  grid.text("Gender:", x = 0.1, y = 0.9, just = "left", gp = gpar(fontsize = 12))
  
  # 性别图例
  grid.rect(x = 0.15, y = 0.8, width = 0.1, height = 0.05, 
            gp = gpar(fill = gender_colors["Female"], col = NA))
  grid.text("Female", x = 0.3, y = 0.8, just = "left", gp = gpar(fontsize = 10))
  
  grid.rect(x = 0.15, y = 0.7, width = 0.1, height = 0.05, 
            gp = gpar(fill = gender_colors["Male"], col = NA))
  grid.text("Male", x = 0.3, y = 0.7, just = "left", gp = gpar(fontsize = 10))
  
  # 疾病状态图例
  grid.text("Disease Status:", x = 0.1, y = 0.5, just = "left", gp = gpar(fontsize = 12))
  
  grid.rect(x = 0.15, y = 0.4, width = 0.1, height = 0.05, 
            gp = gpar(fill = disease_colors["No"], col = NA))
  grid.text("No", x = 0.3, y = 0.4, just = "left", gp = gpar(fontsize = 10))
  
  grid.rect(x = 0.15, y = 0.3, width = 0.1, height = 0.05, 
            gp = gpar(fill = disease_colors["Yes"], col = NA))
  grid.text("Yes", x = 0.3, y = 0.3, just = "left", gp = gpar(fontsize = 10))
  
  popViewport()
}

# 绘制图例
draw_legend()

# 重置circos参数
circos.clear()

# 保存图表
dev.copy(pdf, "circular_demography_plot.pdf", width = 14, height = 8)
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






