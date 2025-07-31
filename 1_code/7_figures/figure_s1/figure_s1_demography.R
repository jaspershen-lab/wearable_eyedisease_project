# ================================================================================
# OCTA聚类人群和可穿戴设备聚类人群的年龄性别BMI和HbA1c分布图
# 基于代码三的风格，为两种聚类结果分别创建增强版人口统计学可视化
# ================================================================================

library(tidyverse)
library(ggstatsplot)
library(ggpie)
library(ggplot2)
library(ggsci)
library(patchwork)
library(circlize)
library(r4projects)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# ================== 1. 数据加载和准备 ==================

# 加载基础信息
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 加载OCTA聚类结果（代码一）
octa_clusters <- read.csv("3_data_analysis/6_clustering_modeling/mfuzz/WF_only_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv")

# 加载可穿戴设备聚类结果（代码二）
wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/time_window_clustering/time_window_2cluster_membership_data.csv")

# 创建输出目录
dir.create("3_data_analysis/7_figures/figure_s1/clustering_demographics", recursive = TRUE)
setwd("3_data_analysis/7_figures/figure_s1/clustering_demographics")

# ================== 2. 数据处理函数 ==================

# 准备OCTA聚类人群数据
prepare_octa_demographics <- function(baseline_info, octa_clusters) {
  
  cat("准备OCTA聚类人群数据...\n")
  
  # 合并基础信息和聚类结果
  octa_demo <- baseline_info %>%
    inner_join(octa_clusters, by = c("ID" = "subject_id")) %>%
    filter(!is.na(max_cluster)) %>%
    mutate(
      # 转换性别为因子
      gender_factor = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
      # 确保数值型变量
      age = as.numeric(age),
      bmi = as.numeric(bmi),
      hba1c = as.numeric(hba1c),
      # 添加队列标识
      cohort = "OCTA Analysis"
    ) %>%
    # 按年龄排序以便环形图展示
    arrange(age)
  
  cat("OCTA聚类人群统计:\n")
  cat("总人数:", nrow(octa_demo), "\n")
  cat("性别分布:\n")
  print(table(octa_demo$gender_factor))
  cat("年龄范围:", min(octa_demo$age, na.rm=TRUE), "-", max(octa_demo$age, na.rm=TRUE), "岁\n")
  cat("平均年龄:", round(mean(octa_demo$age, na.rm=TRUE), 1), "±", round(sd(octa_demo$age, na.rm=TRUE), 1), "\n")
  cat("BMI范围:", min(octa_demo$bmi, na.rm=TRUE), "-", max(octa_demo$bmi, na.rm=TRUE), "\n")
  cat("平均BMI:", round(mean(octa_demo$bmi, na.rm=TRUE), 1), "±", round(sd(octa_demo$bmi, na.rm=TRUE), 1), "\n")
  cat("HbA1c范围:", min(octa_demo$hba1c, na.rm=TRUE), "-", max(octa_demo$hba1c, na.rm=TRUE), "%\n")
  cat("平均HbA1c:", round(mean(octa_demo$hba1c, na.rm=TRUE), 1), "±", round(sd(octa_demo$hba1c, na.rm=TRUE), 1), "\n\n")
  
  return(octa_demo)
}

# 准备可穿戴设备聚类人群数据
prepare_wearable_demographics <- function(baseline_info, wearable_clusters) {
  
  cat("准备可穿戴设备聚类人群数据...\n")
  
  # 选择一个代表性时间窗口的聚类结果
  cluster_cols <- names(wearable_clusters)[grep("^cluster_", names(wearable_clusters))]
  
  if("cluster_early_recovery" %in% cluster_cols) {
    selected_window <- "early_recovery"
    cluster_col <- "cluster_early_recovery"
  } else if("cluster_baseline" %in% cluster_cols) {
    selected_window <- "baseline"
    cluster_col <- "cluster_baseline"
  } else {
    cluster_col <- cluster_cols[1]
    selected_window <- gsub("cluster_", "", cluster_col)
  }
  
  cat("使用时间窗口:", selected_window, "\n")
  
  # 合并基础信息和聚类结果
  wearable_demo <- baseline_info %>%
    inner_join(wearable_clusters, by = c("ID" = "subject_id")) %>%
    filter(!is.na(.data[[cluster_col]])) %>%
    mutate(
      # 转换性别为因子
      gender_factor = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
      # 确保数值型变量
      age = as.numeric(age),
      bmi = as.numeric(bmi),
      hba1c = as.numeric(hba1c),
      # 添加队列标识
      cohort = "Wearable Analysis"
    ) %>%
    # 按年龄排序以便环形图展示
    arrange(age)
  
  cat("可穿戴设备聚类人群统计:\n")
  cat("总人数:", nrow(wearable_demo), "\n")
  cat("性别分布:\n")
  print(table(wearable_demo$gender_factor))
  cat("年龄范围:", min(wearable_demo$age, na.rm=TRUE), "-", max(wearable_demo$age, na.rm=TRUE), "岁\n")
  cat("平均年龄:", round(mean(wearable_demo$age, na.rm=TRUE), 1), "±", round(sd(wearable_demo$age, na.rm=TRUE), 1), "\n")
  cat("BMI范围:", min(wearable_demo$bmi, na.rm=TRUE), "-", max(wearable_demo$bmi, na.rm=TRUE), "\n")
  cat("平均BMI:", round(mean(wearable_demo$bmi, na.rm=TRUE), 1), "±", round(sd(wearable_demo$bmi, na.rm=TRUE), 1), "\n")
  cat("HbA1c范围:", min(wearable_demo$hba1c, na.rm=TRUE), "-", max(wearable_demo$hba1c, na.rm=TRUE), "%\n")
  cat("平均HbA1c:", round(mean(wearable_demo$hba1c, na.rm=TRUE), 1), "±", round(sd(wearable_demo$hba1c, na.rm=TRUE), 1), "\n\n")
  
  return(wearable_demo)
}

# ================== 3. 颜色方案定义 ==================

define_color_schemes <- function() {
  
  # 性别颜色（保持与代码三一致）
  sex_color <- c("Female" = "#eac4d5", "Male" = "#95b8d1")
  
  return(list(
    sex = sex_color
  ))
}

# ================== 4. 增强版环形热图创建函数 ==================

create_enhanced_circular_heatmap <- function(demo_data, data_type, colors) {
  
  cat(sprintf("创建%s增强版环形热图（含BMI和HbA1c）...\n", data_type))
  
  # 准备环形图数据
  df_circular <- demo_data %>% 
    arrange(age) %>%  # 按年龄排序
    mutate(
      factors = factor(ID, levels = ID),  # 使用ID作为因子
      x = 1,
      y = 1
    )
  
  # 创建新的图形设备
  pdf(paste0(tolower(gsub(" ", "_", data_type)), "_enhanced_circular_heatmap.pdf"), width = 12, height = 12)
  
  # 设置circos参数 - 调整以适应更多轨道
  circos.par(
    "track.height" = 0.12,  # 减小轨道高度以容纳更多轨道
    start.degree = 85,
    clock.wise = TRUE,
    gap.after = c(rep(0, nrow(df_circular) - 1), 15),
    cell.padding = c(0, 0, 0, 0)
  )
  
  # 初始化circos
  circos.initialize(factors = df_circular$factors,
                    x = df_circular$x,
                    xlim = c(0.5, 1.5))
  
  ## 轨道1: 年龄轨道（与代码三完全相同的渐变色）
  temp_age <- df_circular$age
  
  circos.track(
    factors = df_circular$factors,
    y = temp_age,
    ylim = c(0.8 * min(temp_age, na.rm = TRUE), 1.1 * max(temp_age, na.rm = TRUE)),
    bg.border = "black",
    track.height = 0.12,
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      
      # 添加y轴（仅第一个扇区）
      if(i == 1) {
        circos.yaxis(
          side = "left",
          at = c(
            round(min(temp_age, na.rm = TRUE), 0),
            round((min(temp_age, na.rm = TRUE) + max(temp_age, na.rm = TRUE)) / 2, 0),
            round(max(temp_age, na.rm = TRUE), 0)
          ),
          sector.index = get.all.sector.index()[1],
          labels.cex = 0.35,
          labels.niceFacing = FALSE
        )
      }
      
      # 计算年龄渐变色（完全基于代码三的#769f4a色彩）
      current_age <- temp_age[i]
      if(!is.na(current_age)) {
        age_min <- min(temp_age, na.rm = TRUE)
        age_max <- max(temp_age, na.rm = TRUE)
        age_normalize <- (current_age - age_min) / (age_max - age_min)
        age_normalize <- max(0, min(1, age_normalize))
        age_color <- colorRampPalette(c("#c8d6b0", "#769f4a", "#5a7836"))(100)[ceiling(age_normalize * 99) + 1]
        
        circos.lines(
          x = mean(xlim, na.rm = TRUE),
          y = temp_age[i],
          pch = 16,
          cex = 6,  # 略微减小以适应更多轨道
          type = "h",
          col = age_color,
          lwd = 3
        )
      }
    }
  )
  
  ## 轨道2: BMI轨道
  temp_bmi <- df_circular$bmi
  
  circos.track(
    factors = df_circular$factors,
    y = temp_bmi,
    ylim = c(0.8 * min(temp_bmi, na.rm = TRUE), 1.1 * max(temp_bmi, na.rm = TRUE)),
    bg.border = "black",
    track.height = 0.12,
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      
      # 添加y轴（仅第一个扇区）
      if(i == 1) {
        circos.yaxis(
          side = "left",
          at = c(
            round(min(temp_bmi, na.rm = TRUE), 1),
            round((min(temp_bmi, na.rm = TRUE) + max(temp_bmi, na.rm = TRUE)) / 2, 1),
            round(max(temp_bmi, na.rm = TRUE), 1)
          ),
          sector.index = get.all.sector.index()[1],
          labels.cex = 0.35,
          labels.niceFacing = FALSE
        )
      }
      
      # 计算BMI渐变色
      current_bmi <- temp_bmi[i]
      if(!is.na(current_bmi)) {
        bmi_min <- min(temp_bmi, na.rm = TRUE)
        bmi_max <- max(temp_bmi, na.rm = TRUE)
        bmi_normalize <- (current_bmi - bmi_min) / (bmi_max - bmi_min)
        bmi_normalize <- max(0, min(1, bmi_normalize))
        bmi_color <- colorRampPalette(c("#ffcccc", "#ff6666", "#cc0000"))(100)[ceiling(bmi_normalize * 99) + 1]
        
        circos.lines(
          x = mean(xlim, na.rm = TRUE),
          y = temp_bmi[i],
          pch = 16,
          cex = 6,
          type = "h",
          col = bmi_color,
          lwd = 3
        )
      }
    }
  )
  
  ## 轨道3: HbA1c轨道
  temp_hba1c <- df_circular$hba1c
  
  circos.track(
    factors = df_circular$factors,
    y = temp_hba1c,
    ylim = c(0.8 * min(temp_hba1c, na.rm = TRUE), 1.1 * max(temp_hba1c, na.rm = TRUE)),
    bg.border = "black",
    track.height = 0.12,
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      
      # 添加y轴（仅第一个扇区）
      if(i == 1) {
        circos.yaxis(
          side = "left",
          at = c(
            round(min(temp_hba1c, na.rm = TRUE), 1),
            round((min(temp_hba1c, na.rm = TRUE) + max(temp_hba1c, na.rm = TRUE)) / 2, 1),
            round(max(temp_hba1c, na.rm = TRUE), 1)
          ),
          sector.index = get.all.sector.index()[1],
          labels.cex = 0.35,
          labels.niceFacing = FALSE
        )
      }
      
      # 计算HbA1c渐变色
      current_hba1c <- temp_hba1c[i]
      if(!is.na(current_hba1c)) {
        hba1c_min <- min(temp_hba1c, na.rm = TRUE)
        hba1c_max <- max(temp_hba1c, na.rm = TRUE)
        hba1c_normalize <- (current_hba1c - hba1c_min) / (hba1c_max - hba1c_min)
        hba1c_normalize <- max(0, min(1, hba1c_normalize))
        hba1c_color <- colorRampPalette(c("#ccccff", "#6666ff", "#0000cc"))(100)[ceiling(hba1c_normalize * 99) + 1]
        
        circos.lines(
          x = mean(xlim, na.rm = TRUE),
          y = temp_hba1c[i],
          pch = 16,
          cex = 6,
          type = "h",
          col = hba1c_color,
          lwd = 3
        )
      }
    }
  )
  
  ## 轨道4: 性别轨道
  temp_gender <- as.character(df_circular$gender_factor)
  gender_colors_vec <- rep("grey", length(temp_gender))
  gender_colors_vec[temp_gender == "Female"] <- colors$sex["Female"]
  gender_colors_vec[temp_gender == "Male"] <- colors$sex["Male"]
  
  circos.track(
    factors = df_circular$factors,
    y = df_circular$y,
    ylim = c(0, 1),
    bg.border = "black",
    bg.col = "white",
    track.height = 0.08,  # 性别轨道稍小
    panel.fun = function(x, y) {
      name = get.cell.meta.data("sector.index")
      i = get.cell.meta.data("sector.numeric.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      
      circos.rect(
        xleft = xlim[1],
        ybottom = ylim[1],
        xright = xlim[2],
        ytop = ylim[2],
        col = gender_colors_vec[i],
        border = "white",
        lwd = 0.5
      )
    }
  )
  
  # 添加中心白色圆形
  draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360, 
              rou1 = 0.2, rou2 = 0, col = "white", border = "white")
  
  # 添加中心标签
  text(0, 0.05, nrow(df_circular), cex = 2, font = 2)
  text(0, -0.05, data_type, cex = 1)
  text(0, -0.15, "Participants", cex = 1)
  
  # 添加综合图例系统
  legend_x <- 0.75
  legend_y_start <- 0.9
  
  # 年龄图例
  pushViewport(viewport(x = legend_x, y = legend_y_start, width = 0.18, height = 0.08))
  age_color_bar <- colorRampPalette(c("#c8d6b0", "#769f4a", "#5a7836"))(100)
  for(i in 1:100) {
    grid.rect(x = 0.1 + (i-1)*0.8/100, y = 0.5, width = 0.8/100, height = 0.3,
              gp = gpar(fill = age_color_bar[i], col = NA))
  }
  age_breaks <- seq(min(df_circular$age, na.rm = TRUE), max(df_circular$age, na.rm = TRUE), length.out = 3)
  for(i in 1:3) {
    grid.text(round(age_breaks[i]), x = 0.1 + (i-1)*0.8/2, y = 0.2, just = "center", gp = gpar(fontsize = 8))
  }
  grid.text("Age", x = 0.5, y = 0.85, just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
  popViewport()
  
  # BMI图例
  pushViewport(viewport(x = legend_x, y = legend_y_start - 0.15, width = 0.18, height = 0.08))
  bmi_color_bar <- colorRampPalette(c("#ffcccc", "#ff6666", "#cc0000"))(100)
  for(i in 1:100) {
    grid.rect(x = 0.1 + (i-1)*0.8/100, y = 0.5, width = 0.8/100, height = 0.3,
              gp = gpar(fill = bmi_color_bar[i], col = NA))
  }
  bmi_breaks <- seq(min(df_circular$bmi, na.rm = TRUE), max(df_circular$bmi, na.rm = TRUE), length.out = 3)
  for(i in 1:3) {
    grid.text(round(bmi_breaks[i], 1), x = 0.1 + (i-1)*0.8/2, y = 0.2, just = "center", gp = gpar(fontsize = 8))
  }
  grid.text("BMI", x = 0.5, y = 0.85, just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
  popViewport()
  
  # HbA1c图例
  pushViewport(viewport(x = legend_x, y = legend_y_start - 0.3, width = 0.18, height = 0.08))
  hba1c_color_bar <- colorRampPalette(c("#ccccff", "#6666ff", "#0000cc"))(100)
  for(i in 1:100) {
    grid.rect(x = 0.1 + (i-1)*0.8/100, y = 0.5, width = 0.8/100, height = 0.3,
              gp = gpar(fill = hba1c_color_bar[i], col = NA))
  }
  hba1c_breaks <- seq(min(df_circular$hba1c, na.rm = TRUE), max(df_circular$hba1c, na.rm = TRUE), length.out = 3)
  for(i in 1:3) {
    grid.text(round(hba1c_breaks[i], 1), x = 0.1 + (i-1)*0.8/2, y = 0.2, just = "center", gp = gpar(fontsize = 8))
  }
  grid.text("HbA1c (%)", x = 0.5, y = 0.85, just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
  popViewport()
  
  # 性别图例
  text(legend_x + 0.02, legend_y_start - 0.45, "Gender", cex = 1.2, font = 2, adj = 0)
  
  # Female
  rect(legend_x, legend_y_start - 0.52, legend_x + 0.03, legend_y_start - 0.49, 
       col = colors$sex["Female"], border = NA)
  text(legend_x + 0.05, legend_y_start - 0.505, "Female", cex = 0.8, adj = 0)
  
  # Male
  rect(legend_x, legend_y_start - 0.58, legend_x + 0.03, legend_y_start - 0.55, 
       col = colors$sex["Male"], border = NA)
  text(legend_x + 0.05, legend_y_start - 0.565, "Male", cex = 0.8, adj = 0)
  
  # 清除circos
  circos.clear()
  
  # 关闭图形设备
  dev.off()
  
  cat(sprintf("%s增强版环形热图创建完成！\n", data_type))
  
  return(paste0(tolower(gsub(" ", "_", data_type)), "_enhanced_circular_heatmap.pdf"))
}

# ================== 5. 增强统计总结函数 ==================

create_enhanced_demographic_summary <- function(demo_data, cohort_name) {
  
  # 基本统计（包含BMI和HbA1c）
  summary_stats <- demo_data %>%
    summarise(
      cohort = cohort_name,
      n_total = n(),
      n_female = sum(gender_factor == "Female"),
      n_male = sum(gender_factor == "Male"),
      female_percent = round(sum(gender_factor == "Female") / n() * 100, 1),
      male_percent = round(sum(gender_factor == "Male") / n() * 100, 1),
      age_mean = round(mean(age, na.rm = TRUE), 1),
      age_sd = round(sd(age, na.rm = TRUE), 1),
      age_median = round(median(age, na.rm = TRUE), 1),
      age_min = min(age, na.rm = TRUE),
      age_max = max(age, na.rm = TRUE),
      bmi_mean = round(mean(bmi, na.rm = TRUE), 1),
      bmi_sd = round(sd(bmi, na.rm = TRUE), 1),
      bmi_median = round(median(bmi, na.rm = TRUE), 1),
      bmi_min = round(min(bmi, na.rm = TRUE), 1),
      bmi_max = round(max(bmi, na.rm = TRUE), 1),
      hba1c_mean = round(mean(hba1c, na.rm = TRUE), 1),
      hba1c_sd = round(sd(hba1c, na.rm = TRUE), 1),
      hba1c_median = round(median(hba1c, na.rm = TRUE), 1),
      hba1c_min = round(min(hba1c, na.rm = TRUE), 1),
      hba1c_max = round(max(hba1c, na.rm = TRUE), 1)
    )
  
  return(summary_stats)
}

# ================== 6. 主执行函数 ==================

main_enhanced_circular_heatmaps <- function() {
  
  cat("========================================\n")
  cat("开始创建增强版聚类人群环形热图\n")
  cat("（包含年龄、性别、BMI和HbA1c）\n")
  cat("========================================\n")
  
  # 1. 数据准备
  cat("1. 数据加载和准备...\n")
  octa_demo <- prepare_octa_demographics(baseline_info, octa_clusters)
  wearable_demo <- prepare_wearable_demographics(baseline_info, wearable_clusters)
  
  # 2. 定义颜色方案
  colors <- define_color_schemes()
  
  # 3. 创建OCTA增强版环形热图
  cat("2. 创建OCTA聚类人群增强版环形热图...\n")
  octa_file <- create_enhanced_circular_heatmap(octa_demo, "OCTA Cohort", colors)
  
  # 4. 创建可穿戴设备增强版环形热图
  cat("3. 创建可穿戴设备聚类人群增强版环形热图...\n")
  wearable_file <- create_enhanced_circular_heatmap(wearable_demo, "Wearable Cohort", colors)
  
  # 5. 生成增强统计总结
  cat("4. 生成增强统计总结...\n")
  octa_stats <- create_enhanced_demographic_summary(octa_demo, "OCTA Cohort")
  wearable_stats <- create_enhanced_demographic_summary(wearable_demo, "Wearable Cohort")
  
  # 合并统计数据
  combined_stats <- bind_rows(octa_stats, wearable_stats)
  
  # 保存统计总结
  write.csv(combined_stats, "enhanced_demographic_summary_statistics.csv", row.names = FALSE)
  
  cat("========================================\n")
  cat("增强版聚类人群环形热图创建完成！\n")
  cat("========================================\n")
  
  # 打印总结信息
  cat("\n📊 生成的增强版环形热图文件:\n")
  cat("1. OCTA聚类人群:\n")
  cat("   -", octa_file, "\n")
  cat("2. 可穿戴设备聚类人群:\n")
  cat("   -", wearable_file, "\n")
  cat("3. 统计数据:\n")
  cat("   - enhanced_demographic_summary_statistics.csv\n\n")
  
  cat("🎨 增强版环形热图特点:\n")
  cat("✅ 年龄轨道：基于#769f4a的渐变色（与代码三一致）\n")
  cat("✅ BMI轨道：红色渐变色系统\n")
  cat("✅ HbA1c轨道：蓝色渐变色系统\n")
  cat("✅ 性别轨道：粉色(女性) + 蓝色(男性)\n")
  cat("✅ 按年龄排序：便于观察年龄分布模式\n")
  cat("✅ 中心标签：显示总人数和队列名称\n")
  cat("✅ 完整图例：所有4个变量的渐变条和色块\n\n")
  
  cat("📈 数据概览:\n")
  print(combined_stats)
  
  # 返回结果
  return(list(
    octa_demo = octa_demo,
    wearable_demo = wearable_demo,
    octa_file = octa_file,
    wearable_file = wearable_file,
    combined_stats = combined_stats
  ))
}

# ================== 7. 执行主函数 ==================
results <- main_enhanced_circular_heatmaps()
