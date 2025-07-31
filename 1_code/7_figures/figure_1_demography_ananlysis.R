library(tidyverse)
library(tidymass)
library(r4projects)
library(ggstatsplot)
library(ggpie)
library(ggplot2)
library(ggsci)

setwd(get_project_wd())
rm(list = ls())

# Load data
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Get unique subject IDs from heart_rate_data
heart_rate_ids <- unique(heart_rate_data@sample_info$subject_id)

# Create output directory
dir.create("3_data_analysis/7_figures/figure1/demography_analysis", recursive = TRUE)
setwd("3_data_analysis/7_figures/figure1//demography_analysis")

# Add surgery type labels
baseline_info <- baseline_info %>%
  mutate(
    surgery_type = case_when(
      surgery_1..0.PI.1.other. == 0 ~ "Anterior (Cataract)",
      surgery_1..0.PI.1.other. == 1 ~ "Posterior (PPV)",
      TRUE ~ NA_character_
    )
  )

# Filter to only include PPV patients with diabetes
baseline_info_filtered <- baseline_info %>% 
  filter(ID %in% heart_rate_ids) %>%
  filter(surgery_type == "Posterior (PPV)") %>%
  filter(diabetes_history == 1)  # Only patients with diabetes

# Prepare data for analysis
subject_info <- baseline_info_filtered %>%
  mutate(
    # Convert gender to factor with proper labels  
    gender_factor = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
    # Ensure age is numeric
    age = as.numeric(age),
    # Ensure BMI and HbA1c are numeric
    bmi = as.numeric(bmi),
    hba1c = as.numeric(hba1c),
    # Create group variable (all are diabetic PPV, so create a single group)
    group = "Diabetic PPV",
    # Add cohort variable
    cohort = "discovery"
  )

# Define color schemes
sex_color <- c("Female" = "#eac4d5", "Male" = "#95b8d1")
group_color <- c("Diabetic PPV" = "#2E8B57")

# Define theme
theme_base <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 13)
  )

# Print sample info including BMI and HbA1c
cat("\n============ SAMPLE INFO ============\n")
cat("Total diabetic PPV participants:", nrow(subject_info), "\n")
cat("Female:", sum(subject_info$gender_factor == "Female"), "\n")
cat("Male:", sum(subject_info$gender_factor == "Male"), "\n")
cat("Age range:", min(subject_info$age, na.rm=TRUE), "-", max(subject_info$age, na.rm=TRUE), "years\n")
cat("Mean age (SD):", round(mean(subject_info$age, na.rm=TRUE), 1), 
    "(", round(sd(subject_info$age, na.rm=TRUE), 1), ")\n")
cat("BMI range:", min(subject_info$bmi, na.rm=TRUE), "-", max(subject_info$bmi, na.rm=TRUE), "\n")
cat("Mean BMI (SD):", round(mean(subject_info$bmi, na.rm=TRUE), 1), 
    "(", round(sd(subject_info$bmi, na.rm=TRUE), 1), ")\n")
cat("HbA1c range:", min(subject_info$hba1c, na.rm=TRUE), "-", max(subject_info$hba1c, na.rm=TRUE), "%\n")
cat("Mean HbA1c (SD):", round(mean(subject_info$hba1c, na.rm=TRUE), 1), 
    "(", round(sd(subject_info$hba1c, na.rm=TRUE), 1), ")\n")

##### Enhanced Circular Heatmap with BMI and HbA1c #####
library(circlize)

# Prepare data for circular plot
df_circular <-
  subject_info %>% 
  dplyr::filter(cohort == "discovery") %>%
  dplyr::arrange(age) %>%  # Sort by age for better visualization
  dplyr::mutate(
    factors = factor(ID, levels = ID),  # Use ID as factors
    x = 1,
    y = 1
  )

# Set circos parameters - adjusted for more tracks
circos.par(
  "track.height" = 0.15,  # Reduced track height to fit more tracks
  start.degree = 85,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df_circular) - 1), 15),
  cell.padding = c(0, 0, 0, 0)
)

# Initialize circos
circos.initialize(factors = df_circular$factors,
                  x = df_circular$x,
                  xlim = c(0.5, 1.5))

## Track 1: Age track with gradient colors
temp_age <- df_circular$age
circos.track(
  factors = df_circular$factors,
  y = temp_age,
  ylim = c(0.8 * min(temp_age, na.rm = TRUE), 1.1 * max(temp_age, na.rm = TRUE)),
  bg.border = "black",
  track.height = 0.15,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # Add y-axis only for the first sector
    if(i == 1) {
      circos.yaxis(
        side = "left",
        at = c(
          round(min(temp_age, na.rm = TRUE), 0),
          round((min(temp_age, na.rm = TRUE) + max(temp_age, na.rm = TRUE)) / 2, 0),
          round(max(temp_age, na.rm = TRUE), 0)
        ),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.4,
        labels.niceFacing = FALSE
      )
    }
    
    # Calculate gradient color for age
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
        cex = 8,
        type = "h",
        col = age_color,
        lwd = 4
      )
    }
  }
)

## Track 2: BMI track with gradient colors
temp_bmi <- df_circular$bmi
circos.track(
  factors = df_circular$factors,
  y = temp_bmi,
  ylim = c(0.8 * min(temp_bmi, na.rm = TRUE), 1.1 * max(temp_bmi, na.rm = TRUE)),
  bg.border = "black",
  track.height = 0.15,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # Add y-axis only for the first sector
    if(i == 1) {
      circos.yaxis(
        side = "left",
        at = c(
          round(min(temp_bmi, na.rm = TRUE), 1),
          round((min(temp_bmi, na.rm = TRUE) + max(temp_bmi, na.rm = TRUE)) / 2, 1),
          round(max(temp_bmi, na.rm = TRUE), 1)
        ),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.4,
        labels.niceFacing = FALSE
      )
    }
    
    # Calculate gradient color for BMI
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
        cex = 8,
        type = "h",
        col = bmi_color,
        lwd = 4
      )
    }
  }
)

## Track 3: HbA1c track with gradient colors
temp_hba1c <- df_circular$hba1c
circos.track(
  factors = df_circular$factors,
  y = temp_hba1c,
  ylim = c(0.8 * min(temp_hba1c, na.rm = TRUE), 1.1 * max(temp_hba1c, na.rm = TRUE)),
  bg.border = "black",
  track.height = 0.15,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # Add y-axis only for the first sector
    if(i == 1) {
      circos.yaxis(
        side = "left",
        at = c(
          round(min(temp_hba1c, na.rm = TRUE), 1),
          round((min(temp_hba1c, na.rm = TRUE) + max(temp_hba1c, na.rm = TRUE)) / 2, 1),
          round(max(temp_hba1c, na.rm = TRUE), 1)
        ),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.4,
        labels.niceFacing = FALSE
      )
    }
    
    # Calculate gradient color for HbA1c
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
        cex = 8,
        type = "h",
        col = hba1c_color,
        lwd = 4
      )
    }
  }
)

## Track 4: Gender track
temp_gender <- as.character(df_circular$gender_factor)
gender_colors_vec <- rep("grey", length(temp_gender))
gender_colors_vec[temp_gender == "Female"] <- sex_color["Female"]
gender_colors_vec[temp_gender == "Male"] <- sex_color["Male"]

circos.track(
  factors = df_circular$factors,
  y = df_circular$y,
  ylim = c(0, 1),
  bg.border = "black",
  bg.col = "white",
  track.height = 0.1,  # Slightly smaller for gender track
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

# Add center circle
draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360, 
            rou1 = 0.2, rou2 = 0, col = "white", border = "white")

# Add center labels
text(0, 0.05, nrow(df_circular), cex = 2, font = 2)
text(0, -0.05, "Diabetic PPV", cex = 1)
text(0, -0.15, "Participants", cex = 1)

# Add comprehensive legends
legend_x <- 0.78
legend_y_start <- 0.85

# Age legend
pushViewport(viewport(x = legend_x, y = legend_y_start, width = 0.15, height = 0.08))
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

# BMI legend
pushViewport(viewport(x = legend_x, y = legend_y_start - 0.15, width = 0.15, height = 0.08))
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

# HbA1c legend
pushViewport(viewport(x = legend_x, y = legend_y_start - 0.3, width = 0.15, height = 0.08))
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

# Gender legend
text(legend_x + 0.02, legend_y_start - 0.45, "Gender", cex = 1.2, font = 2, adj = 0)

# Female
rect(legend_x, legend_y_start - 0.52, legend_x + 0.03, legend_y_start - 0.49, 
     col = sex_color["Female"], border = NA)
text(legend_x + 0.05, legend_y_start - 0.505, "Female", cex = 0.8, adj = 0)

# Male
rect(legend_x, legend_y_start - 0.58, legend_x + 0.03, legend_y_start - 0.55, 
     col = sex_color["Male"], border = NA)
text(legend_x + 0.05, legend_y_start - 0.565, "Male", cex = 0.8, adj = 0)

# Clear circos
circos.clear()

# Save enhanced circular plot
dev.copy(pdf, "enhanced_circular_heatmap_diabetic_ppv.pdf", width = 12, height = 12)
dev.off()

cat("\n============ ENHANCED VISUALIZATION COMPLETE ============\n")
cat("Generated enhanced circular heatmap with Age, BMI, HbA1c, and Gender tracks\n")