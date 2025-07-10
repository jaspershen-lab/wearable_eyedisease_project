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

# Print sample info
cat("\n============ SAMPLE INFO ============\n")
cat("Total diabetic PPV participants:", nrow(subject_info), "\n")
cat("Female:", sum(subject_info$gender_factor == "Female"), "\n")
cat("Male:", sum(subject_info$gender_factor == "Male"), "\n")
cat("Age range:", min(subject_info$age, na.rm=TRUE), "-", max(subject_info$age, na.rm=TRUE), "years\n")
cat("Mean age (SD):", round(mean(subject_info$age, na.rm=TRUE), 1), 
    "(", round(sd(subject_info$age, na.rm=TRUE), 1), ")\n")

##### Age Analysis #####
# Since we only have one group (diabetic PPV), we'll show age distribution by gender
plot_age <-
  ggbetweenstats(
    data  = subject_info %>%
      dplyr::filter(cohort == "discovery"),
    x  = gender_factor,  # Use gender_factor instead of sex
    y  = age,
    type = "nonparametric",
    pairwise.comparisons = TRUE,
    pairwise.display = "all",
    point.args = list(
      alpha = 0.9,
      size = 3,
      position = position_jitter(width = 0.1)
    ),
    title = "Age Distribution by Gender",
    subtitle = paste("Diabetic PPV Patients (N =", nrow(subject_info), ")")
  ) +
  theme_base +
  labs(x = "Gender", y = "Age (years)") +
  theme(legend.position = "none") +
  scale_color_manual(values = sex_color)

print(plot_age)

# Save age plot
ggsave("age_analysis_ggstatsplot.pdf", plot_age, width = 8, height = 6)

##### Gender Distribution #####
# Create donut chart for gender distribution
plot_gender <-
  subject_info %>%
  dplyr::filter(cohort == "discovery") %>%
  ggdonut(
    group_key = "gender_factor",  # Use gender_factor instead of sex
    count_type = "full",
    label_info = "all",
    label_type = "horizon",
    label_split = NULL,
    label_size = 4,
    label_pos = "in",
    label_threshold = 0  # Show all labels regardless of percentage
  ) +
  scale_fill_manual(values = sex_color) +
  labs(
    title = "Gender Distribution",
    subtitle = paste("Diabetic PPV Patients (N =", nrow(subject_info), ")")
  ) +
  theme_base +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(plot_gender)

# Save gender plot
ggsave("gender_distribution_donut.pdf", plot_gender, width = 6, height = 6)

##### Alternative: Age distribution plot (since we have only one group) #####
# Create a more detailed age analysis
plot_age_dist <-
  gghistostats(
    data = subject_info %>% dplyr::filter(cohort == "discovery"),
    x = age,
    type = "parametric",
    test.value = mean(subject_info$age, na.rm = TRUE),
    title = "Age Distribution Analysis",
    subtitle = paste("Diabetic PPV Patients (N =", nrow(subject_info), ")"),
    caption = "Vertical line shows mean age"
  ) +
  theme_base +
  labs(x = "Age (years)", y = "Frequency")

print(plot_age_dist)

# Save age distribution plot
ggsave("age_distribution_detailed.pdf", plot_age_dist, width = 8, height = 6)

##### Combined visualization #####
library(patchwork)

# Create a combined plot
combined_plot <- plot_age | plot_gender

# Save combined plot
ggsave("combined_age_gender_analysis.pdf", combined_plot, width = 14, height = 6)

##### Circular Heatmap (Figure1_b style) #####
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

# Set circos parameters
circos.par(
  "track.height" = 0.2,
  start.degree = 85,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df_circular) - 1), 15),  # Increase gap from 10 to 45 degrees
  cell.padding = c(0, 0, 0, 0)
)

# Initialize circos
circos.initialize(factors = df_circular$factors,
                  x = df_circular$x,
                  xlim = c(0.5, 1.5))

## Age track with gradient colors
range(df_circular$age, na.rm = TRUE)
temp_value <- df_circular$age

circos.track(
  factors = df_circular$factors,
  y = temp_value,
  ylim = c(0.8 * min(temp_value, na.rm = TRUE), 1.1 * max(temp_value, na.rm = TRUE)),
  bg.border = "black",
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # Add y-axis only for the first sector with integer labels
    if(i == 1) {
      circos.yaxis(
        side = "left",
        at = c(
          round(min(temp_value, na.rm = TRUE), 0),
          round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)) / 2, 0),
          round(max(temp_value, na.rm = TRUE), 0)
        ),
        sector.index = get.all.sector.index()[1],
        labels.cex = 0.4,
        labels.niceFacing = FALSE
      )
    }
    
    # Calculate gradient color for each patient based on their age
    current_age <- temp_value[i]
    if(!is.na(current_age)) {
      # Calculate age range
      age_min <- min(temp_value, na.rm = TRUE)
      age_max <- max(temp_value, na.rm = TRUE)
      
      # Normalize age to 0-1 scale
      age_normalize <- (current_age - age_min) / (age_max - age_min)
      age_normalize <- max(0, min(1, age_normalize))
      
      # Create gradient color based on #769f4a
      age_color <- colorRampPalette(c("#c8d6b0", "#769f4a", "#5a7836"))(100)[ceiling(age_normalize * 99) + 1]
      
      # Draw line with gradient color (same style as reference code)
      circos.lines(
        x = mean(xlim, na.rm = TRUE),
        y = temp_value[i],
        pch = 16,
        cex = 8,
        type = "h",
        col = age_color,  # Use gradient color instead of fixed ggsci color
        lwd = 4
      )
    }
  }
)

# 添加年龄图例 - 使用基于#769f4a的渐变色
color_breaks <- seq(min(df_circular$age, na.rm = TRUE), max(df_circular$age, na.rm = TRUE), length.out = 5)
color_labels <- round(color_breaks)
pushViewport(viewport(x = 0.87, y = 0.85, width = 0.2, height = 0.1))
color_bar <- colorRampPalette(c("#c8d6b0", "#769f4a", "#5a7836"))(100)
for(i in 1:100) {
  grid.rect(x = 0.1 + (i-1)*0.8/100, y = 0.5, width = 0.8/100, height = 0.3,
            gp = gpar(fill = color_bar[i], col = NA))
}
for(i in 1:5) {
  grid.text(color_labels[i], x = 0.1 + (i-1)*0.8/4, y = 0.2, just = "center", gp = gpar(fontsize = 8))
}
grid.text("Age", x = 0.5, y = 0.85, just = "center", gp = gpar(fontsize = 10, fontface = "bold"))
popViewport()

## Gender track
# Create color vector for gender
temp_gender <- as.character(df_circular$gender_factor)  # Use gender_factor which has "Female"/"Male" labels
gender_colors_vec <- rep("grey", length(temp_gender))
gender_colors_vec[temp_gender == "Female"] <- sex_color["Female"]  # Keep using sex_color as defined
gender_colors_vec[temp_gender == "Male"] <- sex_color["Male"]

circos.track(
  factors = df_circular$factors,
  y = df_circular$y,
  ylim = c(0, 1),
  bg.border = "black",  # Add border to make it visible
  bg.col = "white",     # Add background color
  track.height = 0.15,  # Make it slightly taller
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
      border = "white",  # Add white border between sectors
      lwd = 0.5
    )
  }
)

# Debug: Check gender distribution
cat("\n============ DEBUG INFO ============\n")
cat("Gender distribution in data:\n")
print(table(df_circular$gender_factor, useNA = "always"))
cat("Gender colors being used:\n")
print(sex_color)  # Keep using sex_color as it's already defined
cat("First few gender values and colors:\n")
for(i in 1:min(5, nrow(df_circular))) {
  cat("Patient", i, "- Gender:", as.character(df_circular$gender_factor[i]), 
      "- Color:", gender_colors_vec[i], "\n")
}

draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360, 
            rou1 = 0.2, rou2 = 0, col = "white", border = "white")

# Add center labels
text(0, 0.05, nrow(df_circular), cex = 2, font = 2)
text(0, -0.05, "Diabetic PPV", cex = 1)
text(0, -0.15, "Participants", cex = 1)

# Add legend
legend_x <- 0.85
legend_y <- 0.6

# Gender legend  
text(legend_x, legend_y, "Gender", cex = 1.2, font = 2, adj = 0)

# Female
rect(legend_x, legend_y - 0.07, legend_x + 0.03, legend_y - 0.04, 
     col = sex_color["Female"], border = NA)
text(legend_x + 0.05, legend_y - 0.055, "Female", cex = 0.8, adj = 0)

# Male
rect(legend_x, legend_y - 0.13, legend_x + 0.03, legend_y - 0.10, 
     col = sex_color["Male"], border = NA)
text(legend_x + 0.05, legend_y - 0.115, "Male", cex = 0.8, adj = 0)

# Clear circos
circos.clear()

# Save circular plot
dev.copy(pdf, "circular_heatmap_diabetic_ppv.pdf", width = 10, height = 10)
dev.off()

cat("\n============ VISUALIZATION COMPLETE ============\n")
cat("Generated files:\n")
cat("1. age_analysis_ggstatsplot.pdf - Age comparison by gender\n")
cat("2. gender_distribution_donut.pdf - Gender distribution donut chart\n")
cat("3. age_distribution_detailed.pdf - Detailed age distribution\n")
cat("4. combined_age_gender_analysis.pdf - Combined visualization\n")
cat("5. circular_heatmap_diabetic_ppv.pdf - Circular heatmap with gradient age track\n")