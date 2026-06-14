# Sebastian Nöth
# 09 May 2026

######################################################################################################################################################################################
#### TABLE ECO-GROUP ####
# Define function for calculating the standard error and pasting it together with the value
mean_se <- function(x) {
  
  mean_x <- mean(x, na.rm = TRUE)
  
  se_x <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  
  paste0(
    sprintf("%.3f", mean_x),
    " ± ",
    sprintf("%.3f", se_x)
  )
}


# Compile necessary data into a data frame
summary_table_eco <- data.frame(
  
  Metric = c(
    "Three-timer sampling completeness",
    "Extinction rate PC (overall)",
    "Extinction rate PC (spinose)",
    "Extinction rate PC (non-spinose)"
  ),
  
  Mixed_symbiotic = c(
    mean_se(divdyn_open_ocean_mixed_symb$samp3t),
    mean_se(divdyn_open_ocean_mixed_symb$extPC),
    mean_se(divdyn_open_ocean_mixed_symb_spinose$extPC),
    mean_se(divdyn_open_ocean_mixed_symb_non$extPC)
  ),
  
  Mixed_non_symbiotic = c(
    mean_se(divdyn_open_ocean_mixed_nosymb$samp3t),
    mean_se(divdyn_open_ocean_mixed_nosymb$extPC),
    mean_se(divdyn_open_ocean_mixed_nosymb_spinose$extPC),
    mean_se(divdyn_open_ocean_mixed_nosymb_non$extPC)
  ),
  
  Thermocline = c(
    mean_se(divdyn_open_ocean_thermocline$samp3t),
    mean_se(divdyn_open_ocean_thermocline$extPC),
    mean_se(divdyn_open_ocean_thermocline_spinose$extPC),
    mean_se(divdyn_open_ocean_thermocline_non$extPC)
  ),
  
  Subthermocline = c(
    mean_se(divdyn_open_ocean_subthermocline$samp3t),
    mean_se(divdyn_open_ocean_subthermocline$extPC),
    mean_se(divdyn_open_ocean_subthermocline_spinose$extPC),
    mean_se(divdyn_open_ocean_subthermocline_non$extPC)
  ),
  
  High_latitude = c(
    mean_se(divdyn_Highlatitude$samp3t),
    mean_se(divdyn_Highlatitude$extPC),
    "NA",
    "NA"
  ),
  
  Upwelling_High_productivity = c(
    mean_se(divdyn_Upwelling$samp3t),
    mean_se(divdyn_Upwelling$extPC),
    mean_se(divdyn_Upwelling_spinose$extPC),
    "NA"
  )
)


# install and load necessary packages
#install.packages("gt")
#install.packages("webshot2")

library(gt)
library(webshot2)

# plot/create table with desired style 
summary_table_eco_gt <- summary_table_eco |> 
  gt() |>
  
  tab_options(
    table.width = px(800),
    data_row.padding = px(6),
    summary_row.padding = px(6),
    table.font.size = px(20)
  ) |>
  
  cols_label(
    Metric = "Metric",
    Mixed_symbiotic = md("Mixed <br> symbiotic"),
    Mixed_non_symbiotic = md("Mixed <br> non-symbiotic"),
    Thermocline = "Thermocline",
    Subthermocline = md("Sub- <br> thermocline"),
    High_latitude = "High latitude",
    Upwelling_High_productivity = "Upwelling"
  ) |>
  
  cols_align(align = "center") |>
  
  cols_width(everything() ~ px(150)) |>
  
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "black",
      weight = px(1)
    ),
    locations = cells_body()
  ) |>
  
  tab_style(
    style = list(
      cell_fill(color = "grey"),
      cell_borders(
        sides = "all",
        color = "black",
        weight = px(1)
      ),
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels()
  ) |>
  
  tab_style(
    style = cell_fill(color = "lightgrey"),
    locations = cells_body(
      columns = Metric
    )
  )


# download the table
gtsave(
  summary_table_eco_gt,
  "summary_table_eco.png",
  vwidth = 1200,
  vheight = 600
)


######################################################################################################################################################################################
#### TABLE MORPHO-GROUP ####
# Compile necessary data into a data frame
summary_table_morpho <- data.frame(
  
  Metric = c(
    "Three-timer sampling completeness",
    "Extinction rate (per capita)",
    "Extinction rate (3t)",
    "Origination rate (per capita)",
    "Origination rate (3t)",
    "Turnover (per capita)",
    "Turnover (3t)"
  ),
  
  Spinose = c(
    mean_se(divdyn_spinose$samp3t),
    mean_se(divdyn_spinose$extPC),
    mean_se(divdyn_spinose$ext3t),
    mean_se(divdyn_spinose$oriPC),
    mean_se(divdyn_spinose$ori3t),
    mean_se(divdyn_spinose$extPC + divdyn_spinose$oriPC),
    mean_se(divdyn_spinose$ext3t + divdyn_spinose$ori3t)
  ),
  
  Non_spinose = c(
    mean_se(divdyn_non$samp3t),
    mean_se(divdyn_non$extPC),
    mean_se(divdyn_non$ext3t),
    mean_se(divdyn_non$oriPC),
    mean_se(divdyn_non$ori3t),
    mean_se(divdyn_non$extPC + divdyn_non$oriPC),
    mean_se(divdyn_non$ext3t + divdyn_non$ori3t)
  )
)


# plot/create table with desired style 
summary_table_morpho_gt <- summary_table_morpho |>
  gt() |>
  
  tab_options(
    table.width = px(600),
    data_row.padding = px(6),
    table.font.size = px(12)
  ) |>
  
  cols_label(
    Metric = "Metric",
    Spinose = "Spinose",
    Non_spinose = "Non-spinose"
  ) |>
  
  cols_align(
    align = "center"
  ) |>
  
  cols_width(
    everything() ~ px(180)
  ) |>
  
  tab_style(
    style = cell_borders(
      sides = "all",
      color = "black",
      weight = px(1)
    ),
    locations = cells_body()
  ) |>
  
  tab_style(
    style = list(
      cell_fill(color = "grey"),
      cell_text(weight = "bold"),
      cell_borders(
        sides = "all",
        color = "black",
        weight = px(1)
      )
    ),
    locations = cells_column_labels()
  ) |>
  
  tab_style(
    style = list(
      cell_fill(color = "lightgrey")
    ),
    locations = cells_body(
      columns = Metric
    )
  )

# download the table
gtsave(
  summary_table_morpho_gt,
  "summary_table_morpho.png",
  vwidth = 900,
  vheight = 500,
  expand = 10
)


######################################################################################################################################################################################
#### HISTOGRAM SPIONSE TO NON-SPINOSE PROPORTIONS ####
# download necessary packages
library(ggplot2)
library(tidyr)

# convert data into long form
sum_foraminifer_long <- pivot_longer(
  sum_foraminifer_df,
  cols = c("Spinose_Percent", "Nonspinose_Percent"),
  names_to = "Type",
  values_to = "Percent"
)


# set up the download
png("Morpho_proportions.png", width = 2640, height = 1848, res = 300)


# plot the data with desired parameters/style
ggplot(sum_foraminifer_long,
       aes(x = Ecogroup,
           y = Percent,
           fill = Type)) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_discrete(
    name = "Morpho-group:",
    labels = c("Non-Spinose", "Spinose")
  )+
  
  geom_text(
    aes(
      y = ifelse(Type == "Spinose_Percent", 25, 75),
      label = paste0(round(Percent, 1), "%"),
    )
  )+
  
  labs(
    x = "Eco-group",
    y = "Percent",
    
  )+
  
  theme_minimal() +
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  
)

# close of the download
dev.off()


######################################################################################################################################################################################
#### LINE DIAGRAM OF SPIONSE AND NON-SPINOSE EXTINCTION RATES THROUGHOUT THE CENOZOIC ####
# set up the download
png("Morpho_timeseries.png", width = 2640, height = 1848, res = 300)

# load in the stage data
data(stages)

# plot the diagram base/background and axis
tsplot(stages, boxes=c("short","system"), shading="short",
       xlim=82:95,
       ylim=c(0,1.4),
       boxes.col=c("col","systemCol"), 
       labels.args=list(cex=1))

# add the line for spinose
lines(divdyn_spinose$age_bins,
      divdyn_spinose$extPC,
      col = "#00BFC4",
      lwd = 4)

# add the line for non-spinose
lines(divdyn_non$age_bins,
      divdyn_non$extPC,
      col = "#F8766D",
      lwd = 4)

# add the legend
legend("topright",
       legend = c("Spinose", "Non-spinose"),
       col = c("#00BFC4", "#F8766D"),
       lwd = 4,
       bty = "o",
       bg = "white")

# close of the download
dev.off()


######################################################################################################################################################################################
####BAR CHART FOR THE EXTINCTION RATES (WITH STANDARD ERROR) FOR EG 1-4 FOR BOTH NON-SPINOSE AND SPINOSE ####
# load necessary packages
library(dplyr)

# summarize values into a data frame
extinction_summary <- data.frame(
  ecogroup = c(
    "subthermocline", "subthermocline",
    "thermocline", "thermocline",
    "mixed_nosymb", "mixed_nosymb",
    "mixed_symb", "mixed_symb"
  ),
  
  spinose = c(
    "spinose", "non-spinose",
    "spinose", "non-spinose",
    "spinose", "non-spinose",
    "spinose", "non-spinose"
  ),
  
  values = c(
    mean(divdyn_open_ocean_subthermocline_spinose$extPC, na.rm = T),
    mean(divdyn_open_ocean_subthermocline_non$extPC, na.rm = T),
    
    mean(divdyn_open_ocean_thermocline_spinose$extPC, na.rm = T),
    mean(divdyn_open_ocean_thermocline_non$extPC, na.rm = T),
    
    mean(divdyn_open_ocean_mixed_nosymb_spinose$extPC, na.rm = T),
    mean(divdyn_open_ocean_mixed_nosymb_non$extPC, na.rm = T),
    
    mean(divdyn_open_ocean_mixed_symb_spinose$extPC, na.rm = T),
    mean(divdyn_open_ocean_mixed_symb_non$extPC, na.rm = T)
  ),
  error = c(
    sd(divdyn_open_ocean_subthermocline_spinose$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    sd(divdyn_open_ocean_subthermocline_non$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    
    sd(divdyn_open_ocean_thermocline_spinose$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    sd(divdyn_open_ocean_thermocline_non$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    
    sd(divdyn_open_ocean_mixed_nosymb_spinose$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    sd(divdyn_open_ocean_mixed_nosymb_non$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    
    sd(divdyn_open_ocean_mixed_symb_spinose$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC))),
    sd(divdyn_open_ocean_mixed_symb_non$extPC, na.rm = TRUE) / sqrt(sum(!is.na(divdyn_open_ocean_subthermocline_spinose$extPC)))
  )
)

summary_df <- extinction_summary %>%
  group_by(ecogroup, spinose) %>%
  summarise(
    mean_ext = values,
    se_ext = error,
    .groups = "drop"
  )

# plot the data with desired parameters/style
png("Extiction_values.png", width = 2500, height = 2100, res = 300)

# plot data with desired 
ggplot(summary_df, aes(x = ecogroup, y = mean_ext, fill = spinose)) +
  
  geom_col(position = position_dodge(0.7)) +
  
  geom_errorbar(
    aes(ymin = mean_ext - se_ext,
        ymax = mean_ext + se_ext),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(0.85, 0.90)
  )+
    
  labs(
    x = "Eco-group",
    y = "Mean extinction rate (extPC)",
    fill = "Morpho-group:"
  ) +
  
  scale_x_discrete(labels = c(
    "subthermocline" = "Sub-thermocline",
    "thermocline" = "Thermocline",
    "mixed_nosymb" = "Mixed \n (no symbionts)",
    "mixed_symb" = "Mixed \n (symbiotic)"
  )
  )

# close off the download
dev.off()
