# Code to replicate figure 3 and supplemental gifs after running CDMetaPOP
# Disease spread of the Myotis velifer spatial example metapopulation following disease introduction

# Set-up -----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
#For mapping and animations:
library(sf) 
library(gganimate)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

# Set base directory where the input summary Rds and points files are stored
input_dir <- "Fig5_summary_data"

# Define where figure outputs should go
output_dir <- "figure_outputs"

# Define batches
# For Fig 6 we ran CDMetaPop for 4 scenarios:
#   0 = Neutral
#   1 = Resistance
#   2 = Tolerance
#   3 = Resistance + Tolerance
batch_nums_spatial <- 0:3

# Static spatial spread (Fig 5) --------------------------------------------

# Read in data
# plot_data<- readRDS(paste0(input_dir, "/Fig_5_data.Rds"))

#Get map outlines
countries <- ne_countries(scale = "medium", returnclass = "sf")
states <- states50

#Define the patches where the disease was introduced
highlight_patches <- c(4,60,78)

highlight_sf <- plot_data %>%
  filter(patch %in% highlight_patches, year == 0) %>%
  mutate(Season = case_when(Season == "winter" ~ "Disease"))

# Plot
ggplot(plot_data) +
  geom_sf(data = countries, fill = NA, color = "gray30", size = 0.3) +
  geom_sf(data = states, fill = NA, color = "gray50", size = 0.2) +
  geom_sf(aes(color = Value, size = Value), alpha = 0.8) +
  geom_sf(data = highlight_sf, aes(shape = Season), fill = "red", color = "black", size = 3) +
  scale_shape_manual(name = "", values = 23) +
  scale_color_viridis_c(option = "magma") +
  scale_size(range = c(1, 6)) +
  facet_wrap(~year, ncol = 2, 
             labeller = labeller(year = c("0" = "Year 0", "50" = "Year 50", 
                                          "200" = "Year 200", "349" = "Year 350"))) +
  coord_sf(xlim = c(-120, -95), ylim = c(15, 40)) +
  scale_y_continuous(breaks = seq(15, 40, by = 10)) +
  scale_x_continuous(breaks = seq(-120, -95, by = 10)) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Infection",
    size = "Infection"
  ) +
  guides(color = guide_legend(title = "No. infected"),
         size = guide_legend(title = "No. infected")) +
  theme_minimal(base_size = 12)

ggsave(filename = paste0(output_dir, "/Figure_5.png"), 
       dpi = 300, width = 5, height = 4.5, units = "in", bg = "white")


# Spatial spread GIF (supplemental)---------------------------------------------

# Read in data
Fig_5_supplemental_data<- readRDS(paste0(input_dir, "/Fig_5_supplemental_data.Rds"))

# Get map outlines
countries <- ne_countries(scale = "medium", returnclass = "sf")
#states <- ne_states(country = "United States of America", returnclass = "sf")
states <- states50

# Generate plot gifs
for (batch in batch_nums_spatial) {
  combined_sf<- Fig_5_supplemental_data[[as.character(batch)]]

  # Build plots
  plot_obj <- ggplot(combined_sf) +
    geom_sf(data = countries, fill = NA, color = "gray30", size = 0.3) +
    geom_sf(data = states, fill = NA, color = "gray50", size = 0.2) +
    geom_sf(aes(color = Value_scaled, size = Value_scaled), alpha = 0.8) +
    scale_color_viridis_c(option = "magma") +
    scale_size(range = c(1, 6)) +
    facet_wrap(~State) +
    coord_sf(xlim = c(-120, -95), ylim = c(15, 40)) +
    scale_y_continuous(breaks = seq(15, 40, by = 10)) +
    scale_x_continuous(breaks = seq(-120, -95, by = 10)) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Latitude",
      y = "Longitude",
      color = "Relative value",
      size = "Relative value",
      title = paste("Batch", batch, "- Year: {frame_time}")
    ) +
    transition_time(year) +
    ease_aes("linear") +
    guides(color = guide_legend(), size = guide_legend())
  
  # Animate and save
  anim_save(
    filename = file.path(paste0(output_dir, "/batch_", batch, "_disease.gif")),
    animation = 
      plot_obj,
    fps = 5
  )
}
