# Code to replicate figure 5 and supplemental gifs after running CDMetaPOP
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

# Set base directory where the CDMetaPop raw output files are stored
# This should contain directories named run0batch0mc0species0, run0batch0mc1species0, etc
base_dir_spatial <- "Spatial_inputs"

# Set directory where the provided points files are stored
base_dir_points <- "Points_files"

# Define where figure outputs should go
output_dir <- "Figure_outputs"

# Define batches
# For Fig 6 we ran CDMetaPop for 4 scenarios:
#   0 = Neutral
#   1 = Resistance
#   2 = Tolerance
#   3 = Resistance + Tolerance
batch_nums_spatial <- 0:3
batches_spatial <- paste0("batch", 0:3)

#Define the number of MCs used:
num_mcs_spatial <- 0:24

# Static spatial spread (Fig 5) --------------------------------------------

# Load spatial patch points and background map
pts <- read.csv(paste0(base_dir_points, "/points.csv"))
names(pts)[1] <- "patch"

countries <- ne_countries(scale = "medium", returnclass = "sf")
#states <- ne_states(country = "United States of America", returnclass = "sf")
states <- states50
full <- rast(paste0(base_dir_points, "/FullENM_MedianTSS.tif"))

pts_sf <- st_as_sf(pts, coords = c(2,3), crs = crs(full))
pts_sf <- st_transform(pts_sf, crs = crs(countries))


# Define batches for neutral-only image
batch_vals <- 0 # Only neutral
target_years <- c(0, 50, 200, 349)
base_dir   <- base_dir_spatial

all_batches <- list()

for (batch in batch_vals) {
  message("Processing static spread batch: ", batch)
  all_runs <- list()
    
    for (mc in num_mcs_spatial) {
      mc_paths <- file.path(
        base_dir_spatial,
        paste0("run0batch", batch, "mc", mc, "species0/summary_popAllTime_DiseaseStates.csv")
      )
    
    if (!file.exists(mc_paths)) next
    
    disease_sum <- read.csv(mc_paths)
    
    # Parse I values only
    state_vals <- strsplit(disease_sum$States_SecondUpdate, "\\|") %>%
      lapply(\(x) x[-1])
    
    split_I_vals <- lapply(state_vals, function(row) {
      do.call(rbind, lapply(row, function(x) {
        vals <- as.numeric(strsplit(x, ";")[[1]])
        vals[2]  # Only Infection (I) is index 2
      }))
    })
    
    i_df <- bind_rows(lapply(seq_along(split_I_vals), function(i) {
      data.frame(
        patch = 1:length(split_I_vals[[i]]),
        year = i - 1,  # Adjust for zero-indexing
        Value = split_I_vals[[i]]
      )
    }))
    
    i_df$mc <- mc
    all_runs[[length(all_runs) + 1]] <- i_df
  }
  
  # Combine all MCs for this batch
  batch_df <- bind_rows(all_runs) %>%
    filter(year %in% target_years) %>%
    group_by(patch, year) %>%
    summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(batch = batch)
  
  all_batches[[length(all_batches) + 1]] <- batch_df
}

# Combine all batches
final_df <- bind_rows(all_batches)

# Join with spatial data
plot_data <- right_join(pts_sf, final_df, by = "patch")

# Convert year to factor for faceting
plot_data$year <- factor(plot_data$year, levels = target_years)

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

ggsave(filename = paste0(output_dir, "/from_source/Figure_5.png"), 
       dpi = 300, width = 5, height = 4.5, units = "in", bg = "white")


# Spatial spread GIF (supplemental)---------------------------------------------

# Load spatial patch points and background map
pts <- read.csv(paste0(base_dir_points, "/points.csv"))
names(pts)[1] <- "patch"

countries <- ne_countries(scale = "medium", returnclass = "sf")
#states <- ne_states(country = "United States of America", returnclass = "sf")
states <- states50
full <- rast(paste0(base_dir_points, "/FullENM_MedianTSS.tif"))

pts_sf <- st_as_sf(pts, coords = c(2,3), crs = crs(full))
pts_sf <- st_transform(pts_sf, crs = crs(countries))

# .....Generate plot gifs---------------------
batch_vals <- batch_nums_spatial
base_dir   <- base_dir_spatial

for (batch in batch_nums_spatial) {
  message("Processing spread gif batch: ", batch)
  
  all_runs <- list()
  
  for (mc in num_mcs_spatial) {
    mc_paths <- file.path(
      base_dir_spatial,
      paste0("run0batch", batch, "mc", mc, "species0/summary_popAllTime_DiseaseStates.csv")
    )
    
    if (!file.exists(mc_paths)) next
    
    disease_sum <- read.csv(mc_paths)
    
    # ---- .....S, I, D ----
    state_vals <- strsplit(disease_sum$States_SecondUpdate, "\\|") %>% 
      lapply(\(x) x[-1])
    
    split_state_vals <- lapply(state_vals, function(row) {
      do.call(rbind, lapply(row, function(x) {
        vals <- as.numeric(strsplit(x, ";")[[1]])
        names(vals) <- c("S", "I", "D")
        vals
      }))
    })
    
    state_df <- bind_rows(lapply(seq_along(split_state_vals), function(i) {
      df <- as.data.frame(split_state_vals[[i]])
      df$patch <- 1:nrow(df)
      df$year <- i
      df
    }))
    
    state_long <- pivot_longer(state_df, cols = c("S", "I", "D"), names_to = "State", values_to = "Value")
    
    # --- .....P ----
    res_vals_all <- strsplit(disease_sum$EnvResivoir_GetMetrics, "\\|") %>%
      lapply(\(x) as.numeric(x[-1]))
    
    res_df <- bind_rows(lapply(seq_along(res_vals_all), function(i) {
      data.frame(
        patch = 1:length(res_vals_all[[i]]),
        year = i,
        State = "P",
        Value = res_vals_all[[i]]
      )
    }))
    
    run_df <- bind_rows(state_long, res_df)
    run_df$mc <- mc
    
    all_runs[[length(all_runs) + 1]] <- run_df
  }
  
  # Combine and average over MC replicates
  batch_df <- bind_rows(all_runs) %>%
    group_by(patch, year, State) %>%
    summarize(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Value_scaled = log1p(Value))
  
  # Join with spatial info
  combined_sf <- right_join(pts_sf, batch_df, by = "patch")
  combined_sf$State <- factor(combined_sf$State, levels = c("S", "I", "D", "P"))
  
  # Build plot
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
    filename = file.path(paste0(output_dir, "/from_source/batch_", batch, "_disease.gif")),
    animation = 
      plot_obj,
    fps = 5
  )
  
}
