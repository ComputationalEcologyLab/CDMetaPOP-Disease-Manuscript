library(tidyverse)

# allele frequencies ------------------------------------------------------

# Set base directory
base_dir <- "Adaptive_Run_08/Adaptive_Run_08"
batches <- paste0("batch", 0:3)

# Function to calculate allele frequencies
calc_allele_freq <- function(file, year, run) {
  df <- read.csv(file)[, 21:28]
  freqs <- colSums(df) / (2 * nrow(df))
  tibble(
    Year = year,
    Run = run,
    Locus = rep(paste0("Locus ", 1:4), each = 2),
    Allele = rep(c("p", "q"), times = 4),
    Frequency = freqs
  )
}

# Collect summary allele frequencies for all batches
all_batches_summary <- map_dfr(batches, function(batch) {
  run_folders <- paste0("run0", batch, "mc", 0:9, "species0")
  
  all_freqs <- map_dfr(run_folders, function(run) {
    map_dfr(0:199, function(i) {
      file <- file.path(base_dir, run, paste0("ind", i, ".csv"))
      if (file.exists(file)) {
        calc_allele_freq(file, i, run)
      } else {
        warning(paste("Missing file:", file))
        return(NULL)
      }
    })
  })
  
  # Summarize across runs
  summary_df <- all_freqs %>%
    group_by(Year, Locus, Allele) %>%
    summarise(
      mean = mean(Frequency),
      sd = sd(Frequency),
      .groups = "drop"
    ) %>%
    mutate(
      Batch = paste0(batch),
      Allele = factor(Allele, levels = c("p", "q"))
    )
  
  return(summary_df)
})

## Clean data...

all_batches_summary <- all_batches_summary %>%
  rename(Scenario = Batch) %>%
  mutate(Scenario = recode(Scenario,
                           #"batch0" = "Null",
                           "batch0" = "Neutral",
                           "batch1" = "Resistance",
                           "batch2" = "Tolerance",
                           "batch3" = "Res. + Tol.")) %>%
  mutate(Scenario = factor(Scenario, levels = c(#"Null",
                                                "Neutral", "Resistance",
                                                "Tolerance", "Res. + Tol."))) %>%
  mutate(Locus = recode(Locus,
                        "Locus 1" = "Neutral 1",
                        "Locus 2" = "Neutral 2",
                        "Locus 3" = "Resistant",
                        "Locus 4" = "Tolerant"))

## Plot...

Allele_Freq <- ggplot(all_batches_summary, aes(x = Year, y = mean, color = Allele, fill = Allele)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  facet_grid(Locus ~ Scenario) +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0,200, by = 50),
                    name = "Year", sec.axis = sec_axis(~ ., name = "Scenario"))+
  scale_y_continuous(name = "Allele Frequency", sec.axis = sec_axis(~ ., name = "Locus")) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.x.bottom = element_text(size= 8),
        axis.text.y.left = element_text(size = 8),
        panel.border = element_rect(color = "black", 
                                    fill = NA,  linewidth = 0.3)) +
  labs(
    title = "b",
    y = "Allele Frequency",
    x = "Year",
    color = "Allele",
    fill = "Allele"
  )

Allele_Freq

ggsave(
  filename = "/mnt/DataDrive1/data/akwilliams/output_figures/Adaptive_Run_10/_allele_freq.png",
  plot = Allele_Freq,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


# population size ---------------------------------------------------------

batch_nums <- 0:3

# Collect all plot data
all_plot_data <- map_dfr(batch_nums, function(batch) {
  
  mc_paths <- paste0("/mnt/DataDrive1/data/akwilliams/outputs/Adaptive_Run_10/run0batch", 
                     batch, "mc", 0:24, "species0/summary_popAllTime_DiseaseStates.csv")
  
  all_reps <- lapply(mc_paths, function(path){
    df <- read_csv(path, show_col_types = FALSE)
    
    df %>%
      mutate(SID = sub("\\|.*", "", States_SecondUpdate)) %>%
      separate(SID, into = c("S", "I", "D"), sep = ";", convert = TRUE) #%>%
    #mutate(D = cumsum(D))  # cumulative deaths
  })
  
  combined_df <- bind_rows(all_reps, .id = "replicate")
  
  summary_df <- combined_df %>%
    mutate(N = S + I) %>%  # Total population size
    group_by(Year) %>%
    summarise(
      S_mean = mean(S), S_sd = sd(S),
      I_mean = mean(I), I_sd = sd(I),
      D_mean = mean(D), D_sd = sd(D),
      N_mean = mean(N), N_sd = sd(N),
      .groups = "drop"
    )
  
  plot_df <- summary_df %>%
    pivot_longer(cols = -Year, names_to = c("State", ".value"),
                 names_pattern = "(.)_(mean|sd)") %>%
    mutate(
      Batch = paste0("Batch ", batch),
      State = factor(State, levels = c("S", "I", "D", "N"))
    )
  
  return(plot_df)
})

## Cleaning data...

# remove N - not needed for final figure: 
all_plot_data <- all_plot_data %>%
  #filter(!State == "N") %>%
  rename(Scenario = Batch) %>%
  mutate(Scenario = recode(Scenario,
                           #Batch 0" = "Null",
                           "Batch 0" = "Neutral",
                           "Batch 1" = "Resistance",
                           "Batch 2" = "Tolerance",
                           "Batch 3" = "Res. + Tol.")) %>%
  mutate(Scenario = factor(Scenario, levels = c("Null", "Neutral", "Resistance",
                                                "Tolerance", "Res. + Tol.")))

## Plot...

SIDN_fig <- ggplot(all_plot_data, aes(x = Year, y = mean, 
                                     color = Scenario, shape = Scenario, fill = Scenario 
                                  )) +
  geom_line(linewidth = 0.5) +
  geom_point(data = all_plot_data %>% filter(Year %% 20 == 0), 
             size = 2, color = "black") +
  # REPLACE 0 WITH sd
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, color = NA) +
  facet_wrap(~ State, ncol = 2, scales = "free_y") +  # 4 panels: S, I, R, D
  scale_shape_manual(values = c(0,1,2,4,5)) +
  scale_x_continuous(limits = c(0,250)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA,  linewidth = 0.5)) +
  labs(
    #title = "a",
    x = "Year",
    y = "Population Size",
    color = "Scenario",
    fill = "Scenario",
    linetype = "Scenario"
  )
SIDN_fig
ggsave("/mnt/DataDrive1/data/akwilliams/output_figures/Adaptive_Run_10/_SIDN_by_state.png", SIDN_fig,
        width = 12, height = 7, dpi = 300, bg = "white")


# SIDP-specific figure ----------------------------------------------------

library(tidyverse)

batch_nums <- 0:3

# Collect all plot data
all_plot_data <- map_dfr(batch_nums, function(batch) {
  
  mc_paths <- paste0("Adaptive_Run_08/Adaptive_Run_08/run0batch", 
                     batch, "mc", 0:9, "species0/summary_popAllTime_DiseaseStates.csv")
  
  all_reps <- lapply(mc_paths, function(path){
    df <- read_csv(path, show_col_types = FALSE)
    
    df %>%
      mutate(SID = sub("\\|.*", "", States_SecondUpdate)) %>%
      separate(SID, into = c("S", "I", "D"), sep = ";", convert = TRUE) %>%
      mutate(
        # sum values in EnvReservoir_GetMetrics for P
        P = sapply(strsplit(EnvResivoir_GetMetrics, "\\|"), 
                   function(x) sum(as.numeric(x), na.rm = TRUE))
      )
  })
  
  combined_df <- bind_rows(all_reps, .id = "replicate")
  
  summary_df <- combined_df %>%
    mutate(N = S + I) %>%  # Total pop (optional now)
    group_by(Year) %>%
    summarise(
      S_mean = mean(S), S_sd = sd(S),
      I_mean = mean(I), I_sd = sd(I),
      D_mean = mean(D), D_sd = sd(D),
      P_mean = mean(P), P_sd = sd(P),
      .groups = "drop"
    )
  
  plot_df <- summary_df %>%
    pivot_longer(cols = -Year, names_to = c("State", ".value"),
                 names_pattern = "(.)_(mean|sd)") %>%
    mutate(
      Batch = paste0("Batch ", batch),
      State = factor(State, levels = c("S", "I", "D", "P"))
    )
  
  return(plot_df)
})

# Rename, recode, filter out N
all_plot_data <- all_plot_data %>%
  rename(Scenario = Batch) %>%
  mutate(Scenario = recode(Scenario,
                           "Batch 0" = "Neutral",
                           "Batch 1" = "Resistance",
                           "Batch 2" = "Tolerance",
                           "Batch 3" = "Res. + Tol.")) %>%
  mutate(Scenario = factor(Scenario, levels = c("Null", "Neutral", "Resistance",
                                                "Tolerance", "Res. + Tol.")))


# Final plot
SIDP_fig <- ggplot(all_plot_data, aes(x = Year, y = mean, 
                                      color = Scenario, shape = Scenario, fill = Scenario)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = all_plot_data %>% filter(Year %% 50 == 0), 
             size = 2, color = "black") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, color = NA) +
  facet_wrap(~ State, ncol = 2, scales = "free_y") +  # Now includes P
  scale_shape_manual(values = c(1,2,4,5)) +
  scale_color_manual(values = c("#a3a500", "#00bf7d", "#00b0f6", "#e76bf3")) +
  scale_fill_manual(values = c("#a3a500", "#00bf7d", "#00b0f6", "#e76bf3")) +
  #scale_x_continuous(limits = c(0,200)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  labs(
    x = "Year",
    y = "Total",
    color = "Scenario",
    fill = "Scenario",
    linetype = "Scenario",
    title = "a"
  )

SIDP_fig

ggsave("SIDP_by_state.png", SIDP_fig,
       width = 12, height = 7, dpi = 300, bg = "white")

# combine plots -----------------------------------------------------------

library(patchwork)

png(filename = "/mnt/DataDrive1/data/akwilliams/output_figures/Adaptive_Run_10/combined.png",
    width = 8, height = 9, units = "in", res = 300)
SIDP_fig / Allele_Freq 
dev.off()

# spatial spread v2 w looping ---------------------------------------------

library(tidyverse)
library(sf)
library(gganimate)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

# Load spatial patch points and background map
pts <- read.csv("figure_data/pts_reassigned.csv")
names(pts)[1] <- "patch"

countries <- ne_countries(scale = "medium", returnclass = "sf")
states <- ne_states(country = "United States of America", returnclass = "sf")
full <- rast("figure_data/FullENM_MedianTSS.tif")

pts_sf <- st_as_sf(pts, coords = c(2,3), crs = crs(full))
pts_sf <- st_transform(pts_sf, crs = crs(countries))

# Define base path and batches
base_path <- "/mnt/DataDrive1/data/akwilliams/outputs/Adaptive_Run_10"
batch_ids <- 0:3  # Update to match the number of batches you have
mc_ids <- 0:24     # Number of MC replicates per batch

for (batch in batch_ids) {
  message("Processing batch: ", batch)
  
  all_runs <- list()
  
  for (mc in mc_ids) {
    run_path <- file.path(
      base_path,
      paste0("run0batch", batch, "mc", mc, "species0/summary_popAllTime_DiseaseStates.csv")
    )
    
    if (!file.exists(run_path)) next
    
    disease_sum <- read.csv(run_path)
    
    # --- S, I, D ---
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
    
    # --- P ---
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
    filename = file.path("/mnt/DataDrive1/data/akwilliams/output_figures/Adaptive_Run_10", paste0("batch_", batch, "_disease.gif")),
    animation = 
      plot_obj,
    fps = 5
  )
  
}


# static spatial spread figure --------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

# Load spatial patch points and background map
pts <- read.csv("figure_data/pts_reassigned.csv")
names(pts)[1] <- "patch"

countries <- ne_countries(scale = "medium", returnclass = "sf")
states <- ne_states(country = "United States of America", returnclass = "sf")
full <- rast("figure_data/FullENM_MedianTSS.tif")

pts_sf <- st_as_sf(pts, coords = c(2, 3), crs = crs(full))
pts_sf <- st_transform(pts_sf, crs = crs(countries))

# Define base path and batches
base_path <- "/mnt/DataDrive1/data/akwilliams/outputs/Adaptive_Run_10"
batch_ids <- 0 # Only neutral
mc_ids <- 0:24
target_years <- c(0, 50, 200, 349)

all_batches <- list()

for (batch in batch_ids) {
  message("Processing batch: ", batch)
  all_runs <- list()
  
  for (mc in mc_ids) {
    run_path <- file.path(
      base_path,
      paste0("run0batch", batch, "mc", mc, "species0/summary_popAllTime_DiseaseStates.csv")
    )
    if (!file.exists(run_path)) next
    
    disease_sum <- read.csv(run_path)
    
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

ggsave(filename = "/mnt/DataDrive1/data/akwilliams/output_figures/Adaptive_Run_10/static_spread.png", 
       dpi = 300, width = 5, height = 4.5, units = "in", bg = "white")
