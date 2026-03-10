# Code to replicate figure 4 after running CDMetaPOP
# Population dynamics and allele frequencies of the Myotis velifer spatial example metapopulation following disease introduction
# Note that particularly the summarizing step for plotting allele frequencies 
#   can take a long time, so when possible, it helps to run on an HPC cluster or high performance machine
# If you run as a batch script that way, change the set-up info as needed and 
#   it will automatically write out the plots to the specified output folder

# Set-up -----------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Set base directory where the CDMetaPop raw output files are stored
# This should contain directories named run0batch0mc0species0, run0batch0mc1species0, etc
base_dir_spatial <- "Spatial_inputs"

# Define where figure outputs should go
output_dir <- "figure_outputs"

# For Fig 6 we ran CDMetaPop for 4 scenarios:
#   0 = Neutral
#   1 = Resistance
#   2 = Tolerance
#   3 = Resistance + Tolerance
batch_nums_spatial <- 0:3
batches_spatial <- paste0("batch", 0:3)

#Define the number of MCs used:
num_mcs_spatial <- 0:24

# SIDP Pop Size (Fig 6a) ----------------------------------------------------

# .....Collect all plot data----
message("Collecting plot data")
batch_vals <- batch_nums_spatial
mc_vals    <- num_mcs_spatial
base_dir   <- base_dir_spatial

all_disease_state_data <- map_dfr(batch_nums_spatial, function(batch) {
  
  mc_paths <- paste0(base_dir_spatial, "/run0batch", batch, "mc", num_mcs_spatial, "species0/summary_popAllTime_DiseaseStates.csv")
  
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

# .....Rename, recode, filter out N ----
all_disease_state_data <- all_disease_state_data %>%
  rename(Scenario = Batch) %>%
  mutate(Scenario = recode(Scenario,
                           "Batch 0" = "Neutral",
                           "Batch 1" = "Resistance",
                           "Batch 2" = "Tolerance",
                           "Batch 3" = "Res. + Tol.")) %>%
  mutate(Scenario = factor(Scenario, levels = c("Null", "Neutral", "Resistance",
                                                "Tolerance", "Res. + Tol.")))

# .....Plot SIDP ----
SIDP_fig <- ggplot(all_disease_state_data, aes(x = Year, y = mean, 
                                               color = Scenario, shape = Scenario, fill = Scenario)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = all_disease_state_data %>% filter(Year %% 50 == 0), 
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

ggsave(paste0(output_dir, "/SIDP_by_state.png"), SIDP_fig,
       width = 12, height = 7, dpi = 300, bg = "white")


# Allele frequencies (Fig 6b) ------------------------------------------------------

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

#### .....Collect summary allele frequencies for all batches ####
# This summarizes across scenarios and MCs
message("Collecting summary allele frequencies")
batch <- batches_spatial
mc_vals    <- num_mcs_spatial
base_dir   <- base_dir_spatial

all_batches_summary <- map_dfr(batch, function(batch) {
  mc_paths <- paste0("run0", batch, "mc", num_mcs_spatial, "species0")
  
  all_freqs <- map_dfr(mc_paths, function(run) {
    map_dfr(0:299, function(i) {
      file <- file.path(base_dir_spatial, run, paste0("ind", i, ".csv"))
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

#### .....Clean data ####
# This names the scenarios and the loci that are tracked

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

#### .....Plot ####

Allele_Freq <- ggplot(all_batches_summary, aes(x = Year, y = mean, color = Allele, fill = Allele)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  facet_grid(Locus ~ Scenario) +
  scale_x_continuous(limits = c(0, 300), breaks = seq(0,300, by = 50),
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
  filename = paste0(output_dir, "/spatial_allele_freq.png"),
  plot = Allele_Freq,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


# Combine Fig 6a and 6b  -----------------------------------------------------------

library(patchwork)

png(filename = paste0(output_dir, "/fig6_combined.png"),
    width = 8, height = 9, units = "in", res = 300)
SIDP_fig / Allele_Freq 
dev.off()

