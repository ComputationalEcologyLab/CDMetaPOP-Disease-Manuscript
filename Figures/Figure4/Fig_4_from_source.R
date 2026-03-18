# Code to replicate figure 2 after running CDMetaPOP
# Population dynamics and allele frequencies of the Myotis velifer aspatial example following disease introduction
# Note that particularly the summarizing step for plotting allele frequencies 
#   can take a long time, so when possible, it helps to run on an HPC cluster or high performance machine
# If you run as a batch script that way, change the set-up info as needed and 
#   it will automatically write out the plots to the specified output folder

# Set-up -----------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Set base directory where the CDMetaPop raw output files are stored
# This should contain directories named run0batch0mc0species0, run0batch0mc1species0, etc
base_dir_aspatial <- "Fig4_from _source_data/Aspatial_inputs"

# Define where figure outputs should go
output_dir <- "figure_outputs"

# Define batches for each fig
# For Fig 4 we ran CDMetaPop for 5 scenarios:
# Numbering of these scenarios starts at 0 since CDMetaPop using python numbering syntax
#   0 = Null
#   1 = Neutral
#   2 = Resistance
#   3 = Tolerance
#   4 = Resistance + Tolerance
batch_nums_aspatial <- 0:4
batches_aspatial <- paste0("batch", 0:4)

#Define the number of MCs used:
num_mcs_aspatial <- 0:99

# SIRD Pop Size (Fig 4a) ----------------------------------------------------

# .....Collect all plot data----
message("Collecting plot data")

batch_vals <- batch_nums_aspatial
mc_vals    <- num_mcs_aspatial
base_dir   <- base_dir_aspatial

all_disease_state_data <- map_dfr(batch_vals, function(batch) {
  
  mc_paths <- paste0(base_dir_aspatial, "/run0batch", batch, "mc", num_mcs_aspatial, "species0/summary_popAllTime_DiseaseStates.csv")
  
  all_reps <- lapply(mc_paths, function(path){
    df <- read_csv(path, show_col_types = FALSE)
    
    df %>%
      mutate(SIRD = sub("\\|.*", "", States_SecondUpdate)) %>%
      separate(SIRD, into = c("S", "I", "R", "D"), sep = ";", convert = TRUE)
  })
  
  combined_df <- bind_rows(all_reps, .id = "replicate")
  
  summary_df <- combined_df %>%
    mutate(N = S + I + R) %>%  # Total population size
    group_by(Year) %>%
    summarise(
      S_mean = mean(S), S_sd = sd(S),
      I_mean = mean(I), I_sd = sd(I),
      R_mean = mean(R), R_sd = sd(R),
      D_mean = mean(D), D_sd = sd(D),
      N_mean = mean(N), N_sd = sd(N),
      .groups = "drop"
    )
  
  plot_df <- summary_df %>%
    pivot_longer(cols = -Year, names_to = c("State", ".value"),
                 names_pattern = "(.)_(mean|sd)") %>%
    mutate(
      Batch = paste0("Batch ", batch),
      State = factor(State, levels = c("S", "I", "R", "D", "N"))
    )
  
  return(plot_df)
})

# .....Rename, recode, filter out N ----
all_disease_state_data <- all_disease_state_data %>%
  filter(!State == "N") %>%
  rename(Scenario = Batch) %>%
  mutate(Scenario = recode(Scenario,
                           "Batch 0" = "Null",
                           "Batch 1" = "Neutral",
                           "Batch 2" = "Resistance",
                           "Batch 3" = "Tolerance",
                           "Batch 4" = "Res. + Tol.")) %>%
  mutate(Scenario = factor(Scenario, levels = c("Null", "Neutral", "Resistance",
                                                "Tolerance", "Res. + Tol.")))

# .....Plot SIRD ----
SIRD_fig <- ggplot(all_disease_state_data, aes(x = Year, y = mean, 
                                               color = Scenario, shape = Scenario, fill = Scenario)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = all_disease_state_data %>% filter(Year %% 10 == 0), 
             size = 2, color = "black") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15, color = NA) +
  facet_wrap(~ State, ncol = 2, scales = "free_y") +  # 4 panels: S, I, R, D
  scale_shape_manual(values = c(0,1,2,4,5)) +
  #scale_x_continuous(limits = c(0,200)) +
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

SIRD_fig

ggsave(paste0(output_dir, "/SIRD_by_state.png"), SIRD_fig,
       width = 12, height = 7, dpi = 300, bg = "white")


# Allele frequencies (Fig 4b) ------------------------------------------------------

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
batch <- batches_aspatial
mc_vals    <- num_mcs_aspatial
base_dir   <- base_dir_aspatial

all_batches_summary <- map_dfr(batch, function(batch) {
  mc_paths <- paste0("run0", batch, "mc", num_mcs_aspatial, "species0")
  
  all_freqs <- map_dfr(mc_paths, function(run) {
    map_dfr(0:99, function(i) {
      file <- file.path(base_dir_aspatial, run, paste0("ind", i, ".csv"))
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
                           "batch0" = "Null",
                           "batch1" = "Neutral",
                           "batch2" = "Resistance",
                           "batch3" = "Tolerance",
                           "batch4" = "Res. + Tol.")) %>%
  mutate(Scenario = factor(Scenario, levels = c("Null", "Neutral", "Resistance",
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
  scale_x_continuous(limits = c(0, 100), breaks = seq(0,100, by = 50),
                     name = "Year", sec.axis = sec_axis(~ ., name = "Scenario"))+
  scale_y_continuous(name = "Allele Frequency", sec.axis = sec_axis(~ ., name = "Locus")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
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
  filename = paste0(output_dir, "/aspatial_allele_freq.png"),
  plot = Allele_Freq,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


# Combine Fig 4a and 4b  -----------------------------------------------------------

library(patchwork)

png(filename = paste0(output_dir, "/fig4_combined.png"),
    width = 8, height = 9, units = "in", res = 300)
SIRD_fig / Allele_Freq 
dev.off()

