# Code to replicate figure 2 from provided RDS data files
# Population dynamics and allele frequencies of the Myotis velifer aspatial example following disease introduction

# Set-up -----------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Set base directory where the input summary Rds files are stored
input_dir <- "Fig4_summary_data"

# Define where figure outputs should go
output_dir <- "figure_outputs"


# SIRD Pop Size (Fig 4a) ----------------------------------------------------

# Read in data
Fig_4a_data<- readRDS(paste0(input_dir, "/Fig_4a_data.Rds"))

# .....Plot SIRD
SIRD_fig <- ggplot(Fig_4a_data, aes(x = Year, y = mean, 
                                    color = Scenario, shape = Scenario, fill = Scenario)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = Fig_4a_data %>% filter(Year %% 10 == 0), 
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
# Read in data
Fig_4b_data<- readRDS(paste0(input_dir, "/Fig_4b_data.Rds"))

# Plot

Allele_Freq <- ggplot(Fig_4b_data, aes(x = Year, y = mean, color = Allele, fill = Allele)) +
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

