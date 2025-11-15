# Simulation Data and Code for CDMetaPOP-Disease Manuscript

This repository contains the input files, analysis scripts, summarized output data, and supplementary documentation required to replicate all analyses and figures in the manuscript:

**Simulating the Landscape Eco-Evolution of Host-Pathogen Systems with CDMetaPOP**  A new module for CDMetaPOP for simulating disease transmission and spread in realistic landscapes (Manuscript submitted to Molecular Biology and Evolution, 2025)

**Authors:** Erin L. Landguth, Allison Williams, Marissa Roseman,  Miracle Amadi, Marcel Kouete, Rhys Farrer, Amy Haeseler, Orly Razgour, Chris Richardson, Byron Weckworth, Julie Weckworth, Flora Whiting-Fawcett, Duncan Wilson, Casey Day

**Permanent DOI for this repository:** TBD

**About This Repository:** This repository provides all necessary materials for the replication of our manuscript. It is organized by the figures presented in the paper.

**Software:** The simulations were run using CDMetaPOP v3.08. The main software repository (including the installer) is available at: https://github.com/ComputationalEcologyLab/CDMetaPOP

**Analysis:** Figures were generated using R or Python.

**Repository Structure**

- /Theoretical_Simulated_Expectations/: These files reproduce Figure 1. Zipfiles for /Input_Files/ and /Output_Files/ for validating/verifying the module against classical SIR, SIRS, SEIR, and SIRP models. Also includes the Python script (run_verification.py) to generate the deterministic ODE solutions with the simulated expectations.

- /Selection_Scenarios/: These files reproduce Figure 2 for the aspatial eco-evolutionary demonstration (Null, Neutral, Resistance, Tolerance, and Res+Tol scenarios). Contains input files (/Inputs/) for the aspatial eco-evolutionary demonstration (Null, Neutral, Resistance, Tolerance, and Res+Tol scenarios). Zipfile containing inputs, script for generating figure, and supplemental figure for the example of emergence of a disease response gennetic strategy and its impact on disease dynamics

- /Example_Application/: These files reproduce Figue 3 & 4 and contain all necessary input files to replicate the spatial example of novel disease emergence and spread through a bat system model. Contains /Input_Files/, /Scripts/ to run CDMetaPOP and plot /Figures/ for the spatial example.

- /Supplementary_Materials/: Contains all materials formerly in the manuscript Appendix.

- - CDMetaPOP_Disease_Module_Documentation.md: A detailed user guide (formerly Appendix Text A.1 and Table A.1).

- - CDMetaPOP_Life_Cycle_Diagram.png: (Formerly Figure A.1).

- - /videos/: Contains supplementary videos (formerly Video A.1-A.4) showing spatial simulation results.

How to Replicate Manuscript Figures

Option 1: Replicate Figures from Archived Data (Recommended)

This is the fastest method and does not require running the CDMetaPOP simulations.

Clone this repository.

Open any R script in /analysis_scripts/ (e.g., plot_figure_4.R).

The script will load the corresponding data file from /archived_outputs/.

Run the R script to generate the figure.

Option 2: Full Replication by Re-running Simulations (Advanced)

This method re-runs the CDMetaPOP simulations from scratch to generate the raw data.

Install CDMetaPOP: Follow the installation instructions at the main CDMetaPOP repository.

Navigate to a Figure Directory: cd figure_4_selection_demo/inputs/

Run Simulation: Run the CDMetaPOP executable using the input files in this directory (e.g., CDMetaPOP.exe PopVars_Neutral.csv).

Repeat: Repeat this process for all simulation scenarios (Neutral, Resistance, etc.)

Analyze Outputs: Once simulations are complete, use the scripts in /analysis_scripts/ to process the new raw output files and generate the figures.

Supplementary Materials & Documentation

For detailed parameter definitions, user manual instructions, and supplementary videos, please see the /supplementary_materials/ folder in this repository.

## Theoretical_Validation
####   Input_Files:
- Zipfiles containing all necessary input files to replicate the theoretical validation of disease modules
####   Output_files:
- Zipfiles contining simulation output from theoretical validation of disease modules

## Spatial_Example
####   Input_Files:
- Zipfile containing all necessary input files to replicate the spatial example of novel disease emergence and spread through a bat system model
####   Scripts:
- Scripts to run CDMetaPop and plot figures for the spatial example
####   Figures:
- Supplemental figures for the spatial example

## Disease_Response_Example
- Zipfile containing inputs, script for generating figure, and supplemental figure for the example of emergence of a disease response gennetic strategy and its impact on disease dynamics 
