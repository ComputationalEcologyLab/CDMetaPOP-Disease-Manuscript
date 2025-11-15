# Simulation Data and Code for CDMetaPOP-Disease Manuscript

This repository contains the input files, analysis scripts, summarized output data, and supplementary documentation required to replicate all analyses and figures in the manuscript: **Simulating the Landscape Eco-Evolution of Host-Pathogen Systems with CDMetaPOP**, a new module for CDMetaPOP for simulating disease transmission and spread in realistic landscapes (Manuscript submitted to Molecular Biology and Evolution, 2025)

**Authors:** Erin L. Landguth, Allison Williams, Marissa Roseman,  Miracle Amadi, Marcel Kouete, Rhys Farrer, Amy Haeseler, Orly Razgour, Chris Richardson, Byron Weckworth, Julie Weckworth, Flora Whiting-Fawcett, Duncan Wilson, Casey Day

**Permanent DOI for this repository:** TBD

**About This Repository:** This repository provides all necessary materials for the replication of our manuscript. It is organized by the figures presented in the paper.

**Software:** The simulations were run using CDMetaPOP v3.08. The main software repository (including the installer) is available at: https://github.com/ComputationalEcologyLab/CDMetaPOP

**Analysis:** Figures were generated using R or Python.

## Repository Structure

- **/Theoretical_Simulated_Expectations/:** These files reproduce Figure 1. Zipfiles for **/Input_Files/** and **/Output_Files/** for validating/verifying the module against classical SIR, SIRS, SEIR, and SIRP models. Also includes the Python script (run_verification.py) to generate the deterministic ODE solutions with the simulated expectations.

- **/Selection_Scenarios/:** These files reproduce Figure 2 for the aspatial eco-evolutionary demonstration (Null, Neutral, Resistance, Tolerance, and Res+Tol scenarios). Zipfile containing inputs, script for generating figure, and supplemental figure for the example of emergence of a disease response gennetic strategy and its impact on disease dynamics

- **/Example_Application/:** These files reproduce Figue 3 & 4 and contain all necessary input files to replicate the spatial example of novel disease emergence and spread through a bat system model. Contains: 
  - /Input_Files/ 
  - /Scripts/ to run CDMetaPOP and plot 
  - /Figures/ supplemental figures for the spatial example.

- **/Supplementary_Materials/:** Contains all materials formerly in the manuscript Appendix.
  - CDMetaPOP_Disease_Module_Documentation.docx: A detailed user guide (formerly Appendix Text A.1 and Table A.1).
  - CDMetaPOP_Life_Cycle_Diagram.png: (Formerly Figure A.1).
  - /videos/: Contains supplementary videos (formerly Video A.1-A.4) showing spatial simulation results.

## How to Replicate Manuscript Figures

**Option 1: Replicate Figures from Archived Data (Recommended)**

This is the fastest method and does not require running the CDMetaPOP simulations.

1. Clone this repository.

2. Open any R script provided in each folder.

3. The script will load the corresponding data file from /archived_outputs/.

4. Run the R script to generate the figure.

**Option 2: Full Replication by Re-running Simulations (Advanced)**

This method re-runs the CDMetaPOP simulations from scratch to generate the raw data.

1. Install CDMetaPOP: Follow the installation instructions at the main CDMetaPOP repository.

2. Navigate to a Figure Directory: cd figure_4_selection_demo/inputs/

3. Run Simulation: Run the CDMetaPOP executable using the input files in this directory (e.g., CDMetaPOP.exe PopVars_Neutral.csv).

4. Repeat: Repeat this process for all simulation scenarios (Neutral, Resistance, etc.)

## Supplementary Materials & Documentation

For detailed parameter definitions, user manual instructions, and supplementary videos, please see the /Supplementary_Materials/ folder in this repository.

