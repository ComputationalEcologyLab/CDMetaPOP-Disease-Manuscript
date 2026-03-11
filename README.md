# CDMetaPOP-Disease Manuscript Figure Reproduction

### Table of Contents  
- [Introduction](#introduction)
- [Software Install](#software-install)  
- [Figure 1](#figure-1)
- [Figure 2](#figure-2)
- [Figure 3](#figure-3)
- [Figure 4](#figure-4)
- [Figure 5](#figure-5)
- [Figure 6](#figure-6)
<a name="headers"/>

## Introduction

This repository contains the files to recreate and reproduce the figures and anaylses in the manuscript:
> [!NOTE]
> Landguth, E. L., et al. (Submitted). "Simulating the Landscape Eco-Evolution of Host-Pathogen Systems with CDMetaPOP"

>**Full Authorship:** Erin L. Landguth, Allison Williams, Marissa Roseman, Miracle Amadi, Baylor Fain, Marcel Kouete, Rhys A. Farrer, Amy J. Haeseler, Orly Razgour, Chris Richardson, Byron Weckworth, Julie Weckworth, Flora Whiting-Fawcett, Duncan Wilson, Casey C. Day

## Software Install

## Figure 1

## Figure 2

## Figure 3

## Figure 4

## Figure 5

## Figure 6



## Repository Structure

- **/Theoretical_Validation/:** These files reproduce Figure 1. Zipfiles for **/Input_Files/** and **/Output_Files/** for validating/verifying the module against classical SIR, SIRS, SEIR, and SIRP models. Also includes the Python script (run_verification.py) to generate the deterministic ODE solutions with the simulated expectations.

- **/Selection_Scenarios/:** These files reproduce Figure 2 for the aspatial eco-evolutionary demonstration (Null, Neutral, Resistance, Tolerance, and Res+Tol scenarios). The folder includes a zip file containing the input files, the script used to generate the figure, and a supplemental figure illustrating the emergence of a disease-response genetic strategy and its impact on disease dynamics.

- **/Example_Application/:** These files reproduce Figue 3 & 4 and contain all necessary input files to replicate the spatial example of novel disease emergence and spread through a bat system model. This folder contains:
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

2. Navigate to a Figure Directory: cd /Theoretical_Validation/Inputs_Files/

3. Run Simulation: Run the CDMetaPOP program using the input files in this directory

4. Repeat: Repeat this process for all simulation scenarios (Neutral, Resistance, etc.)

## Supplementary Materials & Documentation

For detailed parameter definitions, user manual instructions, and supplementary videos, please see the /Supplementary_Materials/ folder in this repository.

