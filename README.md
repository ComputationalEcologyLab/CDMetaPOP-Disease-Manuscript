# CDMetaPOP-Disease Manuscript Figure Reproduction

![Compatibility Check](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/actions/workflows/python-check.yml/badge.svg)


### Table of Contents  
[Introduction](#introduction)  
[Software Install](#software-install)  
[Recreation](#recreation)  
[Reproduction](#reproduction)  
<a name="headers"/>

## Introduction

This repository contains the files to recreate and reproduce the figures and analyses in the manuscript:  
**Landguth, E. L., et al. (Submitted). "Simulating the Landscape Eco-Evolution of Host-Pathogen Systems with CDMetaPOP"**

>**Full Authorship:** Erin L. Landguth, Allison Williams, Marissa Roseman, Miracle Amadi, Baylor Fain, Marcel Kouete, Rhys A. Farrer, Amy J. Haeseler, Orly Razgour, Chris Richardson, Byron Weckworth, Julie Weckworth, Flora Whiting-Fawcett, Duncan Wilson, Casey C. Day

## Software Install

To recreate the figures summary data is provided and all that needs to be done is to run the provided code. 
Figure 1, 2, and 3 use python and Figure 4, 5, and 6 use R.

>[!NOTE]
>The provided codes should work with current releases of python and R packages, but environment setup is also provided.  

### To begin
>[!Important]
> Clone this repo and use it as the working directory

#### Python environment
To setup the python environment make sure either [Anaconda](https://www.anaconda.com/docs/getting-started/anaconda/install), [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install), or your favorite conda provider is installed. 
Then enter the following instructions in the terminal (please choose your favorite terminal: Anaconda prompt, Spyder terminal, bash, etc.):
```
conda create --name CDMetaPOP_Disease_env python
conda activate CDMetaPOP_Disease_env
pip install -r requirements.txt
```

#### R environment
To setup the R environment use [Rstudio](https://posit.co/download/rstudio-desktop/) and make sure `Renv` is installed (`install.packages("renv"`)). 
Then enter the following instructions in the terminal:
```
renv::restore()
```


## Recreation

>[!Important]
> Follow the instructions in [Software Install](#software-install) and begin each figure making process with this repo as your working directory

### Python

#### To recreate Figure 1-3:

- Navigate to the desired directory within [Figures](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures), either [Figure1](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure1), [Figure2](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure2), or [Figure3](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure3).
    - ``` cd Figures/FigureX``` 
- run
    - ``` python Figure_X_from_data.py ```

The image for figure X (Figure_X.png) will be created in the *figure_outputs* directory within *FigureX*.

### R

#### To recreate Figure 4-6:

- Double click the [.Rprofile](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/blob/manuscript_prep/.Rprofile) file, which will open Rstudio
- Navigate to the desired directory within [Figures](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures), either [Figure4](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure4), [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure5), or [Figure6](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure6).
    - ``` setwd("Figures/FigureX") ``` 
- run in the console
    - ``` source("Figure_X_from_data.R") ```

The image for figure X (Figure_X.png) will be created in the *figure_outputs* directory with in *FigureX*.

## Reproduction

>[!Important]
> Follow the instructions in [Software Install](#software-install) and begin each figure making process with this repo as your working directory

>[!Note]
> CDMetaPOP v3.08 is provided in this repo, so CDMetaPOP does not have to be install again unless desired.

### Figure 1

>[!Note]
> ~50 MB will be needed to store the data for this figure

To reproduce the data for figure 1, four cases must be simulated. 
The code to run each simulation is:
- SIR
    - ``` python CDMetaPOP/CDMetaPOP.py Figures/Figure1/CDMetaPOP_inputs/SIR/ RunVars.csv output/output ```
- SEIR
    - ``` python CDMetaPOP/CDMetaPOP.py Figures/Figure1/CDMetaPOP_inputs/SEIR/ RunVars.csv output/output ```
- SIRP
    - ``` python CDMetaPOP/CDMetaPOP.py Figures/Figure1/CDMetaPOP_inputs/SIRP/ RunVars.csv output/output ```
- SIRS
    - ``` python CDMetaPOP/CDMetaPOP.py Figures/Figure1/CDMetaPOP_inputs/SIRS/ RunVars.csv output/output ```

Once the simulations are completed:
- Rename each created directory to the corresponding model type 
    - e.g. The directory that is created in *Figures/Figure1/CDMetaPOP_inputs/SIR/output/* needs to be renamed to SIR
- Move the renamed folder to *Figures/Figure1/Figure_1_from_source_data/*

After all the data has been added to *Figures/Figure1/Figure_1_from_source_data/*. 
- Navigate to [Figure1](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure1).
    - ``` cd Figures/Figure1``` 
- run
    -``` python Figure_1_from_source.py ```

The image for figure 1 (Figure_1.png) will be created in the *figure_outputs/from_source* directory with in [Figure1](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure1).

### Figure 2

>[!Note]
> ~35 GB will be needed to store the data for this figure

To reproduce the data for figure 2 run:
- ``` python Figure_2_from_source.py ```

The image for figure 1 (Figure_1.png) will be created in the *figure_outputs/from_source* directory with in [Figure2](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure2).

### Figure 3

>[!Note]
> ~35 GB will be needed to store the data for this figure

To reproduce the data for figure 3 run:
- ``` python Figure_3_from_source.py ```

The image for figure 3 (Figure_3.png) will be created in the *figure_outputs/from_source* directory with in [Figure3](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure3).

### Figure 4

>[!Note]
> ~ 8 GB will be needed to store the data for this figure

To reproduce the data for figure 4 run:
- ``` python CDMetaPOP/CDMetaPOP.py Figures/Figure4/CDMetaPOP_inputs/ RunVars.csv output/output ```

Once the simulation is completed:
- Rename the created directory to *Aspatial_inputs*. 
    - e.g. The directory that is created in *Figures/Figure4/CDMetaPOP_inputs/output/* needs to be renamed to *Aspatial_inputs*.
- Move the renamed folder to *Figures/Figure4/Figure_4_from_source_data/*

After the data has been added to *Figures/Figure4/Figure_4_from_source_data/*. 
- Navigate to [Figure4](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure1).
    - ``` setwd("Figures/Figure4") ``` 
- run
    -``` source("Figure_4_from_source.R") ```

The image for figure 4 (Figure_4.png) will be created in the *figure_outputs/from_source* directory with in [Figure4](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure4). 

### Figure 5 and 6

>[!Note]
> ~ 25 GB will be needed to store the data for this figure

Figure 5 and 6 use the same data. 
One can either run the simulation in [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure5) or [Figure6](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure6) and then move the data to the other folder. 
These instructions will instruct the reader to create the simulations in [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure5) and the code for figure 6 will point to that data.
To reproduce the data for figure 5 and 6 run:
- ``` python CDMetaPOP/CDMetaPOP.py Figures/Figure5/CDMetaPOP_inputs/ RunVars.csv output/output ```  

Once the simulation is completed:
- Rename the created directory to *Spatial_inputs*. 
    - e.g. The directory that is created in either *Figures/Figure5/CDMetaPOP_inputs/output/* needs to be renamed to *Spatial_inputs*.
- Move the renamed folder to *Figures/Figure5/Figure_5_from_source_data/*

After the data has been added to either *Figures/Figure5/Figure_5_from_source_data/*. 
- For figure 5:
    - Navigate to the [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure5).
        - ``` setwd("Figures/Figure5") ``` 
    - run
        -``` source("Figure_5_from_source.R") ```

    The image for figure 5 (Figure_5.png) will be created in the *figure_outputs/from_source* directory with in [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure5). 
- For figure 6:
    - Navigate to the [Figure6](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure5).
        - ``` setwd("Figures/Figure6") ``` 
    - run
        -``` source("Figure_6_from_source.R") ```

    The image for figure 6 (Figure_6.png) will be created in the *figure_outputs/from_source* directory with in [Figure6](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figures/Figure6). 





























