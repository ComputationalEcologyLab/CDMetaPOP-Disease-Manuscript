# CDMetaPOP-Disease Manuscript Figure Reproduction

![Compatibility Check](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/actions/workflows/python-check.yml/badge.svg)


### Table of Contents  
[Introduction](#introduction)  
[Software Install](#software-install)  
[Figure 1](#figure-1)  
[Figure 2](#figure-2)  
[Figure 3](#figure-3)  
[Figure 4](#figure-4)  
[Figure 5](#figure-5)  
[Figure 6](#figure-6)  
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


## Figure 1

>[!Important]
> Follow the instructions in [Software Install](#software-install), for python and begin with this repo as your working directory

#### Recreation

To recreate Figure 1:
- Navigate to the [Figure1](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure1) directory within Figure_Reproduction.
    - ``` cd Figure_Reproduction/Figure1``` 
- run
    - ``` python Fig_1_from_data.py ```

The image for figure 1 (Figure_1.png) will be created in the *figure_outputs* directory with in [Figure1](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure1).

#### Reproduction

To reproduce Figure 1:

## Figure 2

>[!Important]
> Follow the instructions in [Software Install](#software-install), for python and begin with this repo as your working directory

#### Recreation

To recreate Figure 2:
- Navigate to the [Figure2](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure2) directory within Figure_Reproduction.
    - ``` cd Figure_Reproduction/Figure2``` 
- run
    - ``` python Fig_2_from_data.py ```

The image for figure 2 (Figure_2.png) will be created in the *figure_outputs* directory with in [Figure2](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure2).

#### Reproduction

To reproduce Figure 2:

## Figure 3

>[!Important]
> Follow the instructions in [Software Install](#software-install), for python and begin with this repo as your working directory

#### Recreation

To recreate Figure 3:
- Navigate to the [Figure3](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure3) directory within Figure_Reproduction.
    - ``` cd Figure_Reproduction/Figure3``` 
- run
    - ``` python Fig_3_from_data.py ```

The image for figure 3 (Figure_3.png) will be created in the *figure_outputs* directory with in [Figure3](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure3).

#### Reproduction

To reproduce Figure 3:

## Figure 4

>[!Important]
> Follow the instructions in [Software Install](#software-install), for R and begin with this repo as your working directory

#### Recreation

To recreate Figure 4:
- Double click the [.Rprofile](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/blob/manuscript_prep/.Rprofile) file, which will open Rstudio
- Navigate to the [Figure4](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure4) directory within Figure_Reproduction.
    - ``` setwd("Figure_Reproduction/Figure4") ``` 
- run in the console
    - ``` source("Fig_4_from_data.R") ```

The image for figure 4 (Figure_4.png) will be created in the *figure_outputs* directory with in [Figure4](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure4).

#### Reproduction

To reproduce Figure 4:

## Figure 5

>[!Important]
> Follow the instructions in [Software Install](#software-install), for R and begin with this repo as your working directory

#### Recreation

To recreate Figure 5:
- Double click the [.Rprofile](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/blob/manuscript_prep/.Rprofile) file, which will open Rstudio
- Navigate to the [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure5) directory within Figure_Reproduction.
    - ``` setwd("Figure_Reproduction/Figure5") ``` 
- run in the console
    - ``` source("Fig_5_from_data.R") ```

The image for figure 5 (Figure_5.png) will be created in the *figure_outputs* directory with in [Figure5](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure5).

#### Reproduction

To reproduce Figure 5:

## Figure 6

>[!Important]
> Follow the instructions in [Software Install](#software-install), for R and begin with this repo as your working directory

#### Recreation

To recreate Figure 6:
- Double click the [.Rprofile](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/blob/manuscript_prep/.Rprofile) file, which will open Rstudio
- Navigate to the [Figure6](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure6) directory within Figure_Reproduction.
    - ``` setwd("Figure_Reproduction/Figure6") ```  
- run in the console
    - ``` source("Fig_6_from_data.R") ```

The image for figure 6 (Figure_6.png) will be created in the *figure_outputs* directory with in [Figure6](https://github.com/ComputationalEcologyLab/CDMetaPOP-Disease-Manuscript/tree/manuscript_prep/Figure_Reproduction/Figure6).

#### Reproduction

To reproduce Figure 6:


