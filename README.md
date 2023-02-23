# qRAT - qPCR relative expression analysis tool

## What's qRAT?

qRAT is a R based standalone desktop application to automate the processing of raw Quantification Cycle (Cq) data files exported from virtually any qPCR instrument using well established and state-of-the-art statistical and graphical techniques. The purpose of this tool is to provide a comprehensive, straightforward, and easy-to-use solution for the relative quantification of RT-qPCR data that requires no programming knowledge or additional software installation.

The current implementation allows ΔCq calculation (relative to endogenous control(s)), ΔΔCq calculation (relative to endogenous control(s) and a reference sample) and inter-plate variation correction. Moreover, functionalities for parsing, filtering and visualisation of relative RT-qPCR data are included.

This repository contains the R scripts that are used to create the shiny app.

## How to run qRAT on Windows?

qRAT is distributed as a Windows desktop application with Electron using the RInno Package. The setup can be easily downloaded and installed from the [qRAT Homepage](https://www.uibk.ac.at/microbiology/services/qrat/).

## How to run qRAT on Linux and MacOS?

Users can run the app using R console and RStudio. 
1. Upgrade to the most recent version of R and Rstudio
2. Start RStudio and install all the bioconductor packages by entering the following commands in the console
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HTqPCR")
BiocManager::install("ddCt")
BiocManager::install("limma")
```
3. Install the qRAT package with devtools from github by executing the following commands
```
## install.packages('devtools')
require('devtools')
install_github('DaniFlat/qRAT')
```
4. Load the libraries and start the application with the following commands
```
library(HTqPCR)
library(ddCt)
library(limma)
library(qRAT)
library(scales)
library(data.table)
library(DT)
library(dplyr)
library(waiter)
library(tidyr)
library(stringr)
library(magrittr)
library(shinycssloaders)
library(curl)
library(viridisLite)
library(shinyjqui)
library(plotly)
library(shinyjs)
library(ggplot2)
library(shinyWidgets)
qRAT()
```

## Web Application

qRAT is available as web application on https://qrat.shinyapps.io/qrat/

## Citation

If you use qRAT for any published research, please include the following citation:

Flatschacher, D., Speckbacher, V. & Zeilinger, S. qRAT: an R-based stand-alone application for relative expression analysis of RT-qPCR data. BMC Bioinformatics 23, 286 (2022). https://doi.org/10.1186/s12859-022-04823-7

## Getting Help
[qRAT Homepage](https://www.uibk.ac.at/microbiology/services/qrat/)
