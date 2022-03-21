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
2. Start RStudio and install all the R packages
```
listOfPackages <- c("bslib","DT","data.table","dplyr","ggplot2",
                    "plotly","reshape2","scales","shiny",
                    "shinycssloaders","shinyWidgets","shinyjs","thematic","waiter","xtable")
for (i in listOfPackages){
     if(! i %in% installed.packages()){
         install.packages(i, dependencies = TRUE)
     }
     require(i)
}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HTqPCR")
BiocManager::install("ddCt")
```
3. Launch the app from R/Rstudio, either paste this in the command line:
```
library(shiny)
shiny::runGitHub('qRAT', 'DaniFlat')
```
or download the latest release and place the files into an app directory in your working directory and launch the app with
```
library(shiny)
runApp("qRAT")
```

## Citation

If you use qRAT for any published research, please include the following citation:

## Getting Help
[qRAT Homepage](https://www.uibk.ac.at/microbiology/services/qrat/)
