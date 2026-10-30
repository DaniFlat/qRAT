# qRAT — qPCR Relative Expression Analysis Tool <img src="www/logo.svg" align="right" height="150" />

[![R-shiny](https://img.shields.io/badge/Shiny-v1.14.0-blue?logo=shiny&logoColor=white)](https://shiny.posit.co/)
[![BioConductor](https://img.shields.io/badge/BioConductor-v3.23-green)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.1186/s12859--022--04823--7-blue)](https://doi.org/10.1186/s12859-022-04823-7)

**qRAT** is a comprehensive, R-based application designed to bridge the gap between raw qPCR data and publication-ready results. It automates the entire data processing, from data parsing and quality control to statistical validation and visualization, requiring **zero programming knowledge**.

---

## 🌐 Web Application

The easiest way to use qRAT is via the hosted web interface:
👉 **[Launch qRAT on shinyapps.io](https://qrat.shinyapps.io/qrat/)**

---

## 💻 Local Installation (using R)

To run qRAT locally, follow these steps:

### 1. Install Dependencies
Ensure you have the latest version of [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/) installed. Open RStudio and run:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

pkgs <- c(
  "shiny", "shinyjs", "shinyWidgets", "shinycssloaders", "shinyjqui", "bslib", "bsicons", "waiter",
  "ggplot2", "ggiraph", "svglite", "dplyr", "tidyr", "data.table", "DT", "rlang", "xfun",
  "ggpubr", "scales", "RColorBrewer", "viridisLite", "ctrlGene", "magrittr", "pkgload",
  "HTqPCR", "ddCt", "limma"
)

BiocManager::install(pkgs)
```

### 2. Install and Run qRAT

```r
shiny::runGitHub("qRAT", "DaniFlat")
```

---

## Example Data

Not sure where to start? Download the [Sample Dataset (.zip)](https://github.com/DaniFlat/qRAT/raw/main/example_data/qRAT_example_data.zip) to explore qRAT's functionalities immediately. The single files and a short description are in the [example_data](example_data) folder.

---

## 📝 Citation

If you use qRAT in your research, please cite our paper:

> Flatschacher, D., Speckbacher, V. & Zeilinger, S. **qRAT: an R-based stand-alone application for relative expression analysis of RT-qPCR data.** *BMC Bioinformatics* 23, 286 (2022). [https://doi.org/10.1186/s12859-022-04823-7](https://doi.org/10.1186/s12859-022-04823-7)

---

## 🛠 Support & Documentation

* **Manual:** The [qRAT User Guide (PDF)](https://github.com/DaniFlat/qRAT/blob/main/www/qRAT_User_Guide.pdf) is also available within the app's "Help" tab.
* **Homepage:** [Official qRAT Website](https://www.uibk.ac.at/microbiology/services/qrat/)
* **Issues:** Found a bug? Please open an [Issue on GitHub](https://github.com/DaniFlat/qRAT/issues).

Developed with ❤️ at the University of Innsbruck.
