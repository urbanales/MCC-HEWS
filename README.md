# The effectiveness of heat prevention plans in reducing heat-related mortality across Europe

This repository contains the full reproducible R code pipeline for the study: **"The effectiveness of heat prevention plans in reducing heat-related mortality across Europe."** The analysis is structured in four main scripts, each corresponding to a major stage of the research workflow.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Usage Instructions](#usage-instructions)
- [Script Details](#script-details)
    - [01. Data Preparation](#01-data-preparation)
    - [02. Stage 2 Meta-Regression](#02-stage-2-meta-regression)
    - [03. Stage 3 Attributable Mortality Calculation](#03-stage-3-attributable-mortality-calculation)
    - [04. Final Analysis and Reporting](#04-final-analysis-and-reporting)
- [Outputs](#outputs)
- [Authors](#authors)
- [License](#license)

---

## Overview

This project assesses the impact of heat prevention plans (heat warning systems, HWS) on heat-related mortality across European cities, following a multi-stage, meta-regression based analytical pipeline. Each R script corresponds to a major step in the analysis, from data preparation to final reporting and visualization.

---

## Repository Structure

```
.
├── 01.stage1_data_preparation.R
├── 02.stage2_meta_regression.R
├── 03.stage3_attributable_mortality_calculation.R
├── 04.stage4_final_analysis_and_reporting.R
├── data/
│   └── eu_countries.csv
├── output/
│   └── (tables, figures, and summary files)
└── README.md
```

---

## Prerequisites

- **R version**: ≥ 4.1.0
- **Required R packages**:  
  `dlnm`, `MASS`, `dplyr`, `lubridate`, `pracma`, `readr`, `tidyverse`, `ggplot2`, `scales`, `ggh4x`, `RColorBrewer`
- **Input data**: See `data/` directory and script-specific comments for input requirements.

**Note on Raw Dataset Availability:**  
The main dataset of daily mortality, temperature, and heat warning system status is not publicly available due to privacy and data sharing restrictions. Only aggregated data objects and code are provided here to ensure full reproducibility of all analyses presented in the study.

---

## Usage Instructions

1. **Clone or download this repository.**
2. **Install required R packages** (if not already installed).
3. **Prepare input data** as described in the scripts and place in the `data/` directory.
4. **Run each script in order:**
    - Open each `.R` script in R (preferably RStudio).
    - Source all previous scripts before running the next.
    - Modify scenario parameters as needed (see comments in scripts).
5. **Outputs** (tables, figures, CSVs) will be saved in the `output/` directory.

---

## Script Details

### 01. Data Preparation

**File:** `01.stage1_data_preparation.R`

- Loads and prepares city-level and national-level data.
- Creates data objects used in subsequent stages (`eu_cities`, `dlist`, `hws_ind`, etc.).
- **Note:** The raw daily data are not included in this repository due to privacy restrictions; only pre-processed and aggregated data objects necessary for full reproducibility are provided.

---

### 02. Stage 2 Meta-Regression

**File:** `02.stage2_meta_regression.R`

- Performs first- and second-stage meta-analytical modeling.
- Produces BLUPs (Best Linear Unbiased Predictors) and meta-regression outputs for use in stage 3.

---

### 03. Stage 3 Attributable Mortality Calculation

**File:** `03.stage3_attributable_mortality_calculation.R`

- Uses BLUPs to calculate heat-attributable mortality numbers/fractions (AN, AF) for each city, country, region, and EU as a whole.
- Runs scenario simulations (counterfactual and factual) with uncertainty quantification.
- Outputs simulation results for further aggregation.

---

### 04. Final Analysis and Reporting

**File:** `04.stage4_final_analysis_and_reporting.R`

- Aggregates simulation results to produce summary tables for countries, regions, HWS classes, and the EU.
- Generates publication-ready figures and summary CSVs.
- Provides code for visualizing main results and supplementary figures.

---

## Outputs

- **Tables:** CSV summary tables at country, region, HWS class, and EU level.
- **Figures:** Publication-ready plots (PDF/PNG) illustrating heat-attributable fractions and differences under different scenarios.
- **Intermediate RData files:** For reproducibility and further analysis.

Output files are saved in the `output/` folder. See script comments for file naming conventions.

---

## Authors

- Aleš Urban
- Veronika Huber
- Pierre Masselot
- Antonio Gasparrini

For questions, please contact the corresponding author.

---

## License

This project is released under the MIT License. See `LICENSE` file for details.

---
