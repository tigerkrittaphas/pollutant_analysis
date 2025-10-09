# Heart Failure and PM2.5 Analysis in Thailand

This repository contains code and data for analyzing the association between air pollution (PM2.5 and other pollutants) and heart failure hospitalizations across provinces and subgroups in Thailand.

## Project Structure

- `main.R`: Main analysis script for province-level DLNM modeling and meta-analysis.
- `subgroup.R`: Subgroup analysis pipeline (e.g., by region, NHSO region code, gender, age).
- `subgroup_plot.R`: Visualization of subgroup results (gender, age).
- `plot_basic.R`: Basic exploratory plots and summary statistics.
- `plot_thai_map.R`: Geospatial visualization of hospitalization and PM2.5 levels by province.
- `data/`: Contains input CSV files (weather, hospitalization, NDVI).
- `output/`: Stores all analysis results, plots, and tables.
- `renv/`: R environment management (dependencies).
- `README.md`: This file.

## Analysis Workflow

1. **Data Preparation**  
   - Input data: `data/daily_all_weather_30jul.csv`  
   - Data cleaning and formatting in [`main.R`](main.R) and [`plot_basic.R`](plot_basic.R).

2. **Main Analysis**  
   - Province-specific time-series modeling using DLNM (Distributed Lag Non-linear Models).
   - Meta-analysis to pool province estimates ([`main.R`](main.R)).
   - Outputs: Relative Risk estimates, forest plots, lag-response curves.

3. **Subgroup Analysis**  
   - Analysis by region, NHSO region code, gender, age ([`subgroup.R`](subgroup.R)).
   - Results saved in `output/regional_subgroup_analysis_30jul/`.

4. **Visualization**  
   - Forest plots and lag-response plots for subgroups ([`subgroup_plot.R`](subgroup_plot.R)).
   - Geospatial maps of hospitalization and PM2.5 ([`plot_thai_map.R`](plot_thai_map.R)).
   - Summary statistics and correlation matrices ([`plot_basic.R`](plot_basic.R)).

## How to Run

1. **Set up R environment**
   - Use [renv](https://rstudio.github.io/renv/) for reproducible package management.
   - Run `source("renv/activate.R")` or start R in this directory.

2. **Run analyses**
   - Execute scripts in order:
     - `main.R` for main analysis.
     - `subgroup.R` for subgroup analysis.
     - `subgroup_plot.R` for subgroup visualizations.
     - `plot_basic.R` and `plot_thai_map.R` for exploratory and geospatial plots.

3. **View results**
   - All outputs (plots, tables, CSVs) are saved in the `output/` directory.

## Key Files and Functions

- [`main.R`](main.R): Main DLNM analysis and meta-analysis.
- [`subgroup.R`](subgroup.R): [`run_analysis_for_subgroup`](subgroup.R) function for subgroup DLNM analysis.
- [`subgroup_plot.R`](subgroup_plot.R): [`plot_across_subgroup`](subgroup_plot.R) function for subgroup RR visualization.
- [`plot_basic.R`](plot_basic.R): Summary statistics, correlation, and time series plots.
- [`plot_thai_map.R`](plot_thai_map.R): Province-level choropleth maps.

## Data Sources

- Hospitalization and weather data: `data/daily_all_weather_30jul.csv`
- NDVI data: `data/NDVI2016-2019.csv`

## Output

- All results, plots, and tables are saved in the `output/` directory, organized by analysis type and date.

## Contact

For questions or collaboration, please contact the repository owner.

---
