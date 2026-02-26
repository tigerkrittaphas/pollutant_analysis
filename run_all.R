config_env <- new.env()
source("run_config.R", local = config_env)
config <- config_env$config

if (is.null(config)) {
  stop("run_config.R did not define a `config` list.")
}

if (!is.null(config$project_root)) {
  setwd(config$project_root)
}

options(
  hf_pm_analysis.root = getwd(),
  hf_pm_analysis.data_path = file.path(getwd(), config$data_path)
)

if (length(config$experiments) == 0) {
  stop("No experiments defined in run_config.R")
}

default_exp <- config$experiments[[1]]
options(
  hf_pm_analysis.pollutant = default_exp$pollutant,
  hf_pm_analysis.outcome = default_exp$outcome,
  hf_pm_analysis.lag = default_exp$lag
)

if (isTRUE(config$basic_plots$run)) {
  options(hf_pm_analysis.basic_output_dir = config$basic_plots$output_dir)
  source("plot_basic.R", local = new.env(parent = globalenv()))
}

if (isTRUE(config$thai_map$run)) {
  options(hf_pm_analysis.map_output_dir = config$thai_map$output_dir)
  source("plot_thai_map.R", local = new.env(parent = globalenv()))
}

for (exp in config$experiments) {
  options(
    hf_pm_analysis.pollutant = exp$pollutant,
    hf_pm_analysis.outcome = exp$outcome,
    hf_pm_analysis.lag = exp$lag
  )

  if (isTRUE(config$main$run)) {
    options(hf_pm_analysis.main_output_dir = config$main$output_dir)
    source("main.R", local = new.env(parent = globalenv()))
  }

  if (isTRUE(config$subgroup$run)) {
    subgroup_columns <- config$subgroup$subgroup_columns
    if (is.null(subgroup_columns) || length(subgroup_columns) == 0) {
      subgroup_columns <- "region"
    }

    for (subgroup_col in subgroup_columns) {
      subgroup_prefix <- switch(
        subgroup_col,
        region = "regional",
        nhso_region_code = "nhso_region",
        subgroup_col
      )

      subgroup_output_dir <- file.path(
        config$subgroup$output_dir,
        paste0(subgroup_prefix, "_", exp$pollutant, "-", exp$outcome, "-", exp$lag, "lag")
      )
      overall_results_dir <- file.path(
        "output",
        config$main$output_dir,
        paste0("results_", exp$pollutant, "-", exp$outcome, "-", exp$lag, "lag")
      )

      options(
        hf_pm_analysis.subgroup_output_dir = subgroup_output_dir,
        hf_pm_analysis.subgroup_column = subgroup_col,
        hf_pm_analysis.overall_results_dir = overall_results_dir
      )
      source("subgroup.R", local = new.env(parent = globalenv()))

      if (isTRUE(config$subgroup_plot$run)) {
        input_csv <- file.path(
          subgroup_output_dir,
          paste0("pooled_lagged_RR_by_", subgroup_col, ".csv")
        )
        plot_output_dir <- file.path(
          config$subgroup_plot$output_dir,
          paste0(subgroup_prefix, "_", exp$pollutant, "-", exp$outcome, "-", exp$lag, "lag")
        )

        options(
          hf_pm_analysis.subgroup_plot_input = input_csv,
          hf_pm_analysis.subgroup_plot_col = "Subgroup",
          hf_pm_analysis.subgroup_plot_output_dir = plot_output_dir,
          hf_pm_analysis.subgroup_plot_max_lag = config$subgroup_plot$max_lag,
          hf_pm_analysis.subgroup_plot_title = paste("RR association with a 10-unit increase in", exp$pollutant),
          hf_pm_analysis.subgroup_plot_subtitle = paste("RR association by", subgroup_col),
          hf_pm_analysis.subgroup_plot_prefix = paste0("subgroup_RR_lagged_", subgroup_col)
        )
        source("subgroup_plot.R", local = new.env(parent = globalenv()))
      }
    }
  }
}
