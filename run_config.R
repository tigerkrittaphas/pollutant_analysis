base_output_dir <- "output/16Feb_all"

config <- list(
  project_root = getwd(),
  data_path = "data/daily_all_weather_30jul.csv",
  experiments = list(
    list(pollutant = "PM2.5", outcome = "hf_prim", lag = 7),
    list(pollutant = "PM2.5", outcome = "hf_prim_male", lag = 7),
    list(pollutant = "PM2.5", outcome = "hf_prim_female", lag = 7),
    list(pollutant = "PM2.5", outcome = "hf_prim_under65", lag = 7),
    list(pollutant = "PM2.5", outcome = "hf_prim_over65", lag = 7),
    list(pollutant = "PM2.5", outcome = "af_prim", lag = 7),
    list(pollutant = "PM2.5", outcome = "af_prim_male", lag = 7),
    list(pollutant = "PM2.5", outcome = "af_prim_female", lag = 7),
    list(pollutant = "PM2.5", outcome = "af_prim_under65", lag = 7),
    list(pollutant = "PM2.5", outcome = "af_prim_over65", lag = 7),
    list(pollutant = "PM2.5", outcome = "hf", lag = 7),
    list(pollutant = "PM2.5", outcome = "hf_other", lag = 7),
    list(pollutant = "PM2.5", outcome = "af", lag = 7),
    list(pollutant = "PM2.5", outcome = "af_other", lag = 7)
  ),
  main = list(
    run = TRUE,
    output_dir = basename(base_output_dir)
  ),
  basic_plots = list(
    run = TRUE,
    output_dir = base_output_dir
  ),
  thai_map = list(
    run = TRUE,
    output_dir = base_output_dir
  ),
  subgroup = list(
    run = TRUE,
    output_dir = base_output_dir,
    subgroup_columns = c("region", "nhso_region_code")
  ),
  subgroup_plot = list(
    run = TRUE,
    max_lag = 4,
    output_dir = base_output_dir
  )
)
