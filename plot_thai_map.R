set.seed(123)

project_root <- getOption("hf_pm_analysis.root", getwd())
data_path <- getOption(
  "hf_pm_analysis.data_path",
  file.path(project_root, "data", "daily_all_weather_30jul.csv")
)
pollutant_col_name <- getOption("hf_pm_analysis.pollutant", "PM2.5")

# Load necessary libraries
library(ggplot2)
library(sf)
library(dplyr)
library(geodata)

dat <- read.csv(data_path)

# properly format date
dat$date_start <- as.Date(dat$date_start, "%Y-%m-%d")
dat$dow <- as.factor(weekdays(dat$date_start))
dat$month <- as.factor(months(dat$date_start))
dat$year <- as.factor(format(dat$date_start, format = "%Y"))

# get strata for case crossover analysis
dat$stratum <- as.factor(dat$year:dat$month:dat$dow)

results_dir <- getOption("hf_pm_analysis.map_output_dir", "output/demographic_pollutants")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

invalid_prov <- c("Saraburi", "Samut Sakhon", "Nonthaburi")
dat <- dat %>%
  filter(!prov_name %in% invalid_prov)

hospitalization_vars <- c("hf_prim", "af_prim", "af", "af_other", "hf", "hf_other")
hospitalization_labels <- c(
  hf_prim = "Heart Failure (Primary)",
  af_prim = "Atrial Fibrillation (Primary)",
  af = "Atrial Fibrillation",
  af_other = "Atrial Fibrillation (Other)",
  hf = "Heart Failure",
  hf_other = "Heart Failure (Other)"
)

# =================
# Download province-level (level=1) data for Thailand
thailand_provinces <- gadm(country = "THA", level = 1, path = tempdir()) |>
  st_as_sf()

for (var in hospitalization_vars) {
  var_label <- hospitalization_labels[[var]]
  count_col <- paste0(var, "_count")

  avg_hosp_count <- dat %>%
    group_by(prov_name) %>%
    summarise(
      value = mean(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    )

  thailand_map_with_hosp <- thailand_provinces |>
    left_join(avg_hosp_count, by = c("NAME_1" = "prov_name"))

  ggplot(data = thailand_map_with_hosp) +
    geom_sf(aes(fill = value), color = "white", size = 0.1) +
    scale_fill_gradient(
      low = "white",
      high = "blue",
      na.value = "darkgrey"
    ) +
    labs(
      title = paste("Average Daily Hospitalization by", var_label, "by Province in Thailand"),
      fill = "Hospitalization Count",
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(results_dir, paste0("thai_map_", var, ".png")),
    height = 8, width = 10, dpi = 300, units = "in",
    bg = "white"
  )
  ggsave(
    filename = file.path(results_dir, paste0("thai_map_", var, ".svg")),
    height = 8, width = 10, dpi = 300, units = "in",
    bg = "white"
  )

  thailand_map_with_hosp_save <- thailand_map_with_hosp %>%
    mutate(!!count_col := value) %>%
    select(NAME_1, VARNAME_1, NL_NAME_1, all_of(count_col)) %>%
    select(-geometry)

  write.csv(
    thailand_map_with_hosp_save,
    file = file.path(results_dir, paste0("thai_map_", var, ".csv")),
    row.names = FALSE
  )
}

# =============================================
# Create average map of selected pollutant by provinces
avg_prov_poll <- dat %>%
  group_by(prov_name) %>%
  summarise(
    average_pollutant = mean(.data[[pollutant_col_name]], na.rm = TRUE),
    .groups = "drop"
  )

##  Merge Spatial Data with Your Data
thailand_map_with_pollutant <- thailand_provinces |>
  left_join(avg_prov_poll, by = c("NAME_1" = "prov_name"))

##  Plot the Heat Map (Choropleth)
ggplot(data = thailand_map_with_pollutant) +
  geom_sf(aes(fill = average_pollutant), color = "white", size = 0.1) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    na.value = "darkgrey"
  ) + # Use a nice color scale
  labs(
    title = paste("Average", pollutant_col_name, "levels by Province in Thailand"),
    fill = paste("Average", pollutant_col_name),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(), # Removes both x and y axis titles
    axis.text = element_blank(), # Removes the axis text (numbers)
    axis.ticks = element_blank(), # Removes the axis tick marks
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(results_dir, paste0("thai_map_", pollutant_col_name, ".png")),
  height = 8, width = 10, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, paste0("thai_map_", pollutant_col_name, ".svg")),
  height = 8, width = 10, dpi = 300, units = "in",
  bg = "white"
)

# Save the aggregated data to csv
thailand_map_with_pollutant_save <- thailand_map_with_pollutant %>%
  select(NAME_1, VARNAME_1, NL_NAME_1, average_pollutant) %>%
  select(-geometry)

write.csv(thailand_map_with_pollutant_save,
  file = file.path(results_dir, paste0("thai_map_", pollutant_col_name, ".csv")),
  row.names = FALSE
)
