# This file is for plotting the subgroup analysis results for easier visualization.

library(ggplot2)
library(dplyr)
library(rlang)

plot_across_subgroup <- function(
    df,
    subgroup_col,
    max_lag = 4,
    title = "RR association with a 10-unit increase in PM2.5",
    subtitle = "RR association by subgroup") {
  subgroup_sym <- rlang::sym(subgroup_col)

  rr_line_plot <- df %>%
    mutate(across(
      c("Lag", "RR", "Lower", "Upper"),
      as.numeric
    )) %>%
    filter(Lag < max_lag) %>%
    ggplot(
      aes(
        x = Lag,
        y = RR,
        ymin = Lower,
        ymax = Upper,
        color = !!subgroup_sym,
        fill = !!subgroup_sym,
        group = !!subgroup_sym
      )
    ) +
    geom_errorbar(width = 0.2, linewidth = 0.8, alpha = 0.4, position = "dodge") +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    scale_x_continuous(breaks = 0:(max_lag - 1)) +
    labs(
      subtitle = subtitle,
      title = title,
      x = "Lag",
      y = "RR and 95% CI",
      color = "Subgroup",
      fill = "Subgroup"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  return(rr_line_plot)
}

run_subgroup_plot <- function(
    input_csv,
    subgroup_col,
    output_dir,
    max_lag = 4,
    title = "RR association with a 10-unit increase in PM2.5",
    subtitle = "RR association by subgroup",
    output_prefix = "subgroup_RR_lagged_plot") {
  df <- read.csv(input_csv)
  plot_obj <- plot_across_subgroup(df, subgroup_col, max_lag, title, subtitle)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(output_dir, paste0(output_prefix, ".png")),
    plot = plot_obj,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(output_prefix, ".svg")),
    plot = plot_obj,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )

  invisible(plot_obj)
}

if (sys.nframe() == 0) {
  input_csv <- getOption("hf_pm_analysis.subgroup_plot_input", NA_character_)
  subgroup_col <- getOption("hf_pm_analysis.subgroup_plot_col", "Subgroup")
  output_dir <- getOption("hf_pm_analysis.subgroup_plot_output_dir", "output/subgroup_plots")
  max_lag <- getOption("hf_pm_analysis.subgroup_plot_max_lag", 4)
  title <- getOption("hf_pm_analysis.subgroup_plot_title", "RR association with a 10-unit increase in PM2.5")
  subtitle <- getOption("hf_pm_analysis.subgroup_plot_subtitle", "RR association by subgroup")
  output_prefix <- getOption("hf_pm_analysis.subgroup_plot_prefix", "subgroup_RR_lagged_plot")

  if (!is.na(input_csv) && file.exists(input_csv)) {
    run_subgroup_plot(
      input_csv = input_csv,
      subgroup_col = subgroup_col,
      output_dir = output_dir,
      max_lag = max_lag,
      title = title,
      subtitle = subtitle,
      output_prefix = output_prefix
    )
  } else {
    stop("subgroup_plot.R: input CSV not found. Set hf_pm_analysis.subgroup_plot_input.")
  }
}
