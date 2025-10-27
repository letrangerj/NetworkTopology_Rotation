#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

metrics_path <- "results/phase_04/syn_rec_ratio/pyov_metrics_simplified.csv"
output_dir <- "figures/network_topology/syn_rec_ratio/figures"

if (!file.exists(metrics_path)) {
  stop("Missing syn/rec metrics. Run 06_syn_rec_ratio.R first.")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

highlight_ids <- c(
  "PYO_27", "PYO_17", "PYO_18", "PYO_93",
  "PYO_36", "PYO_37", "PYO_38", "PYO_39", "PYO_40", "PYO_41", "PYO_42", "PYO_43", "PYO_44", "PYO_45",
  "PYO_32", "PYO_33", "PYO_21", "PYO_24", "PYO_26", "PYO_19", "PYO_34", "PYO_25",
  "PYO_12", "PYO_13"
)

highlight_color <- "#d1495b"

metrics <- read_csv(metrics_path, show_col_types = FALSE) %>%
  mutate(is_highlight = pyo_id %in% highlight_ids)

make_scatter <- function(df, x_col, y_col, x_label, y_label, title, subtitle, filename) {
  base_df <- df %>% filter(.data[[x_col]] > 0, .data[[y_col]] > 0)
  highlight_df <- filter(base_df, is_highlight)

  plot <- ggplot(base_df, aes_string(x = x_col, y = y_col)) +
    geom_point(color = "grey75", size = 2, alpha = 0.85) +
    geom_point(data = highlight_df, color = highlight_color, size = 2.8) +
    geom_text(
      data = highlight_df,
      aes(label = pyo_id),
      size = 3,
      vjust = -0.7,
      color = highlight_color,
      check_overlap = TRUE
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    scale_x_log10(name = x_label) +
    scale_y_log10(name = y_label) +
    labs(title = title, subtitle = subtitle) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  outfile <- file.path(output_dir, filename)
  ggsave(outfile, plot, width = 7, height = 5, dpi = 300)
  message("Saved syn/rec scatter: ", outfile)
}

make_scatter(
  metrics,
  x_col = "fg_utilization_count",
  y_col = "fg_production_count",
  x_label = "Functional-group receptors (utilization count)",
  y_label = "Functional-group synthetases (production count)",
  title = "Pyoverdine syn/rec balance (FG counts)",
  subtitle = "x: receptor (utilizer) groups, y: synthetase (producer) groups; dashed line = 1:1",
  filename = "06a_syn_rec_scatter_fg.png"
)

make_scatter(
  metrics,
  x_col = "strain_utilization_count",
  y_col = "strain_production_count",
  x_label = "Strain receptors (utilization count)",
  y_label = "Strain synthetases (production count)",
  title = "Pyoverdine syn/rec balance (strain counts)",
  subtitle = "x: receptor (utilizer) strains, y: synthetase (producer) strains; dashed line = 1:1",
  filename = "06a_syn_rec_scatter_strain.png"
)
