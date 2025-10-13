#!/usr/bin/env Rscript

# Phase 04 - Step 02b: Pyoverdine Degree Distributions
# Generate in-degree (production) and out-degree (utilization) distributions
# Outputs: figures/network_topology/pyov_...pdf/png

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

plot_degree_hist <- function(deg_vec, title, xlabel, ymax = NULL) {
  library(dplyr)
  library(ggplot2)

  # 自定义分箱
  bins <- seq(floor(log10(min(deg_vec))), ceiling(log10(max(deg_vec))), by = 0.2)
  df <- data.frame(degree = deg_vec) %>%
    mutate(bin = cut(log10(degree), breaks = bins, include.lowest = TRUE)) %>%
    count(bin) %>%
    mutate(
      bin_mid = 10 ^ (sapply(strsplit(as.character(bin), ","), function(x)
        mean(as.numeric(gsub("[^0-9.-]", "", x)))))
    )

  p <- ggplot(df, aes(x = bin_mid, y = n)) +
    geom_col(
      width = 0.15,  # 控制每个bar在log尺度下的视觉宽度
      fill = "steelblue",
      color = "white",
      alpha = 0.85
    ) +
    scale_x_log10() +
    labs(
      title = title,
      x = xlabel,
      y = "No. of siderophore group"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )

  if (!is.null(ymax)) p <- p + coord_cartesian(ylim = c(0, ymax))
  return(p)
}



plot_ccdf <- function(deg_vec, title, xlabel) {
  nonzero <- deg_vec[deg_vec > 0]
  if (length(nonzero) == 0) {
    return(ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "No non-zero degrees"))
  }
  sorted <- sort(nonzero, decreasing = TRUE)
  n <- length(sorted)
  ccdf <- data.frame(degree = sorted, prob = (1:n)/n)
  ggplot(ccdf, aes(x = degree, y = prob)) +
    geom_line(color = "red") +
    geom_point(size = 0.6, color = "red", alpha = 0.6) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = title, x = paste0(xlabel, " (log10)"), y = "P(X ≥ degree) (log10)") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
}

safe_dir_create("figures/network_topology")
safe_dir_create("figures/network_topology/pyov_degree")

# Load degree metrics
rds_fg <- "results/phase_04/degree/degree_metrics_fg.rds"
rds_str <- "results/phase_04/degree/degree_metrics_str.rds"
if (!file.exists(rds_fg) || !file.exists(rds_str)) {
  stop("Degree metrics RDS files not found. Run scripts/phase_04_topology/02_degree_and_distributions.R first.")
}
fg <- readRDS(rds_fg)
strd <- readRDS(rds_str)

# Extract pyoverdine degree vectors
pyov_in_fg <- as.integer(fg$pyov_prod_in)
names(pyov_in_fg) <- names(fg$pyov_prod_in)
pyov_out_fg <- as.integer(fg$pyov_util_out)
names(pyov_out_fg) <- names(fg$pyov_util_out)

pyov_in_str <- as.integer(strd$pyov_prod_in)
names(pyov_in_str) <- names(strd$pyov_prod_in)
pyov_out_str <- as.integer(strd$pyov_util_out)
names(pyov_out_str) <- names(strd$pyov_util_out)


# FG-level plots
p_fg_in_hist <- plot_degree_hist(pyov_in_fg, "Pyoverdine Production In-Degree (FG level)", "in degree")
p_fg_out_hist <- plot_degree_hist(pyov_out_fg, "Pyoverdine Utilization Out-Degree (FG level)", "out degree")
p_fg_in_ccdf <- plot_ccdf(pyov_in_fg, "Pyoverdine Production In-Degree CCDF (FG)", "Production in-degree")
p_fg_out_ccdf <- plot_ccdf(pyov_out_fg, "Pyoverdine Utilization Out-Degree CCDF (FG)", "Utilization out-degree")

ggsave("figures/network_topology/pyov_degree/pyov_degree_in_fg.png", p_fg_in_hist, width = 6, height = 5)
ggsave("figures/network_topology/pyov_degree/pyov_degree_out_fg.png", p_fg_out_hist, width = 6, height = 5)
ggsave("figures/network_topology/pyov_degree/pyov_degree_in_fg_ccdf.png", p_fg_in_ccdf, width = 6, height = 5)
ggsave("figures/network_topology/pyov_degree/pyov_degree_out_fg_ccdf.png", p_fg_out_ccdf, width = 6, height = 5)

# STR-level plots
p_str_in_hist <- plot_degree_hist(pyov_in_str, "Pyoverdine Production In-Degree (STR level)", "in degree")
p_str_out_hist <- plot_degree_hist(pyov_out_str, "Pyoverdine Utilization Out-Degree (STR level)", "out degree")
p_str_in_ccdf <- plot_ccdf(pyov_in_str, "Pyoverdine Production In-Degree CCDF (STR)", "Production in-degree")
p_str_out_ccdf <- plot_ccdf(pyov_out_str, "Pyoverdine Utilization Out-Degree CCDF (STR)", "Utilization out-degree")

ggsave("figures/network_topology/pyov_degree/pyov_degree_in_str.png", p_str_in_hist, width = 6, height = 5)
ggsave("figures/network_topology/pyov_degree/pyov_degree_out_str.png", p_str_out_hist, width = 6, height = 5)
ggsave("figures/network_topology/pyov_degree/pyov_degree_in_str_ccdf.png", p_str_in_ccdf, width = 6, height = 5)
ggsave("figures/network_topology/pyov_degree/pyov_degree_out_str_ccdf.png", p_str_out_ccdf, width = 6, height = 5)


# Combined figure: FG CCDF in/out side-by-side (for manuscript panels B/C)
combined_fg_ccdf <- (p_fg_in_ccdf + p_fg_out_ccdf) + plot_layout(ncol = 2) + plot_annotation(title = "Pyoverdine Degree CCDFs (FG level)")
ggsave("figures/network_topology/pyov_degree/pyov_degree_ccdf_fg_combined.png", combined_fg_ccdf, width = 12, height = 5)

combined_str_ccdf <- (p_str_in_ccdf + p_str_out_ccdf) + plot_layout(ncol = 2) + plot_annotation(title = "Pyoverdine Degree CCDFs (STR level)")
ggsave("figures/network_topology/pyov_degree/pyov_degree_ccdf_str_combined.png", combined_str_ccdf, width = 12, height = 5)


# Save a small CSV summary for downstream annotation (pyov name, in_fg, out_fg, in_str, out_str)
pyov_names <- names(pyov_in_fg)
if (is.null(pyov_names)) pyov_names <- seq_along(pyov_in_fg)
summary_df <- data.frame(
  pyov = pyov_names,
  in_fg = pyov_in_fg,
  out_fg = pyov_out_fg,
  in_str = pyov_in_str[pyov_names],
  out_str = pyov_out_str[pyov_names],
  stringsAsFactors = FALSE
)
safe_dir_create("results/phase_04/degree")
write.csv(summary_df, "results/phase_04/degree/pyov_degree_summary.csv", row.names = FALSE)

cat("Pyoverdine degree distribution figures and summary saved under figures/network_topology/pyov_degree and results/phase_04/degree/pyov_degree_summary.csv\n")