#!/usr/bin/env Rscript

# Phase 3 - Step 4: Network Visualization (Fig 3A style)
# Rewrite: shorten edges so arrowheads stop outside nodes; correct draw order and legend placement
# Dependencies: igraph, ggplot2 (no tidyverse piping required)

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
})

set.seed(2025)

# -----------------------------------------------------------------------------
# Paths and I/O
# -----------------------------------------------------------------------------
img_dir   <- "figures/network"
docs_dir  <- "docs/phase_03"

if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(docs_dir)) dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

edges_csv <- "data/interim/edges_functional_groups_conservative.csv"
pyov_csv  <- "data/interim/nodes_pyoverdines_conservative.csv"
fg_csv    <- "data/interim/nodes_functional_groups_conservative.csv"

edges <- read.csv(edges_csv, stringsAsFactors = FALSE)
pyov  <- read.csv(pyov_csv,  stringsAsFactors = FALSE)
fg    <- read.csv(fg_csv,    stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Build vertex and edge tables for igraph
# -----------------------------------------------------------------------------
# Pyoverdine nodes
vertex_pyov <- data.frame(
  name  = pyov$node_id,
  type  = "pyov",
  color = "orange",
  size  = pmax(2, sqrt(pyov$matrix_index) * 1.4),
  stringsAsFactors = FALSE
)

# Functional group nodes: color by strategy class
# - Single-receptor producer: green
# - Multi-producer: yellow
# - Nonproducer: red
# - Fallback/other: gray80
vertex_fg <- data.frame(
  name  = fg$agent_id,
  type  = "fg",
  color = ifelse(
    fg$strategy == "Single-receptor producer", "green",
    ifelse(
      fg$strategy == "Multi-producer", "yellow",
      ifelse(fg$strategy == "Nonproducer", "red", "gray80")
    )
  ),
  size  = pmax(3, sqrt(fg$strain_count) * 2),
  stringsAsFactors = FALSE
)

vertex_df <- rbind(vertex_pyov, vertex_fg)
edge_data <- data.frame(
  from = edges$source,
  to   = edges$target,
  edge_type = edges$edge_type,
  stringsAsFactors = FALSE
)

# Build igraph
g <- graph_from_data_frame(d = edge_data, directed = TRUE, vertices = vertex_df)

# -----------------------------------------------------------------------------
# Layout and node coordinates
# -----------------------------------------------------------------------------
set.seed(42)
coords <- layout_with_fr(g, niter = 500)

node_df <- data.frame(
  name  = V(g)$name,
  x     = coords[, 1],
  y     = coords[, 2],
  type  = as.character(V(g)$type),
  color = as.character(V(g)$color),
  size  = as.numeric(V(g)$size),
  stringsAsFactors = FALSE
)

# Merge edge coords (from/to)
edge_coords <- merge(edge_data, node_df, by.x = "from", by.y = "name")
edge_coords <- merge(
  edge_coords,
  node_df,
  by.x = "to",
  by.y = "name",
  suffixes = c("", "_to")
)
# Now edge_coords has columns:
#  from, to, edge_type, x, y, type, color, size, x_to, y_to, type_to, color_to, size_to

# -----------------------------------------------------------------------------
# Shorten edges so arrowheads stop outside target nodes
# -----------------------------------------------------------------------------
# We shorten both start and end points along the segment by a fraction of the
# center-to-center distance. This avoids arrowheads being covered by node markers.
shorten_segments <- function(df,
                             start_x = "x", start_y = "y",
                             end_x = "x_to", end_y = "y_to",
                             start_shrink = 0.03,
                             end_shrink   = 0.10) {
  dx <- df[[end_x]] - df[[start_x]]
  dy <- df[[end_y]] - df[[start_y]]
  # Distance
  d  <- sqrt(dx * dx + dy * dy)
  d[d == 0] <- 1e-9

  # Shorten: new start/end
  xs <- df[[start_x]] + start_shrink * dx
  ys <- df[[start_y]] + start_shrink * dy
  xe <- df[[end_x]]   - end_shrink   * dx
  ye <- df[[end_y]]   - end_shrink   * dy

  out <- df
  out$x_s <- xs
  out$y_s <- ys
  out$x_e <- xe
  out$y_e <- ye
  out
}

# Optional: slightly vary the end shrink based on target node type
# FG circles often larger; give them a bit more clearance than pyov nodes.
edge_coords$end_shrink <- ifelse(edge_coords$type_to == "fg", 0.12, 0.09)

shorten_segments_vec <- function(df,
                                 start_x = "x", start_y = "y",
                                 end_x = "x_to", end_y = "y_to",
                                 start_shrink = 0.03,
                                 end_shrink_col = "end_shrink") {
  dx <- df[[end_x]] - df[[start_x]]
  dy <- df[[end_y]] - df[[start_y]]
  d  <- sqrt(dx * dx + dy * dy)
  d[d == 0] <- 1e-9
  es <- if (!is.null(end_shrink_col) && end_shrink_col %in% names(df)) df[[end_shrink_col]] else 0.10

  xs <- df[[start_x]] + start_shrink * dx
  ys <- df[[start_y]] + start_shrink * dy
  xe <- df[[end_x]]   - es           * dx
  ye <- df[[end_y]]   - es           * dy

  out <- df
  out$x_s <- xs
  out$y_s <- ys
  out$x_e <- xe
  out$y_e <- ye
  out
}

# Prepare shortened edge coordinates (production/utilization)
edge_prod <- edge_coords[edge_coords$edge_type == "production", , drop = FALSE]
edge_util <- edge_coords[edge_coords$edge_type == "utilization", , drop = FALSE]

edge_prod <- shorten_segments_vec(edge_prod, start_shrink = 0.03, end_shrink_col = "end_shrink")
edge_util <- shorten_segments_vec(edge_util, start_shrink = 0.03, end_shrink_col = "end_shrink")

# -----------------------------------------------------------------------------
# Legend (inset) data
# -----------------------------------------------------------------------------
x_rng <- range(node_df$x, na.rm = TRUE)
y_rng <- range(node_df$y, na.rm = TRUE)
x_span <- x_rng[2] - x_rng[1]
y_span <- y_rng[2] - y_rng[1]

# Place legend near top-left
x_legend <- x_rng[1] + 0.01 * x_span
y_top    <- y_rng[2] - 0.01 * y_span
vsep     <- 0.05 * y_span

legend_labels <- c(
  "FG: Single-receptor producer",
  "FG: Multi-producer",
  "FG: Nonproducer",
  "Pyoverdine (lock-key) group",
  "Other / unclassified (gray)"
)
legend_shapes <- c(21, 21, 21, 23, 21)
legend_fills  <- c("green", "yellow", "red", "orange", "gray80")

legend_df <- data.frame(
  lx    = rep(x_legend, length(legend_labels)),
  ly    = y_top - (0:(length(legend_labels)-1)) * vsep,
  label = factor(legend_labels, levels = legend_labels),
  shape = legend_shapes,
  fill  = legend_fills,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Plot: draw nodes first, then shortened edges on TOP so arrowheads are visible
# -----------------------------------------------------------------------------
p <- ggplot() +
  # Production edges (shortened) - drawn first (beneath nodes)
  geom_segment(
    data = edge_prod,
    aes(x = x_s, y = y_s, xend = x_e, yend = y_e),
    colour = "gray30", alpha = 0.45, linewidth = 0.25,
    arrow = arrow(length = unit(1.8, "mm"), type = "closed", angle = 15),
    lineend = "round"
  ) +
  # Utilization edges (shortened)
  geom_segment(
    data = edge_util,
    aes(x = x_s, y = y_s, xend = x_e, yend = y_e),
    colour = "gray50", alpha = 0.45, linewidth = 0.25,
    arrow = arrow(length = unit(1.8, "mm"), type = "closed", angle = 15),
    lineend = "round"
  ) +
  # Functional group nodes (filled circles) drawn on top of edges
  geom_point(
    data = node_df[node_df$type == "fg", , drop = FALSE],
    aes(x = x, y = y, size = size),
    shape = 21, stroke = 0.2, colour = "black",
    fill = node_df$color[node_df$type == "fg"],
    alpha = 0.9
  ) +
  # Pyoverdine nodes (filled diamond-like markers) drawn on top of edges
  geom_point(
    data = node_df[node_df$type == "pyov", , drop = FALSE],
    aes(x = x, y = y, size = size),
    shape = 23, stroke = 0.2, colour = "black",
    fill = node_df$color[node_df$type == "pyov"],
    alpha = 0.95
  ) +
  scale_size_continuous(range = c(2, 12), guide = guide_legend(title = "Node Size")) +
  labs(
    title    = "Pyoverdine-Mediated Iron Interaction Network",
    subtitle = "Bipartite: Functional Groups ↔ Pyoverdine Groups\nProduction: FG → PYO  |  Utilization: PYO → FG",
    caption  = "Node sizes ∝ strain counts / usage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    plot.caption  = element_text(hjust = 0.5, size = 9),
    panel.grid    = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  # Inset legend (top-left)
  geom_point(
    data = legend_df,
    aes(x = lx, y = ly),
    shape = legend_df$shape, fill = legend_df$fill,
    colour = "black", size = 4, stroke = 0.2
  ) +
  geom_text(
    data = legend_df,
    aes(x = lx + 0.02 * x_span, y = ly, label = label),
    hjust = 0, vjust = 0.5, size = 3.2
  ) +
  # Gray explanatory text (bottom-right)

  coord_fixed()

# -----------------------------------------------------------------------------
# Save outputs and metrics
# -----------------------------------------------------------------------------
pdf_path <- file.path(img_dir, "fg_pyov_network_panel3A.pdf")
png_path <- file.path(img_dir, "fg_pyov_network_panel3A.png")

ggsave(pdf_path, p, device = cairo_pdf, width = 10, height = 8, dpi = 300)
ggsave(png_path, p, device = "png",       width = 10, height = 8, dpi = 300)

# Metrics
metrics <- data.frame(
  vertices   = vcount(g),
  edges      = ecount(g),
  avg_degree = mean(degree(g)),
  density    = edge_density(g),
  components = components(g)$no
)

write.csv(metrics, file.path(docs_dir, "visualization_metrics.csv"), row.names = FALSE)

writeLines(c(
  "Phase 3 Network Visualization Summary",
  paste("Run timestamp:", Sys.time()),
  "",
  "Network metrics:",
  paste(sprintf("  Vertices: %d", metrics$vertices)),
  paste(sprintf("  Edges: %d", metrics$edges)),
  paste(sprintf("  Density: %.4f", metrics$density)),
  "",
  "Files created:",
  paste("  -", pdf_path),
  paste("  -", png_path),
  paste("  -", file.path(docs_dir, "visualization_metrics.csv"))
), file.path(docs_dir, "visualization_summary.txt"))

cat("\n=== Phase 3 Visualization Step 4 Complete ===\n")
cat("Created files:\n")
cat("  ", pdf_path, "\n")
cat("  ", png_path, "\n")
cat("  ", file.path(docs_dir, "visualization_metrics.csv"), "\n")
cat("  ", file.path(docs_dir, "visualization_summary.txt"), "\n")
