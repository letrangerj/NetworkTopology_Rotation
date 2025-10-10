#!/usr/bin/env Rscript

# Phase 3 - Step 4: Network Visualization (Fig 3A style)
# Updated: include all Phase 2 strategy classes in color mapping and legend
# Dependencies: igraph, ggplot2
#
# Generates:
#   - figures/network/fg_pyov_network_panel3A.pdf
#   - figures/network/fg_pyov_network_panel3A.png
#   - docs/phase_03/visualization_metrics.csv
#   - docs/phase_03/visualization_summary.txt
#
# Note: this script expects the following CSVs produced in earlier steps:
#   - data/interim/edges_functional_groups_conservative.csv
#   - data/interim/nodes_pyoverdines_conservative.csv
#   - data/interim/nodes_functional_groups_conservative.csv

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
})

set.seed(2025)

# Paths
img_dir <- "figures/network"
docs_dir <- "docs/phase_03"
data_edges <- "data/interim/edges_functional_groups_conservative.csv"
data_pyov <- "data/interim/nodes_pyoverdines_conservative.csv"
data_fg <- "data/interim/nodes_functional_groups_conservative.csv"

if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(docs_dir)) dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

# Read inputs
edges <- read.csv(data_edges, stringsAsFactors = FALSE)
pyov <- read.csv(data_pyov, stringsAsFactors = FALSE)
fg <- read.csv(data_fg, stringsAsFactors = FALSE)

# Build vertex tables
vertex_pyov <- data.frame(
  name = pyov$node_id,
  type = "pyov",
  color = "orange", # default pyov color
  size = pmax(2, sqrt(pyov$matrix_index) * 1.4),
  stringsAsFactors = FALSE
)

# Map all Phase 2 strategy classes explicitly:
# - Single-receptor producer (green)
# - Multi-producer (orange)
# - Multi-receptor producer (blue)
# - Nonproducer (red)
# - Producer-only (no utilization) (purple)
# Any other label -> fallback gray80
vertex_fg <- data.frame(
  name = fg$agent_id,
  type = "fg",
  color = ifelse(
    fg$strategy == "Single-receptor producer", "green",
    ifelse(
      fg$strategy == "Multi-producer", "magenta",
      ifelse(
        fg$strategy == "Multi-receptor producer", "blue",
        ifelse(
          fg$strategy == "Producer-only (no utilization)", "purple",
          ifelse(fg$strategy == "Nonproducer", "red", "gray80")
        )
      )
    )
  ),
  size = pmax(3, sqrt(fg$strain_count) * 2),
  stringsAsFactors = FALSE
)

vertex_df <- rbind(vertex_pyov, vertex_fg)

# Build edge data frame for igraph
edge_data <- data.frame(
  from = edges$source,
  to = edges$target,
  edge_type = edges$edge_type,
  stringsAsFactors = FALSE
)

# Build igraph
g <- graph_from_data_frame(d = edge_data, directed = TRUE, vertices = vertex_df)

# Compute layout
set.seed(42)
coords <- layout_with_fr(g, niter = 500)

node_df <- data.frame(
  name = V(g)$name,
  x = coords[, 1],
  y = coords[, 2],
  type = as.character(V(g)$type),
  color = as.character(V(g)$color),
  size = as.numeric(V(g)$size),
  stringsAsFactors = FALSE
)

# Merge edge coordinates with node coords
edge_coords <- merge(edge_data, node_df, by.x = "from", by.y = "name")
edge_coords <- merge(edge_coords, node_df, by.x = "to", by.y = "name", suffixes = c("", "_to"))
# columns: from,to,edge_type, x,y,type,color,size, x_to,y_to,type_to,color_to,size_to

# Shorten segments so arrowheads land outside target nodes
shorten_segments_vec <- function(df, start_shrink = 0.03, end_shrink_default = 0.10) {
  dx <- df$x_to - df$x
  dy <- df$y_to - df$y
  d <- sqrt(dx * dx + dy * dy)
  d[d == 0] <- 1e-9
  # choose end_shrink per-target-type (FG targets may need more clearance)
  end_shrink <- ifelse(!is.null(df$type_to) & df$type_to == "fg", 0.12, end_shrink_default)
  xs <- df$x + start_shrink * dx
  ys <- df$y + start_shrink * dy
  xe <- df$x_to - end_shrink * dx
  ye <- df$y_to - end_shrink * dy
  df$x_s <- xs
  df$y_s <- ys
  df$x_e <- xe
  df$y_e <- ye
  df
}

edge_prod <- subset(edge_coords, edge_type == "production")
edge_util <- subset(edge_coords, edge_type == "utilization")

if (nrow(edge_prod) > 0) edge_prod <- shorten_segments_vec(edge_prod, start_shrink = 0.03, end_shrink_default = 0.09)
if (nrow(edge_util) > 0) edge_util <- shorten_segments_vec(edge_util, start_shrink = 0.03, end_shrink_default = 0.09)

# Legend setup (top-left)
x_rng <- range(node_df$x, na.rm = TRUE)
y_rng <- range(node_df$y, na.rm = TRUE)
x_span <- x_rng[2] - x_rng[1]
y_span <- y_rng[2] - y_rng[1]

x_legend <- x_rng[1] + 0.01 * x_span
y_top <- y_rng[2] - 0.01 * y_span
vsep <- 0.05 * y_span

legend_labels <- c(
  "FG: Single-receptor producer",
  "FG: Multi-producer",
  "FG: Multi-receptor producer",
  "FG: Producer-only (no utilization)",
  "FG: Nonproducer",
  "Pyoverdine (lock-key) group"
)
legend_shapes <- c(21, 21, 21, 21, 21, 23)
legend_fills <- c("green", "magenta", "blue", "purple", "red", "orange")

legend_df <- data.frame(
  lx = rep(x_legend, length(legend_labels)),
  ly = y_top - (0:(length(legend_labels) - 1)) * vsep,
  label = factor(legend_labels, levels = legend_labels),
  shape = legend_shapes,
  fill = legend_fills,
  stringsAsFactors = FALSE
)

# Build plot: draw edges first, then nodes (so nodes sit cleanly above edges)
p <- ggplot() +
  # edges drawn FIRST (below nodes) to avoid overlap ambiguity
  {
    if (nrow(edge_prod) > 0) {
      geom_segment(
        data = edge_prod,
        aes(x = x_s, y = y_s, xend = x_e, yend = y_e),
        colour = "gray30", alpha = 0.35, linewidth = 0.25,
        arrow = arrow(length = unit(1.8, "mm"), type = "closed", angle = 15),
        lineend = "round"
      )
    } else {
      NULL
    }
  } +
  {
    if (nrow(edge_util) > 0) {
      geom_segment(
        data = edge_util,
        aes(x = x_s, y = y_s, xend = x_e, yend = y_e),
        colour = "gray50", alpha = 0.35, linewidth = 0.25,
        arrow = arrow(length = unit(1.8, "mm"), type = "closed", angle = 15),
        lineend = "round"
      )
    } else {
      NULL
    }
  } +
  # nodes (FG) drawn AFTER edges
  geom_point(
    data = node_df[node_df$type == "fg", , drop = FALSE],
    aes(x = x, y = y, size = size),
    shape = 21, stroke = 0.2, colour = "black",
    fill = node_df$color[node_df$type == "fg"],
    alpha = 0.9
  ) +
  # nodes (pyov) drawn AFTER edges
  geom_point(
    data = node_df[node_df$type == "pyov", , drop = FALSE],
    aes(x = x, y = y, size = size),
    shape = 23, stroke = 0.2, colour = "black",
    fill = node_df$color[node_df$type == "pyov"],
    alpha = 0.95
  ) +
  scale_size_continuous(range = c(2, 12), guide = guide_legend(title = "Node Size")) +
  labs(
    title = "Pyoverdine-Mediated Iron Interaction Network",
    subtitle = "Bipartite: Functional Groups ↔ Pyoverdine Groups\nEdge direction: FG → PYO (production), PYO → FG (utilization)",
    caption = "Node sizes ∝ strain counts / usage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    plot.caption = element_text(hjust = 0.5, size = 9),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  # inset legend (top-left)
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
  coord_fixed()

# Save outputs
pdf_path <- file.path(img_dir, "fg_pyov_network_panel3A.pdf")
png_path <- file.path(img_dir, "fg_pyov_network_panel3A.png")

ggsave(pdf_path, p, device = cairo_pdf, width = 10, height = 8, dpi = 300)
ggsave(png_path, p, device = "png", width = 10, height = 8, dpi = 300)

# Metrics and provenance
metrics <- data.frame(
  vertices = vcount(g),
  edges = ecount(g),
  avg_degree = mean(degree(g)),
  density = edge_density(g),
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
  "Visualization improvements:",
  "  - Edges drawn first, nodes on top to reduce overlap ambiguity",
  "  - Updated edge direction explanation in subtitle",
  "  - Removed gray node annotation",
  "",
  "Files created:",
  paste("  -", pdf_path),
  paste("  -", png_path),
  paste("  -", file.path(docs_dir, "visualization_metrics.csv"))
), file.path(docs_dir, "visualization_summary.txt"))

cat("\n=== Phase 3 Visualization Step 4 Complete ===\n")
cat("Created files:\n")
cat(" ", pdf_path, "\n")
cat(" ", png_path, "\n")
cat(" ", file.path(docs_dir, "visualization_metrics.csv"), "\n")
cat(" ", file.path(docs_dir, "visualization_summary.txt"), "\n")
