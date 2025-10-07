#!/usr/bin/env Rscript

# Phase 3 - Step 4: Fixed Network Visualization (Fig 3A style)
# Base R compatible version with proper edge attribute handling

library(igraph)
library(ggplot2)

set.seed(2025)
img_dir   <- "figures/network"
docs_dir  <- "docs/phase_03"

dir.create(img_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)

# Load data -------------------------------------------------------------------
edges  <- read.csv("data/interim/edges_functional_groups_conservative.csv", stringsAsFactors = FALSE)
pyov   <- read.csv("data/interim/nodes_pyoverdines_conservative.csv",       stringsAsFactors = FALSE)
fg     <- read.csv("data/interim/nodes_functional_groups_conservative.csv", stringsAsFactors = FALSE)

# Build vertex table ----------------------------------------------------------
# Pyov nodes
vertex_pyov <- data.frame(
  name  = pyov$node_id,
  type  = "pyov",
  color = "orange",
  size  = pmax(2, sqrt(pyov$matrix_index) * 1.4),
  stringsAsFactors = FALSE
)

# FG nodes with strategy-based coloring
vertex_fg <- data.frame(
  name  = fg$agent_id,
  type  = "fg",
  color = ifelse(fg$strategy == "Multi-producer",        "yellow",
          ifelse(fg$strategy == "Single-receptor producer", "green",
          ifelse(fg$strategy == "Nonproducer",           "red", "gray80"))),
  size  = pmax(3, sqrt(fg$strain_count) * 2),
  stringsAsFactors = FALSE
)

vertex_df <- rbind(vertex_pyov, vertex_fg)
vertex_df$id <- seq_len(nrow(vertex_df))

# Build igraph with proper edge attributes ------------------------------------
edge_data <- data.frame(
  from = edges$source,
  to   = edges$target,
  edge_type = edges$edge_type,
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(
  d = edge_data,
  directed = TRUE,
  vertices = vertex_df
)

# Layout: Fruchterman-Reingold ------------------------------------------------
set.seed(42)
layout_mat <- layout_with_fr(g, niter = 500)

# Prepare node data for ggplot ------------------------------------------------
node_df <- data.frame(
  name  = V(g)$name,
  x     = layout_mat[,1],
  y     = layout_mat[,2],
  type  = V(g)$type,
  color = as.character(V(g)$color),
  size  = as.numeric(V(g)$size),
  stringsAsFactors = FALSE
)

# Prepare edge data directly from original edge data -------------------------
# Merge node coordinates with edge data
edge_coords <- merge(edge_data, node_df, by.x = "from", by.y = "name")
edge_coords <- merge(edge_coords, node_df, by.x = "to", by.y = "name", suffixes = c("", "_end"))

# Plot -----------------------------------------------------------------------
p <- ggplot() +
  # Production edges (FG -> PYO)
  geom_segment(
    data = edge_coords[edge_coords$edge_type == "production", ],
    aes(x = x, y = y, xend = x_end, yend = y_end),
    colour = "gray30", alpha = 0.4, size = 0.25,
    arrow = arrow(length = unit(0.4, "mm"), type = "closed", angle = 15),
    lineend = "round"
  ) +
  # Utilization edges (PYO -> FG)
  geom_segment(
    data = edge_coords[edge_coords$edge_type == "utilization", ],
    aes(x = x, y = y, xend = x_end, yend = y_end),
    colour = "gray50", alpha = 0.4, size = 0.25,
    arrow = arrow(length = unit(0.4, "mm"), type = "closed", angle = 15),
    lineend = "round"
  ) +
  # Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, size = size, colour = color),
    shape = 19, alpha = 0.9
  ) +
  scale_size_continuous(range = c(2, 12), guide = guide_legend(title = "Node Size")) +
  scale_colour_identity(name = "Strategy/Type") +
  labs(
    title = "Pyoverdine-Mediated Iron Interaction Network",
    subtitle = "Bipartite: Functional Groups ↔ Pyoverdine Groups\nProduction: FG → PYO  |  Utilization: PYO → FG",
    caption = "Node sizes ∝ strain counts / usage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    plot.caption  = element_text(hjust = 0.5, size = 9),
    panel.grid    = element_blank(),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.box.background = element_rect(colour = "black")
  )

# Save outputs ---------------------------------------------------------------
pdf_path <- file.path(img_dir, "fg_pyov_network_panel3A.pdf")
png_path <- file.path(img_dir, "fg_pyov_network_panel3A.png")

ggsave(pdf_path, p, device = pdf, width = 10, height = 8, dpi = 300)
ggsave(png_path, p, device = png, width = 10, height = 8, dpi = 300)

# Metrics and summary ---------------------------------------------------------
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
