# =============================================================================
# benchmarks/common.R - shared helpers for every benchmark experiment
#
#   * the plot theme/palette (validated default dataviz palette, fixed order)
#   * log-log slope estimation (empirical exponent vs theory)
#   * the benchmark GRAPH: the cached Chennai graph, or (only if the cache is
#     missing / CITYSHIELD_SYNTH=1) a synthetic grid "road" graph with the same
#     fields, so the experiments can be developed before the download finishes
# All timing goes through time_it() / bench_sizes() from R/core/timer.R.
# =============================================================================

source(file.path(Sys.getenv("CITYSHIELD_ROOT", unset = getwd()), "R", "core", "load_all.R"))
suppressPackageStartupMessages(library(ggplot2))

RESULTS_DIR <- cs_path("benchmarks", "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- plotting ----------------------------------------------------------------
# Categorical slots 1-4 of the validated palette, always in this order.
SERIES <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")

bench_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
          panel.background = element_rect(fill = "#fcfcfb", colour = NA),
          panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          text = element_text(colour = "#0b0b0b"),
          axis.text = element_text(colour = "#52514e"),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "#52514e", size = 9),
          legend.position = "top", legend.title = element_blank(),
          strip.text = element_text(face = "bold"))
}

#' Save a ggplot into benchmarks/results/<name>.png (compact size).
save_plot <- function(p, name, w = 7.5, h = 4.4) {
  path <- file.path(RESULTS_DIR, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 130, bg = "#fcfcfb")
  invisible(path)
}

#' Empirical growth exponent: slope of log(y) against log(x) (least squares).
#' Used to compare observed runtime with the theoretical exponent.
loglog_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  if (sum(ok) < 3L) return(NA_real_)
  unname(stats::coef(stats::lm(log(y[ok]) ~ log(x[ok])))[2])
}

#' Runtime-vs-size line chart on log-log axes; `color` column = algorithm.
runtime_plot <- function(df, title, subtitle = NULL, facet = NULL, xlab = "input size n",
                         ycol = "median_ms", ylab = "median time (ms, log scale)") {
  p <- ggplot(df, aes(n, .data[[ycol]], colour = algorithm, group = algorithm)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = SERIES) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) + bench_theme()
  if (!is.null(facet)) p <- p + facet_wrap(stats::as.formula(paste("~", facet)), scales = "free")
  p
}

#' Write the CSV for an experiment (overwrite) and return the data frame.
write_exp <- function(df, name) { save_bench(df, name, results_dir = RESULTS_DIR, overwrite = TRUE); df }

# ---- the benchmark graph ---------------------------------------------------------

.bench_graph_cache <- new.env()
#' The graph experiments run on (cached in memory).
bench_graph <- function() {
  if (is.null(.bench_graph_cache$g)) {
    use_real <- !SYNTH_MODE && file.exists(data_path("osm", "chennai_graph.rds"))
    .bench_graph_cache$g <- if (use_real) load_chennai_graph() else grid_road_graph()
    .bench_graph_cache$source <- if (use_real) "Chennai OSM" else "synthetic grid"
    message("[bench] graph = ", .bench_graph_cache$source, " (", .bench_graph_cache$g$n, " nodes)")
  }
  .bench_graph_cache$g
}
bench_graph_source <- function() { bench_graph(); .bench_graph_cache$source }

#' BFS subgraph of k nodes from a (seeded) random root. Used for all the
#' 50-150 node instances and for the growing flow instances.
sub_graph <- function(k, seed = 1L) {
  g <- bench_graph(); set.seed(seed)
  csr_bfs_subgraph(g, sample.int(g$n, 1L), k)
}

#' Strongly connected k-node instance from the benchmark graph (core helper).
#' Returns list(g, D) with D the symmetrised metric closure.
connected_instance <- function(k, seed) study_instance(bench_graph(), k, seed)

#' Median of a numeric vector (benchmark analysis helper). O(n log n).
med <- function(x) stats::median(x, na.rm = TRUE)
