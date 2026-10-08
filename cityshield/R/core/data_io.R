# =============================================================================
# core/data_io.R  -  Loading cached datasets (no downloading here, ever)
#
# The OSM download happens exactly once, in scripts/02_download_osm.R, and is
# cached with saveRDS(). Everything else reads the cache through these helpers.
# =============================================================================

# DEVELOPMENT MODE: with CITYSHIELD_SYNTH=1 every dataset lives in data_synth/ and
# the "road network" is a synthetic grid, so the pipeline/dashboard can be tested
# without touching (or waiting for) the real OSM cache in data/.
SYNTH_MODE <- Sys.getenv("CITYSHIELD_SYNTH", "") == "1"
DATA_DIR_NAME <- if (SYNTH_MODE) "data_synth" else "data"

#' Path inside the active data directory: data_path("synthetic", "incidents.rds")
data_path <- function(...) cs_path(DATA_DIR_NAME, ...)

#' Load the cached Chennai road graph (a csr_graph + attributes).
#' Fields: n, m, offset, target, weight (travel time, s), x, y (lon/lat),
#'         extra$len_m, extra$cap (veh/h), extra$hw (road class code).
#' Stops with a helpful message if the cache has not been built yet.
load_chennai_graph <- function(path = data_path("osm", "chennai_graph.rds")) {
  if (SYNTH_MODE) return(grid_road_graph())
  if (!file.exists(path))
    stop("Graph cache not found: ", path,
         "\nRun scripts/02_download_osm.R then scripts/03_build_graph.R once.")
  readRDS(path)
}

#' Load synthetic CityShield data: list(incidents, calls, stations, hospitals).
load_synthetic <- function(dir = data_path("synthetic")) {
  files <- c(incidents = "incidents.rds", calls = "call_logs.rds",
             stations = "stations.rds", hospitals = "hospitals.rds")
  out <- list()
  for (nm in names(files)) {
    p <- file.path(dir, files[[nm]])
    if (!file.exists(p)) stop("Missing ", p, " - run scripts/04_generate_data.R")
    out[[nm]] <- readRDS(p)
  }
  out
}

# Road class codes used everywhere (index = code). Speeds are free-flow km/h and
# capacities are vehicles/hour/direction - simple planning figures, documented
# here so the report can cite them.
ROAD_CLASSES <- data.frame(
  code     = 1:8,
  highway  = c("motorway", "trunk", "primary", "secondary", "tertiary",
               "unclassified", "residential", "service"),
  speed_kmh = c(60, 50, 40, 35, 30, 25, 20, 15),
  cap_vph  = c(4000, 3600, 2400, 1800, 1200, 800, 600, 300),
  stringsAsFactors = FALSE
)

#' Nearest graph node to a (lon, lat) point. Linear scan, O(V).
#' Only used for snapping synthetic stations/incidents to the network during
#' data generation and for UI clicks, never inside a graded algorithm.
nearest_node <- function(g, lon, lat) {
  best <- 1L
  bd <- Inf
  for (i in seq_len(g$n)) {
    d <- (g$x[i] - lon)^2 + (g$y[i] - lat)^2
    if (d < bd) { bd <- d; best <- i }
  }
  best
}

#' Synthetic grid road graph (side x side), two-way edges, random travel times and
#' capacities. Same fields as the Chennai graph. DEVELOPMENT / FALLBACK ONLY.
grid_road_graph <- function(side = 60L, seed = 1L) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  id <- function(i, j) (i - 1L) * side + j
  fr <- integer(0); to <- integer(0)
  for (i in seq_len(side)) for (j in seq_len(side)) {
    if (j < side) { fr <- c(fr, id(i, j), id(i, j + 1L)); to <- c(to, id(i, j + 1L), id(i, j)) }
    if (i < side) { fr <- c(fr, id(i, j), id(i + 1L, j)); to <- c(to, id(i + 1L, j), id(i, j)) }
  }
  m <- length(fr)
  len <- runif(m / 2, 80, 400); len <- rep(len, each = 2)
  cls <- sample(3:7, m / 2, TRUE, prob = c(.1, .2, .3, .3, .1)); cls <- rep(cls, each = 2)
  x <- numeric(side * side); y <- numeric(side * side)
  for (i in seq_len(side)) for (j in seq_len(side)) { x[id(i, j)] <- 80.15 + j * 0.002; y[id(i, j)] <- 12.95 + i * 0.002 }
  g <- csr_build(side * side, fr, to, len / (ROAD_CLASSES$speed_kmh[cls] / 3.6), x = x, y = y,
                 extra = list(len_m = len, cap = ROAD_CLASSES$cap_vph[cls], hw = cls))
  g$meta <- list(source = "synthetic grid")
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  g
}
