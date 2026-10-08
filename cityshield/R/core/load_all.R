# =============================================================================
# core/load_all.R  -  THE central load script
#
# Usage (from anywhere):   source("<path to>/cityshield/R/core/load_all.R")
# or, from inside cityshield/:  source("R/core/load_all.R")
#
# It
#   1. finds the project root (the folder that contains R/ and data/),
#   2. defines PROJECT_ROOT and path helpers,
#   3. sources core files first (order matters: structures before graph),
#   4. then every file in each module folder m1_ ... m7_ (alphabetical order).
# Module files must therefore only define functions (no side effects).
# =============================================================================

# ---- 1. locate project root --------------------------------------------------
.find_root <- function() {
  # Prefer an explicit override, then walk up from the working directory.
  env <- Sys.getenv("CITYSHIELD_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/"))
  d <- normalizePath(getwd(), winslash = "/")
  for (i in 1:6) {
    if (file.exists(file.path(d, "R", "core", "load_all.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  stop("Cannot find CityShield root. Set CITYSHIELD_ROOT or setwd() into the project.")
}

PROJECT_ROOT <- .find_root()

#' Path helper: cs_path("data", "synthetic", "incidents.rds")
cs_path <- function(...) file.path(PROJECT_ROOT, ...)

# ---- 2. source core (explicit order) ----------------------------------------
for (f in c("structures.R", "sort_utils.R", "graph.R", "geometry_types.R", "timer.R", "study_area.R")) {
  source(cs_path("R", "core", f), local = FALSE)
}
source(cs_path("R", "core", "data_io.R"), local = FALSE)

# ---- 3. source every module folder ------------------------------------------
for (mod in c("m1_greedy_dc", "m2_dp_bt_bb", "m3_strings", "m4_graphs",
              "m5_geometry", "m6_randomized", "m7_complexity")) {
  files <- list.files(cs_path("R", mod), pattern = "\\.[Rr]$", full.names = TRUE)
  for (f in files) source(f, local = FALSE)       # list.files() is sorted
}
invisible(NULL)
