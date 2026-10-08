# Phase 1 - one-off environment setup.
#
# Step 1 (once):   Rscript -e "renv::init(bare = TRUE, restart = FALSE)"
#                  -> creates .Rprofile + renv/ (activate.R). Do it from cityshield/.
# Step 2 (this):   Rscript scripts/01_setup_renv.R
#                  -> installs every package and writes renv.lock.
#
# Design note: .Renviron points renv's cache and library OUTSIDE OneDrive
# (C:/renv-cache, C:/renv-lib). Only renv.lock (a small text file) is stored in
# the project, so OneDrive never syncs thousands of package files.
#
# (We split init and install because renv::init() inside a single Rscript call
#  tries to restart the session and hangs on Windows.)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(pkgType = "binary")
options(timeout = 3600)               # default 60 s truncates big binaries (terra, sf) -> corrupt zips
           # Windows binaries: no Rtools needed

pkgs <- c(
  # data loading / geometry (allowed ONLY for loading + validating in tests)
  "osmdata", "sf", "igraph",
  # testing + benchmarking
  "testthat", "bench", "ggplot2",
  # data structures
  "R6",
  # dashboard
  "shiny", "bslib", "leaflet", "plotly",
  # report
  "rmarkdown", "knitr", "jsonlite", "dplyr"
)

# Slow/flaky networks sometimes leave a truncated zip in the download cache and
# renv rolls the whole install back. Retrying re-downloads only what is missing.
for (attempt in 1:4) {
  ok <- tryCatch({ renv::install(pkgs, prompt = FALSE); TRUE },
                 error = function(e) { message("attempt ", attempt, " failed: ", conditionMessage(e)); FALSE })
  if (ok) break
  zips <- list.files(file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "renv", "binary"),
                     pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
  for (z in zips) if (inherits(try(utils::unzip(z, list = TRUE), silent = TRUE), "try-error")) unlink(z)
}
if (!ok) stop("package installation failed after retries")
renv::snapshot(type = "all", prompt = FALSE)
cat("renv setup complete.\n")
