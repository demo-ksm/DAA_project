# Phase 1 - Download Chennai's road network from OpenStreetMap, ONCE.
#
#   Rscript scripts/02_download_osm.R        (safe to re-run: it resumes)
#
# The merged raw result is cached to data/osm/chennai_roads_raw.rds. If that
# file already exists this script does nothing, so the download never repeats.
# Delete the file manually if you really want a fresh download.
#
# Overpass often times out on dense areas, so we download a grid of tiles and
# cache each tile on disk. A tile that keeps failing is SPLIT into four smaller
# tiles (down to MIN_SIZE degrees) - sparse areas stay cheap, dense ones get
# small queries. Servers are rotated when a connection is reset.

source("R/core/load_all.R")   # PROJECT_ROOT / cs_path (module folders may be empty)
suppressPackageStartupMessages({ library(osmdata); library(sf) })

out_final <- cs_path("data", "osm", "chennai_roads_raw.rds")
tile_dir  <- cs_path("data", "osm", "tiles")
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(out_final)) {
  cat("Raw OSM cache already exists - nothing to download:\n  ", out_final, "\n")
  quit(save = "no")
}

# Greater Chennai city core. (lon_min, lat_min, lon_max, lat_max)
bb <- c(80.12, 12.90, 80.31, 13.19)
nx <- 6; ny <- 8
xs <- seq(bb[1], bb[3], length.out = nx + 1)
ys <- seq(bb[2], bb[4], length.out = ny + 1)
MIN_SIZE <- 0.012                                   # ~1.3 km: do not split below this

hw_values <- c("motorway", "trunk", "primary", "secondary", "tertiary",
               "unclassified", "residential",
               "motorway_link", "trunk_link", "primary_link",
               "secondary_link", "tertiary_link")

servers <- c("https://overpass-api.de/api/interpreter",
             "https://lz4.overpass-api.de/api/interpreter",
             "https://z.overpass-api.de/api/interpreter",
             "https://maps.mail.ru/osm/tools/overpass/api/interpreter")

# one attempt-series for a box; returns an sf object, or NULL after `tries` failures
fetch_box <- function(box, tries = 3) {
  for (k in seq_len(tries)) {
    osmdata::set_overpass_url(servers[(k - 1L) %% length(servers) + 1L])
    res <- try({
      q <- opq(bbox = box, timeout = 120) |> add_osm_feature(key = "highway", value = hw_values)
      lines <- osmdata_sf(q)$osm_lines
      if (is.null(lines)) lines <- sf::st_sf(osm_id = character(0), highway = character(0),
                                              geometry = sf::st_sfc(crs = 4326))
      keep <- intersect(c("osm_id", "highway", "oneway", "junction"), names(lines))
      lines[, keep]
    }, silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    cat("   retry", k, "-", substr(conditionMessage(attr(res, "condition")), 1, 60), "\n")
    Sys.sleep(min(45, 8 * k))
  }
  NULL
}

box_key <- function(b) sprintf("box_%.4f_%.4f_%.4f_%.4f.rds", b[1], b[2], b[3], b[4])

# work list of boxes; initial grid first (legacy file names kept for resuming)
work <- list()
for (i in rev(seq_len(nx))) for (j in rev(seq_len(ny)))
  work[[length(work) + 1L]] <- list(box = c(xs[i], ys[j], xs[i + 1], ys[j + 1]),
                                    legacy = sprintf("tile_%d_%d.rds", i, j))
done_files <- character(0)
while (length(work) > 0L) {
  w <- work[[length(work)]]; work[[length(work)]] <- NULL
  f <- file.path(tile_dir, box_key(w$box))
  fl <- if (!is.null(w$legacy)) file.path(tile_dir, w$legacy) else ""
  if (file.exists(fl)) { done_files <- c(done_files, fl); next }
  if (file.exists(f))  { done_files <- c(done_files, f);  next }
  cat(sprintf("box %.3f %.3f %.3f %.3f\n", w$box[1], w$box[2], w$box[3], w$box[4]))
  t <- fetch_box(w$box)
  if (!is.null(t)) {
    saveRDS(t, f); done_files <- c(done_files, f); Sys.sleep(2)
  } else if (w$box[3] - w$box[1] > MIN_SIZE) {      # split into 4 and retry smaller
    cat("   splitting\n")
    mx <- (w$box[1] + w$box[3]) / 2; my <- (w$box[2] + w$box[4]) / 2
    for (q in list(c(w$box[1], w$box[2], mx, my), c(mx, w$box[2], w$box[3], my),
                   c(w$box[1], my, mx, w$box[4]), c(mx, my, w$box[3], w$box[4])))
      work[[length(work) + 1L]] <- list(box = q)
  } else stop("box failed even at minimum size: ", paste(w$box, collapse = " "))
}

tiles <- list()
for (f in done_files) tiles[[length(tiles) + 1L]] <- readRDS(f)
roads <- do.call(rbind, tiles)
roads <- roads[!duplicated(roads$osm_id), ]       # ways returned by two neighbouring tiles
cat("ways:", nrow(roads), "\n")
saveRDS(roads, out_final)
cat("Saved", out_final, "\n")
