# Phase 1 (alternative to 02_download_osm.R) - Chennai roads from a Geofabrik extract.
#
# Why: the public Overpass servers were overloaded (HTTP 504 / connection resets)
# while building this project, so the road network is taken from the regional
# Geofabrik extract of southern India (ONE download, ~560 MB, kept OUTSIDE the
# OneDrive folder in C:/cityshield-raw and not needed again afterwards).
#
#   1. curl -L -C - -o C:/cityshield-raw/southern-zone.osm.pbf \
#        https://download.geofabrik.de/asia/india/southern-zone-latest.osm.pbf
#   2. Rscript scripts/02b_extract_osm_geofabrik.R
#
# Output: data/osm/chennai_roads_raw.rds - the SAME format 02_download_osm.R
# writes (sf lines with osm_id, highway, oneway, junction), so
# 03_build_graph.R works with either source. GDAL's OSM driver reads the .pbf;
# this is data loading only (no graded algorithm involved).

source("R/core/load_all.R")
suppressPackageStartupMessages(library(sf))

pbf <- Sys.getenv("CITYSHIELD_PBF", "C:/cityshield-raw/southern-zone.osm.pbf")
out_final <- cs_path("data", "osm", "chennai_roads_raw.rds")
if (file.exists(out_final)) { cat("Already cached:", out_final, "\n"); quit(save = "no") }
stopifnot(file.exists(pbf))
dir.create(dirname(out_final), recursive = TRUE, showWarnings = FALSE)

# Greater Chennai core (lon_min, lat_min, lon_max, lat_max)
bb <- c(80.12, 12.90, 80.31, 13.19)
wkt <- sprintf("POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",
               bb[1], bb[2], bb[3], bb[2], bb[3], bb[4], bb[1], bb[4], bb[1], bb[2])
hw <- c("motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential",
        "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link")
q <- sprintf("SELECT osm_id, highway, other_tags FROM lines WHERE highway IN (%s)",
             paste0("'", hw, "'", collapse = ","))
Sys.setenv(OSM_MAX_TMPFILE_SIZE = "2000")        # MB of temp storage GDAL may use
cat("reading", pbf, "- this takes several minutes ...\n"); flush.console()
t0 <- proc.time()[["elapsed"]]
roads <- sf::st_read(pbf, layer = "lines", query = q, wkt_filter = wkt, quiet = TRUE)
cat(sprintf("read %d ways in %.0f s\n", nrow(roads), proc.time()[["elapsed"]] - t0))

# other_tags is hstore text: "oneway"=>"yes","junction"=>"roundabout"
tag <- function(x, key) {
  m <- regmatches(x, regexpr(sprintf('"%s"=>"[^"]*"', key), x))
  out <- rep(NA_character_, length(x)); has <- grepl(sprintf('"%s"=>"', key), x, fixed = TRUE)
  out[has] <- sub(sprintf('^"%s"=>"([^"]*)"$', key), "\\1", m)
  out
}
roads$oneway <- tag(roads$other_tags, "oneway")
roads$junction <- tag(roads$other_tags, "junction")
roads <- roads[, c("osm_id", "highway", "oneway", "junction")]
roads <- roads[!duplicated(roads$osm_id), ]
cat("ways kept:", nrow(roads), "\n")
saveRDS(roads, out_final)
cat("Saved", out_final, "\n")
