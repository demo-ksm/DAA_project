# Phase 1 - Turn the cached OSM ways into the CityShield CSR graph.
#
#   Rscript scripts/03_build_graph.R
#
# Pipeline (all preprocessing; vectorised base R + sf/igraph are allowed here
# because this is DATA LOADING, not a graded algorithm):
#   1. read every way's vertex list
#   2. give each distinct coordinate an id (rounded to ~0.1 m)
#   3. a vertex is a JUNCTION if it is used >= 2 times or ends a way
#   4. cut ways at junctions -> one graph edge per stretch of road
#      (this "contraction" removes the thousands of degree-2 shape points)
#   5. direction from the `oneway` tag, travel time from road class
#   6. keep the largest strongly connected component (so every node can reach
#      every other: routing queries never hit dead ends)
#   7. build the CSR graph with our own csr_build() and saveRDS()

source("R/core/load_all.R")
suppressPackageStartupMessages({ library(sf); library(igraph) })

roads <- readRDS(cs_path("data", "osm", "chennai_roads_raw.rds"))
roads <- roads[!sf::st_is_empty(roads), ]
cat("ways:", nrow(roads), "\n")

# ---- 1-2. vertex table ------------------------------------------------------
xy <- sf::st_coordinates(roads)              # columns X, Y, L1 (= way index)
way <- as.integer(xy[, "L1"])
lon <- xy[, "X"]; lat <- xy[, "Y"]
R <- nrow(xy)
cat("vertex rows:", R, "\n")

key <- round((lon - 80) * 1e6) * 4e6 + round((lat - 12) * 1e6)   # exact in double
vid <- match(key, unique(key))               # vertex id of each row
nv  <- max(vid)

first_row <- c(TRUE, way[-1] != way[-R])     # first vertex of its way
last_row  <- c(way[-1] != way[-R], TRUE)     # last vertex of its way

# ---- 3. junctions -----------------------------------------------------------
use_count <- tabulate(vid, nbins = nv)
is_junction <- use_count >= 2L
is_junction[vid[first_row | last_row]] <- TRUE
brk <- is_junction[vid]                      # row is a cut point

# ---- 4. cut ways -> edges ---------------------------------------------------
haversine <- function(lon1, lat1, lon2, lat2) {
  p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
       cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * 6371000 * asin(sqrt(a))
}
seg <- numeric(R)                            # length of segment ending at row r
seg[-1] <- haversine(lon[-R], lat[-R], lon[-1], lat[-1])
seg[first_row] <- 0
cum <- cumsum(seg)

idx <- seq_len(R)
last_break <- cummax(ifelse(brk, idx, 0L))   # last break row <= r
prev_break <- c(0L, last_break[-R])          # last break row <  r
is_end <- brk & !first_row                   # rows that terminate an edge
r_end   <- idx[is_end]
r_start <- prev_break[is_end]
e_from <- vid[r_start]; e_to <- vid[r_end]
e_len  <- cum[r_end] - cum[r_start]
e_way  <- way[r_end]
keep <- e_from != e_to & e_len > 0           # drop self-loops / zero length
e_from <- e_from[keep]; e_to <- e_to[keep]; e_len <- e_len[keep]; e_way <- e_way[keep]

# ---- 5. direction, speed, capacity -----------------------------------------
hw <- sub("_link$", "", roads$highway[e_way])
code <- match(hw, ROAD_CLASSES$highway)
code[is.na(code)] <- 7L                      # unknown -> residential
ow <- if ("oneway" %in% names(roads)) roads$oneway[e_way] else rep(NA_character_, length(e_way))
jn <- if ("junction" %in% names(roads)) roads$junction[e_way] else rep(NA_character_, length(e_way))
fwd_only <- (!is.na(ow) & ow %in% c("yes", "1", "true")) | (!is.na(jn) & jn == "roundabout")
bwd_only <- !is.na(ow) & ow == "-1"
both <- !(fwd_only | bwd_only)

# forward copies (all but backward-only) + reverse copies (both / backward-only)
f_idx <- which(!bwd_only)
b_idx <- which(both | bwd_only)
from <- c(e_from[f_idx], e_to[b_idx])
to   <- c(e_to[f_idx],   e_from[b_idx])
len  <- c(e_len[f_idx],  e_len[b_idx])
cls  <- c(code[f_idx],   code[b_idx])
cat("directed edges before cleaning:", length(from), " nodes:", nv, "\n")

# ---- 6. largest strongly connected component --------------------------------
gi <- igraph::make_graph(rbind(from, to), n = nv, directed = TRUE)
comp <- igraph::components(gi, mode = "strong")
big <- which.max(comp$csize)
in_big <- comp$membership == big
new_id <- cumsum(in_big)                     # 1..n for kept nodes
keepe <- in_big[from] & in_big[to]
from <- new_id[from[keepe]]; to <- new_id[to[keepe]]
len <- len[keepe]; cls <- cls[keepe]
n <- sum(in_big)

# node coordinates (first row where each vertex appears)
first_of <- match(seq_len(nv), vid)
nx <- lon[first_of][in_big]; ny <- lat[first_of][in_big]

# ---- 7. build our own CSR graph --------------------------------------------
speed_ms <- ROAD_CLASSES$speed_kmh[cls] / 3.6
tt  <- len / speed_ms                        # travel time in seconds
cap <- ROAD_CLASSES$cap_vph[cls]             # vehicles / hour
g <- csr_build(n, from, to, tt, x = nx, y = ny,
               extra = list(len_m = len, cap = cap, hw = cls))
g$meta <- list(source = "OpenStreetMap via osmdata", built = Sys.time(),
               bbox = c(80.12, 12.90, 80.31, 13.19),
               note = "largest SCC; degree-2 vertices contracted")
saveRDS(g, cs_path("data", "osm", "chennai_graph.rds"))
cat(sprintf("Saved graph: %d nodes, %d directed edges\n", g$n, g$m))
