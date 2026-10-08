# Phase 1 - Synthetic CityShield data generator (fully reproducible, seeded).
#
#   Rscript scripts/04_generate_data.R
#
# Produces in data/synthetic/ :
#   stations.rds   id, node, lon, lat, type, units
#   hospitals.rds  id, node, lon, lat, beds
#   incidents.rds  id, time_min, node, lon, lat, type, severity, demand, zone
#   call_logs.rds  data.frame(id, incident_id, text, is_duplicate, dup_of)
#                  (+ call_logs.txt: one big text for the string-matching module)
#
# Realism choices (kept simple on purpose):
#   * incidents cluster around a few "hotspot" nodes (Gaussian falloff),
#   * incident arrival rate has a day/night cycle plus ONE storm surge window,
#     so "worst continuous incident-load period" (Max Subarray) has a real answer,
#   * ~15% of call logs are noisy re-reports of an earlier incident, so the
#     LCS de-duplication module has ground truth to be scored against.

source("R/core/load_all.R")
g <- load_chennai_graph()
set.seed(2026)                                   # reproducibility
out_dir <- data_path("synthetic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----------------------------------------------------------------
pick_nodes <- function(k, weights = NULL) sample.int(g$n, k, replace = FALSE, prob = weights)

# ---- hospitals and stations -------------------------------------------------
n_hosp <- 20L; n_stat <- 30L
hn <- pick_nodes(n_hosp)
hospitals <- data.frame(id = seq_len(n_hosp), node = hn,
                        lon = g$x[hn], lat = g$y[hn],
                        beds = sample(c(50, 100, 200, 400, 800), n_hosp, replace = TRUE,
                                      prob = c(.1, .25, .3, .25, .1)))
sn <- pick_nodes(n_stat)
stations <- data.frame(id = seq_len(n_stat), node = sn,
                       lon = g$x[sn], lat = g$y[sn],
                       type = sample(c("fire", "police", "ambulance"), n_stat, replace = TRUE),
                       units = sample(2:8, n_stat, replace = TRUE))

# ---- incidents --------------------------------------------------------------
n_inc <- 3000L
days  <- 14L
n_min <- days * 24L * 60L

# hotspot weights: sum of Gaussians centred on 6 random nodes
centres <- sample.int(g$n, 6L)
sigma <- 0.025                                   # ~2.7 km in degrees
w <- rep(0.02, g$n)
for (c0 in centres) w <- w + exp(-((g$x - g$x[c0])^2 + (g$y - g$y[c0])^2) / (2 * sigma^2))

# time of each incident: rejection sampling from a rate profile
rate_at <- function(t_min) {
  hour <- (t_min %/% 60) %% 24
  base <- 0.5 + 0.5 * sin((hour - 14) / 24 * 2 * pi)        # peaks mid-afternoon
  day  <- t_min %/% (24 * 60)
  storm <- if (day >= 8 && day <= 9) 3 else 0               # storm surge, days 8-9
  base + storm + 0.2
}
times <- integer(n_inc); k <- 0L
while (k < n_inc) {
  t <- sample.int(n_min, 1L) - 1L
  if (runif(1) < rate_at(t) / 4.2) { k <- k + 1L; times[k] <- t }
}
times <- sort(times)                             # data generation only (not a graded sort)

nodes <- sample.int(g$n, n_inc, replace = TRUE, prob = w)
types <- sample(c("fire", "flood", "accident", "medical", "collapse"), n_inc, replace = TRUE,
                prob = c(.15, .15, .30, .30, .10))
# floods are more likely during the storm
storm_idx <- (times %/% (24 * 60)) %in% 8:9
types[storm_idx & runif(n_inc) < 0.5] <- "flood"
severity <- sample(1:5, n_inc, replace = TRUE, prob = c(.3, .3, .2, .15, .05))
demand   <- severity * sample(1:3, n_inc, replace = TRUE)             # resource units
incidents <- data.frame(id = seq_len(n_inc), time_min = times, node = nodes,
                        lon = g$x[nodes], lat = g$y[nodes],
                        type = types, severity = severity, demand = demand)
# coarse zone = 3x3 grid cell of the bounding box (used by matrix-chain / evacuation)
zx <- pmin(3L, 1L + floor(3 * (incidents$lon - min(g$x)) / (max(g$x) - min(g$x) + 1e-9)))
zy <- pmin(3L, 1L + floor(3 * (incidents$lat - min(g$y)) / (max(g$y) - min(g$y) + 1e-9)))
incidents$zone <- (zy - 1L) * 3L + zx

# ---- call logs --------------------------------------------------------------
areas <- c("Anna Nagar", "T Nagar", "Adyar", "Velachery", "Tambaram", "Guindy", "Egmore",
           "Mylapore", "Besant Nagar", "Perambur", "Kodambakkam", "Royapuram", "Porur",
           "Thiruvanmiyur", "Ambattur", "Saidapet", "Chromepet", "Nungambakkam")
roads_ <- c("Mount Road", "OMR", "ECR", "GST Road", "Poonamallee High Road", "Inner Ring Road",
            "Cathedral Road", "Anna Salai", "100 Feet Road", "Kamarajar Salai")
phr <- list(
  fire     = c("fire reported at a shop", "smoke visible from a building", "gas cylinder blast"),
  flood    = c("water logging knee deep", "road submerged after heavy rain", "houses flooded"),
  accident = c("two vehicles collided", "bus overturned", "bike skidded and rider injured"),
  medical  = c("person collapsed and unconscious", "chest pain emergency", "pregnant woman needs ambulance"),
  collapse = c("wall collapsed", "old building partially collapsed", "tree fell on houses")
)
mk_text <- function(i) {
  ty <- incidents$type[i]
  paste0("INC", sprintf("%04d", i), " ", phr[[ty]][sample.int(3, 1)], " near ",
         areas[sample.int(length(areas), 1)], " on ", roads_[sample.int(length(roads_), 1)],
         " severity ", incidents$severity[i], " units needed ", incidents$demand[i])
}
noisy <- function(s) {                           # a second caller re-reports the same incident
  w <- strsplit(s, " ", fixed = TRUE)[[1]][-1]   # callers do not know the incident id: drop token 1
  for (j in seq_len(max(1L, length(w) %/% 8L))) { # one or two forgotten words, some SHOUTING
    p <- sample.int(length(w), 1L)
    if (runif(1) < 0.6) w[p] <- "" else w[p] <- toupper(w[p])
  }
  paste(w[nzchar(w)], collapse = " ")
}
txt <- character(n_inc); dup <- logical(n_inc); dup_of <- integer(n_inc)
for (i in seq_len(n_inc)) txt[i] <- mk_text(i)
n_dup <- round(0.15 * n_inc)
dsrc <- sample.int(n_inc, n_dup)
extra_txt <- character(n_dup)
for (j in seq_len(n_dup)) extra_txt[j] <- noisy(txt[dsrc[j]])
call_logs <- data.frame(
  id = seq_len(n_inc + n_dup),
  incident_id = c(seq_len(n_inc), dsrc),
  text = c(txt, extra_txt),
  is_duplicate = c(rep(FALSE, n_inc), rep(TRUE, n_dup)),
  dup_of = c(rep(NA_integer_, n_inc), dsrc),
  stringsAsFactors = FALSE
)

# ---- save -------------------------------------------------------------------
saveRDS(stations,  file.path(out_dir, "stations.rds"))
saveRDS(hospitals, file.path(out_dir, "hospitals.rds"))
saveRDS(incidents, file.path(out_dir, "incidents.rds"))
saveRDS(call_logs, file.path(out_dir, "call_logs.rds"))
writeLines(paste(call_logs$text, collapse = "\n"), file.path(out_dir, "call_logs.txt"))
cat(sprintf("incidents=%d  calls=%d  stations=%d  hospitals=%d\n",
            nrow(incidents), nrow(call_logs), nrow(stations), nrow(hospitals)))
