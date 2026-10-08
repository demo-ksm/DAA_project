# =============================================================================
# app/helpers.R - small helpers shared by all dashboard tabs
# =============================================================================

# Validated categorical palette (fixed order) + a few semantic colours.
PAL <- list(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
            magenta = "#e87ba4", violet = "#4a3aa7", red = "#d03b3b", grey = "#898781", ink = "#0b0b0b")

#' Leaflet base map (light tiles) centred on the study region.
base_map <- function(lon, lat, zoom = 11) {
  leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
    addProviderTiles(providers$Esri.WorldGrayCanvas) |>      # free grey tiles, no API key (CartoDB tiles now need one)
    setView(lon, lat, zoom)
}

#' Source node of every CSR edge (the CSR stores only targets). O(E).
csr_from <- function(g) {
  from <- integer(g$m)
  for (u in seq_len(g$n)) for (k in seq_len(g$offset[u + 1L] - g$offset[u])) from[g$offset[u] + k - 1L] <- u
  from
}

#' Coordinates of a set of CSR edges as ONE polyline with NA breaks (fast to draw).
edge_xy <- function(g, edges, from = csr_from(g)) {
  m <- length(edges)
  lng <- numeric(3L * m); lat <- numeric(3L * m)
  for (i in seq_len(m)) {
    e <- edges[i]; a <- from[e]; b <- g$target[e]; k <- 3L * (i - 1L)
    lng[k + 1L] <- g$x[a]; lng[k + 2L] <- g$x[b]; lng[k + 3L] <- NA
    lat[k + 1L] <- g$y[a]; lat[k + 2L] <- g$y[b]; lat[k + 3L] <- NA
  }
  list(lng = lng, lat = lat)
}

#' Same for an undirected edge list (eu, ev) on a graph with node coordinates.
pair_xy <- function(g, eu, ev) {
  m <- length(eu); lng <- numeric(3L * m); lat <- numeric(3L * m)
  for (i in seq_len(m)) {
    k <- 3L * (i - 1L)
    lng[k + 1L] <- g$x[eu[i]]; lng[k + 2L] <- g$x[ev[i]]; lng[k + 3L] <- NA
    lat[k + 1L] <- g$y[eu[i]]; lat[k + 2L] <- g$y[ev[i]]; lat[k + 3L] <- NA
  }
  list(lng = lng, lat = lat)
}

#' Run `expr` with a fixed random seed, then restore the caller's random-number state.
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  expr
}

#' "1 road" / "72 roads" and similar plain counts.
count_noun <- function(n, singular, plural = paste0(singular, "s")) paste(format(n, big.mark = ","), if (n == 1) singular else plural)

#' Seconds as plain minutes: "about 13 minutes".
about_minutes <- function(sec) {
  m <- round(sec / 60)
  if (m < 1) "less than a minute" else paste("about", count_noun(m, "minute"))
}

#' A compact KPI tile.
kpi <- function(title, value, theme = "light") {
  bslib::value_box(title = title, value = value, theme = theme, height = "100px", class = "cs-kpi")
}

#' Small "problem / approach / time / space" note for a sidebar (the only place algorithm names appear).
viva_box <- function(problem, approach, time, space) {
  tags$div(class = "small border-start ps-2 my-2",
           tags$div(tags$b("Problem: "), problem), tags$div(tags$b("Approach: "), approach),
           tags$div(tags$b("Time: "), time, "  ", tags$b("Space: "), space))
}

#' Plotly layout defaults (transparent paper, palette fonts).
plot_theme <- function(p, ...) {
  plotly::layout(p, paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
                 font = list(color = "#52514e", size = 12), margin = list(l = 50, r = 10, t = 30, b = 40), ...)
}
