# =============================================================================
# Tab 2 - EVACUATION
#   Push-relabel max flow on the road capacity network, sweep line (Bentley-
#   Ottmann) for the roads crossing the flood boundary, Graham scan for the
#   hull of the affected area. All precomputed.
# =============================================================================

evacuation_ui <- function() {
  nav_panel("Evacuation",
    layout_sidebar(
      sidebar = sidebar(width = 310,
        selectInput("ev_scen", "Flood scenario", choices = NULL),
        checkboxGroupInput("ev_layers", "Map layers", c("Evacuation flow" = "flow", "Roads cut by the flood" = "cross", "Affected area" = "hull"),
                           selected = c("flow", "cross", "hull")),
        viva_box("move people from the flooded zone to shelters through roads with limited capacity",
                 "Push-relabel max flow: super-source -> danger nodes, road edges = capacity (vehicles per hour x 3), shelters -> super-sink; max flow = min cut",
                 "O(V^3) for FIFO push-relabel", "O(V + E)"),
        viva_box("find the roads that cross the flood boundary",
                 "Sweep line (Bentley-Ottmann): sweep left to right, keep the roads cut by the line in a balanced tree and test only neighbours",
                 "O((n + k) log n) for k crossings", "O(n + k)"),
        viva_box("smallest convex area that contains every affected incident",
                 "Graham scan: sort points by angle around the lowest point, keep only left turns",
                 "O(n log n)", "O(n)")),
      layout_columns(col_widths = c(3, 3, 3, 3), uiOutput("ev_k1"), uiOutput("ev_k2"), uiOutput("ev_k3"), uiOutput("ev_k4")),
      layout_columns(col_widths = c(8, 4),
        card(card_header("Flood zone, evacuation flows and exit roads"), leafletOutput("ev_map", height = 540)),
        card(card_header("Summary"), uiOutput("ev_text")))))
}

evacuation_server <- function(input, output, session) {
  E <- PC$evac
  updateSelectInput(session, "ev_scen", choices = setNames(seq_along(E), paste("Scenario", seq_along(E))))
  sc <- reactive({ req(input$ev_scen); E[[as.integer(input$ev_scen)]] })
  from_cache <- new.env()

  output$ev_k1 <- renderUI({ s <- sc(); kpi("People to evacuate", format(sum(s$people), big.mark = ","), "primary") })
  output$ev_k2 <- renderUI({
    s <- sc(); kpi("Evacuated per hour", sprintf("%s (%.0f%%)", format(round(s$flow_total), big.mark = ","), 100 * s$flow_total / sum(s$people)))
  })
  output$ev_k3 <- renderUI(kpi("Exit roads cut by the flood", format(length(sc()$cross_edges), big.mark = ",")))
  output$ev_k4 <- renderUI(kpi("Affected area (km2)", sprintf("%.1f", sc()$hull_info$area_km2)))

  output$ev_map <- renderLeaflet({
    s <- sc(); g <- s$graph; layers <- input$ev_layers
    key <- paste0("from", s$id); if (is.null(from_cache[[key]])) from_cache[[key]] <- csr_from(g)
    from <- from_cache[[key]]
    map <- base_map(s$lon0, s$lat0, 14) |>
      addCircles(lng = s$lon0, lat = s$lat0, radius = s$r_flood, color = PAL$blue, weight = 2, fillOpacity = .07, label = "Edge of the flood zone")
    if ("hull" %in% layers) map <- addPolygons(map, lng = s$hull$lon, lat = s$hull$lat, color = PAL$magenta, weight = 2, fillOpacity = .12,
                                               label = sprintf("Affected area: %d incidents, %.1f km2", length(s$incidents_inside), s$hull_info$area_km2))
    if ("flow" %in% layers && !is.null(s$flows)) {
      fl <- s$flows; used <- which(fl > 1e-9)
      if (length(used)) {
        brk <- unique(stats::quantile(fl[used], c(0, .5, .85, 1)))
        for (b in seq_len(length(brk) - 1L)) {
          sel <- used[fl[used] >= brk[b] & (fl[used] < brk[b + 1L] | b == length(brk) - 1L)]
          if (!length(sel)) next
          xy <- edge_xy(g, sel, from)
          map <- addPolylines(map, lng = xy$lng, lat = xy$lat, color = PAL$aqua, weight = 1 + 2 * b, opacity = .75,
                              label = sprintf("%.0f to %.0f people per hour", brk[b], brk[b + 1L]))
        }
      }
    }
    if ("cross" %in% layers && length(s$cross_edges)) {
      xy <- edge_xy(g, s$cross_edges, from)
      map <- addPolylines(map, lng = xy$lng, lat = xy$lat, color = PAL$red, weight = 5, opacity = .95, label = "Road cut by the flood")
    }
    map |>
      addCircleMarkers(lng = g$x[s$danger], lat = g$y[s$danger], radius = 3, color = PAL$orange, weight = 1, fillOpacity = .9, label = "25 people at risk") |>
      addCircleMarkers(lng = g$x[s$shelters], lat = g$y[s$shelters], radius = 7, color = "#006300", fillColor = "#0ca30c", fillOpacity = 1, label = "Shelter")
  })

  output$ev_text <- renderUI({
    s <- sc(); people <- sum(s$people); cut_txt <- sprintf("The flood cuts %s.", count_noun(length(s$cross_edges), "road"))
    if (s$flow_total >= people - 1e-6) {
      tags$p(sprintf("Everyone can leave the flood zone within one hour. %s", cut_txt))
    } else {
      tags$p(sprintf("About %s of the %s people (%.0f%%) can leave the flood zone within one hour; road capacity and shelter space are the limit. %s",
                     format(round(s$flow_total), big.mark = ","), format(people, big.mark = ","), 100 * s$flow_total / people, cut_txt))
    }
  })
}
