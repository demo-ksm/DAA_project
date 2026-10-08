# =============================================================================
# Tab 3 - RESILIENCE
#   Karger's randomised min cut on a 120-node road study area (the weakest link of
#   the network; fixed 100 trials and a fixed seed) and the worst continuous
#   incident-load period (maximum subarray, divide and conquer).
# =============================================================================

KARGER_TRIALS <- 100L
KARGER_SEED <- 1L

resilience_ui <- function() {
  nav_panel("Resilience",
    layout_sidebar(
      sidebar = sidebar(width = 310,
        viva_box("weakest link: fewest road links whose loss disconnects the area",
                 "Karger's min cut: pick a random edge, contract its endpoints, repeat until 2 super-nodes remain; keep the best of many trials",
                 "O(m a(n)) per trial; about (n^2/2) ln(1/d) trials for failure chance d", "O(n + m)"),
        viva_box("busiest stretch of hours: the run of consecutive hours with the largest total load above capacity",
                 "Divide and conquer maximum subarray: the best stretch is in the left half, the right half, or crosses the middle",
                 "O(n log n)", "O(log n)")),
      layout_columns(col_widths = c(4, 8), uiOutput("rs_k1"), card(height = "100px", class = "justify-content-center", uiOutput("rs_text"))),
      card(card_header("Weakest link of the road network (red = road to cut)"), leafletOutput("rs_map", height = 430)),
      card(card_header("Worst stretch of incident load"), plotlyOutput("rs_load", height = 260), uiOutput("rs_load_text"))))
}

resilience_server <- function(input, output, session) {
  R <- PC$resil; und <- R$study$edges; sg <- R$study$graph

  weakest <- reactive(with_seed(KARGER_SEED, karger_min_cut(und$n, und$eu, und$ev, KARGER_TRIALS))$best)

  output$rs_k1 <- renderUI(kpi("Weakest link found", count_noun(weakest()$cut, "road"), "primary"))
  output$rs_text <- renderUI({
    k <- weakest()$cut
    tags$p(class = "mb-0", if (k == 1) "Removing this road would cut the area off from the rest of the network."
                           else sprintf("Removing these %d roads together would cut the area off from the rest of the network.", k))
  })

  output$rs_map <- renderLeaflet({
    b <- weakest()
    xy <- pair_xy(sg, und$eu, und$ev)
    map <- base_map(mean(sg$x), mean(sg$y), 14) |>
      addPolylines(lng = xy$lng, lat = xy$lat, color = PAL$grey, weight = 1.5, opacity = .7) |>
      addCircleMarkers(lng = sg$x, lat = sg$y, radius = 3, color = ifelse(b$side, PAL$blue, PAL$orange), weight = 1, fillOpacity = .9,
                       label = ifelse(b$side, "One side of the cut", "Other side of the cut"))
    if (length(b$cut_edges)) {
      cx <- pair_xy(sg, und$eu[b$cut_edges], und$ev[b$cut_edges])
      map <- addPolylines(map, lng = cx$lng, lat = cx$lat, color = PAL$red, weight = 6, opacity = 1, label = "Weakest link")
    }
    fitBounds(map, min(sg$x), min(sg$y), max(sg$x), max(sg$y))
  })

  output$rs_load <- renderPlotly({
    w <- R$window; h <- seq_along(R$load)
    plot_ly() |>
      add_bars(x = h, y = R$load, name = "Incidents per hour", marker = list(color = PAL$blue)) |>
      add_lines(x = c(1, length(h)), y = c(R$capacity, R$capacity), name = "Capacity", line = list(color = PAL$grey, dash = "dash")) |>
      plot_theme(xaxis = list(title = "Hour of the 14-day log"), yaxis = list(title = "Incidents"), legend = list(orientation = "h"),
                 shapes = list(list(type = "rect", x0 = w$start - .5, x1 = w$end + .5, y0 = 0, y1 = max(R$load), fillcolor = PAL$orange,
                                    opacity = .25, line = list(width = 0))))
  })
  output$rs_load_text <- renderUI({
    w <- R$window
    tags$div(class = "small", sprintf("Worst stretch: day %d at %02d:00 to day %d at %02d:00 (%d hours), with %s more incidents than the capacity line allows (shaded).",
      (w$start - 1) %/% 24 + 1, (w$start - 1) %% 24, (w$end - 1) %/% 24 + 1, w$end %% 24, w$end - w$start + 1L, format(round(w$sum), big.mark = ",")))
  })
}
