# =============================================================================
# Tab 4 - PLACEMENT
#   Stations : set cover    (exact branch and bound when it finishes, otherwise greedy)
#   Sensors  : vertex cover (exact branch and bound when it finishes, otherwise the
#                            maximal-matching 2-approximation)
# The exact search was run (with a node cap) by scripts/05_precompute.R; its
# `complete` flag says whether it finished, which is the rule used here.
# =============================================================================

placement_ui <- function() {
  nav_panel("Placement",
    navset_card_pill(
      nav_panel("Stations",
        layout_sidebar(sidebar = sidebar(width = 300,
          selectInput("pl_sc_rad", "How far a station can reach", choices = NULL),
          viva_box("fewest stations so that every place is within reach of one",
                   "Exact branch and bound (branch on the sets covering the rarest uncovered place) when it can finish; otherwise greedy set cover (always pick the station covering most uncovered places)",
                   "exact O(f^k); greedy O(S sum|s|)", "O(U + S)")),
          layout_columns(col_widths = c(6, 6), uiOutput("pl_k1"), uiOutput("pl_k2")),
          layout_columns(col_widths = c(8, 4), card(leafletOutput("pl_sc_map", height = 450)), card(uiOutput("pl_sc_text"))))),
      nav_panel("Sensors",
        layout_sidebar(sidebar = sidebar(width = 300,
          viva_box("every road link needs a sensor at one end; use as few sensors as possible",
                   "Exact branch and bound when it can finish; otherwise take both ends of every edge of a maximal matching (2-approximation)",
                   "exact O(2^k); approximation O(n + m)", "O(n)")),
          layout_columns(col_widths = c(6, 6), uiOutput("pl_vk1"), uiOutput("pl_vk2")),
          layout_columns(col_widths = c(8, 4), card(leafletOutput("pl_vc_map", height = 450)), card(uiOutput("pl_vc_text"))))))
  )
}

placement_server <- function(input, output, session) {
  P <- PC$place; g <- P$graph
  updateSelectInput(session, "pl_sc_rad", choices = setNames(seq_along(P$set_cover),
    sprintf("Within %.0f seconds", vapply(P$set_cover, function(x) x$radius_s, 0))))

  # ---- stations (set cover): exact when its search finished, otherwise greedy
  sc <- reactive({
    req(input$pl_sc_rad); x <- P$set_cover[[as.integer(input$pl_sc_rad)]]
    exact <- isTRUE(x$exact$complete) && !is.na(x$exact$size)
    chosen <- if (exact) x$exact$chosen else x$greedy$chosen
    sets <- coverage_sets(P$D, x$radius_s); cov <- logical(g$n); for (s in chosen) cov[sets[[s]]] <- TRUE
    list(x = x, exact = exact, chosen = chosen, covered = cov)
  })
  output$pl_k1 <- renderUI(kpi("Stations chosen", length(sc()$chosen), "primary"))
  output$pl_k2 <- renderUI({ n <- 0L; for (v in sc()$covered) if (v) n <- n + 1L; kpi("Places covered", sprintf("%d of %d", n, g$n)) })
  output$pl_sc_map <- renderLeaflet({
    s <- sc(); ch <- s$chosen; cov <- s$covered
    base_map(mean(g$x), mean(g$y), 14) |>
      addCircleMarkers(lng = g$x, lat = g$y, radius = 3, color = ifelse(cov, PAL$grey, PAL$red), weight = 1, fillOpacity = .8, label = ifelse(cov, "Covered", "Not covered")) |>
      addCircleMarkers(lng = g$x[ch], lat = g$y[ch], radius = 8, color = PAL$blue, fillColor = PAL$blue, fillOpacity = .9, label = "Station") |>
      fitBounds(min(g$x), min(g$y), max(g$x), max(g$y))
  })
  output$pl_sc_text <- renderUI({
    s <- sc(); n <- length(s$chosen)
    tags$div(class = "small", tags$p(if (s$exact)
      sprintf("The exact method was used because this area is small enough for it to finish. %d stations is the smallest number that covers every place.", n)
    else
      sprintf("The fast approximate method was used because this area is too large for the exact search to finish. It chose %d stations, which is never more than %.1f times the smallest possible number.", n, harmonic_bound(s$x$max_set))))
  })

  # ---- sensors (vertex cover): exact when its search finished, otherwise the 2-approximation
  vc <- P$vertex_cover; e <- P$edges
  vc_exact <- isTRUE(vc$exact$complete)
  vc_cover <- if (vc_exact) vc$exact$cover else vc$approx$cover
  output$pl_vk1 <- renderUI(kpi("Sensors", length(vc_cover), "primary"))
  on_vc <- logical(g$n); on_vc[vc_cover] <- TRUE
  watched <- 0L; for (k in seq_along(e$eu)) if (on_vc[e$eu[k]] || on_vc[e$ev[k]]) watched <- watched + 1L
  output$pl_vk2 <- renderUI(kpi("Road links watched", sprintf("%d of %d", watched, length(e$eu))))
  output$pl_vc_map <- renderLeaflet({
    xy <- pair_xy(g, e$eu, e$ev); on <- logical(g$n); on[vc_cover] <- TRUE
    base_map(mean(g$x), mean(g$y), 14) |>
      addPolylines(lng = xy$lng, lat = xy$lat, color = PAL$grey, weight = 1.5) |>
      addCircleMarkers(lng = g$x[!on], lat = g$y[!on], radius = 2, color = PAL$grey, weight = 1) |>
      addCircleMarkers(lng = g$x[on], lat = g$y[on], radius = 6, color = PAL$violet, fillColor = PAL$violet, fillOpacity = .9, label = "Sensor") |>
      fitBounds(min(g$x), min(g$y), max(g$x), max(g$y))
  })
  output$pl_vc_text <- renderUI(tags$div(class = "small", tags$p(if (vc_exact)
    sprintf("The exact method was used because this area is small enough for it to finish. %d sensors is the smallest number that watches all %d road links.", length(vc_cover), length(e$eu))
  else
    sprintf("The fast approximate method was used because this area is too large for the exact search to finish. It chose %d sensors, which is never more than twice the smallest possible number.", length(vc_cover)))))
}
