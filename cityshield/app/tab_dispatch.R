# =============================================================================
# Tab 1 - DISPATCH
#   Routes (Dijkstra with early exit), ambulance-to-incident assignment (bipartite
#   matching via push-relabel max flow) and incident ranking (quicksort, done
#   silently to pick the ten most urgent incidents).
# =============================================================================

dispatch_ui <- function() {
  nav_panel("Dispatch",
    layout_sidebar(
      sidebar = sidebar(width = 310,
        selectInput("dp_inc", "Choose an incident", choices = NULL),
        sliderInput("dp_limit", "Longest ambulance trip allowed (minutes)", 5, 60, 25, step = 5),
        viva_box("give each incident at most one ambulance that can reach it in time, serving as many incidents as possible",
                 "Push-relabel max flow: ambulances and incidents form a bipartite graph, turned into a unit-capacity flow network; a flow of k is a matching of size k",
                 "O(V^3) for FIFO push-relabel", "O(V + E)"),
        viva_box("fastest road route from the nearest station to the chosen incident",
                 "Dijkstra's algorithm with a binary heap, stopping as soon as the incident is reached",
                 "O((V + E) log V)", "O(V + E)"),
        viva_box("show the most urgent incidents first",
                 "Quicksort: partition around a pivot, most severe first, then largest demand",
                 "expected O(n log n)", "O(log n)")),
      layout_columns(col_widths = c(6, 6), uiOutput("dp_k1"), uiOutput("dp_k2")),
      layout_columns(col_widths = c(8, 4),
        card(card_header("Incidents, stations and the selected route"), leafletOutput("dp_map", height = 520)),
        card(card_header("Result"), uiOutput("dp_text"), tableOutput("dp_table")))))
}

dispatch_server <- function(input, output, session) {
  D <- PC$dispatch; inc <- D$incidents; st <- SYN$stations; hos <- SYN$hospitals
  batch <- D$batch; amb <- D$amb

  # the ten most urgent incidents of the batch (quicksort, most severe first, then largest demand)
  top <- batch[rank_incidents(inc$severity[batch], inc$demand[batch], randomized = FALSE)][seq_len(min(10L, length(batch)))]
  cap1 <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  labs <- sprintf("%s, severity %d (ID %d)", cap1(inc$type[top]), inc$severity[top], inc$id[top])
  updateSelectInput(session, "dp_inc", choices = setNames(seq_along(top), labs))

  # --- ambulance units (one left node per ambulance) vs the batch incidents
  matching <- reactive({
    lim <- input$dp_limit * 60
    unit_station <- rep(seq_len(nrow(amb)), amb$units)
    l <- integer(0); r <- integer(0)
    for (a in seq_along(unit_station)) for (j in seq_along(batch))
      if (D$times[unit_station[a], j] <= lim) { l <- c(l, a); r <- c(r, j) }
    list(res = max_bipartite_matching(length(unit_station), length(batch), l, r), unit_station = unit_station)
  })

  # --- live route: Dijkstra from the nearest station to the chosen incident (early exit at the target)
  route <- reactive({
    req(input$dp_inc); i <- top[as.integer(input$dp_inc)]
    s_node <- st$node[match(inc$station[i], st$id)]
    r <- dijkstra(G, s_node, target = inc$node[i])
    list(path = path_from_parent(r$parent, s_node, inc$node[i], r$dist), eta = r$dist[inc$node[i]], incident = i)
  })

  # average time the assigned ambulances need to reach their incidents
  avg_arrival <- reactive({
    m <- matching(); pr <- m$res$pairs
    if (!nrow(pr)) return(NA_real_)
    t <- 0; for (k in seq_len(nrow(pr))) t <- t + D$times[m$unit_station[pr[k, 1]], pr[k, 2]]
    t / nrow(pr)
  })

  output$dp_k1 <- renderUI(kpi("Incidents served", sprintf("%d of %d", matching()$res$size, length(batch)), "primary"))
  output$dp_k2 <- renderUI({ a <- avg_arrival(); kpi("Average ambulance arrival time", if (is.na(a)) "-" else sprintf("%.1f min", a / 60)) })

  output$dp_map <- renderLeaflet({
    m <- matching(); res <- m$res; r <- route()
    served <- res$match_right > 0L
    map <- base_map(CENTER[["lon"]], CENTER[["lat"]], 11) |>
      addCircleMarkers(lng = st$lon, lat = st$lat, radius = 5, color = PAL$blue, fillOpacity = .8, weight = 1,
                       label = paste0("Station ", st$id, " (", st$type, ", ", st$units, " units)"), group = "Stations") |>
      addCircleMarkers(lng = hos$lon, lat = hos$lat, radius = 5, color = PAL$aqua, fillOpacity = .9, weight = 1,
                       label = paste0("Hospital ", hos$id, " (", hos$beds, " beds)"), group = "Hospitals") |>
      addCircleMarkers(lng = inc$lon[batch], lat = inc$lat[batch], radius = 4 + inc$severity[batch],
                       color = ifelse(served, PAL$orange, PAL$red), fillOpacity = .85, weight = 1,
                       label = paste0("Incident ", inc$id[batch], ": ", inc$type[batch], ", severity ", inc$severity[batch],
                                      ifelse(served, " - ambulance assigned", " - no ambulance within reach")), group = "Incidents")
    # matched pairs
    if (res$size > 0) {
      pr <- res$pairs
      for (k in seq_len(nrow(pr))) {
        s <- amb[m$unit_station[pr[k, 1]], ]; ii <- batch[pr[k, 2]]
        map <- addPolylines(map, lng = c(s$lon, inc$lon[ii]), lat = c(s$lat, inc$lat[ii]), color = PAL$orange, weight = 1.5, opacity = .6)
      }
    }
    if (length(r$path) > 1) {
      map <- addPolylines(map, lng = G$x[r$path], lat = G$y[r$path], color = PAL$violet, weight = 4, opacity = .9,
                          label = sprintf("Fastest route: %.1f min", r$eta / 60))
      map <- setView(map, inc$lon[r$incident], inc$lat[r$incident], 13)
    }
    addLayersControl(map, overlayGroups = c("Stations", "Hospitals", "Incidents"), options = layersControlOptions(collapsed = TRUE))
  })

  output$dp_text <- renderUI({
    k <- matching()$res$size; n <- length(batch); r <- route()
    reach <- if (k == n) sprintf("All %d incidents can be reached within %d minutes.", n, input$dp_limit)
             else sprintf("%d of %d incidents can be reached within %d minutes.", k, n, input$dp_limit)
    tags$p(sprintf("%s Route to incident %d: %s.", reach, inc$id[r$incident], about_minutes(r$eta)))
  })

  output$dp_table <- renderTable({
    data.frame(`Incident ID` = inc$id[top], Type = cap1(inc$type[top]), Severity = inc$severity[top],
               `Minutes to reach` = round(inc$station_time[top] / 60, 1), check.names = FALSE)
  })
}
