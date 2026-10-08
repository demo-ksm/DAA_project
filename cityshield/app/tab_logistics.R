# =============================================================================
# Tab 5 - LOGISTICS (all LIVE: each algorithm is fast at these sizes)
#   Supplies         : fractional knapsack (greedy)  - divisible relief supplies
#   Upgrades         : 0/1 knapsack (dynamic programming)
#   Inspection route : TSP - exact branch and bound when it finishes within its
#                      node budget, otherwise the MST 2-approximation
# =============================================================================

UPGRADES <- data.frame(
  name = c("Raise Adyar bridge", "Pump station T Nagar", "Backup power: Egmore hospital", "Storm-drain Velachery",
           "Sensor mesh OMR", "Radio mast Tambaram", "Flood gates Cooum", "Road widening GST", "Shelter retrofit Porur",
           "Drone depot Guindy", "Culvert Perambur", "Ambulance bay Anna Nagar"),
  cost = c(18, 9, 6, 14, 5, 4, 12, 20, 8, 7, 6, 5), benefit = c(40, 22, 20, 31, 15, 11, 27, 35, 19, 16, 12, 13), stringsAsFactors = FALSE)

TSP_NODE_CAP <- 500L       # the exact search "finishes" if it needs at most this many expansions (keeps every click under a second)

logistics_ui <- function() {
  nav_panel("Logistics",
    navset_card_pill(
      nav_panel("Supplies",
        layout_sidebar(sidebar = sidebar(width = 300,
          sliderInput("lg_cap", "Supply available (% of total demand)", 10, 100, 40, step = 5),
          viva_box("split divisible supplies among zones to get the most urgency served",
                   "Fractional knapsack (greedy): rank zones by urgency per unit and fill in that order; the last zone gets a fraction",
                   "O(n log n)", "O(n)")),
          layout_columns(col_widths = c(6, 6), uiOutput("lg_k1"), uiOutput("lg_k2")),
          card(plotlyOutput("lg_sup_plot", height = 360)))),
      nav_panel("Upgrades",
        layout_sidebar(sidebar = sidebar(width = 300,
          sliderInput("lg_budget", "Budget (crore)", 10, 100, 45, step = 1),
          viva_box("choose all-or-nothing upgrades under a budget to get the most resilience benefit",
                   "0/1 knapsack (dynamic programming) over (upgrade, budget)",
                   "O(nW)", "O(nW)")),
          layout_columns(col_widths = c(6, 6), uiOutput("lg_u1"), uiOutput("lg_u2")),
          card(tableOutput("lg_up_tbl")))),
      nav_panel("Inspection route",
        layout_sidebar(sidebar = sidebar(width = 300,
          sliderInput("lg_k", "Inspection sites", 4, 14, 7),
          viva_box("shortest closed tour through all sites (travel times on the road network)",
                   "Travelling salesman: branch and bound with a reduced-cost-matrix bound when it can finish; otherwise the MST 2-approximation (walk the spanning tree in preorder)",
                   "exact O(n!) worst; approximation O(n^2 log n)", "O(n^2 * frontier)")),
          layout_columns(col_widths = c(6, 6), uiOutput("lg_t1"), uiOutput("lg_t2")),
          layout_columns(col_widths = c(8, 4), card(leafletOutput("lg_tsp_map", height = 430)), card(uiOutput("lg_tsp_text"))))))
  )
}

logistics_server <- function(input, output, session) {
  inc <- PC$dispatch$incidents; P <- PC$place; g <- P$graph

  # ---- supplies: zones from the incident log
  zones <- local({
    z <- data.frame(zone = 1:9, demand = 0, urgency = 0)
    for (i in seq_len(nrow(inc))) { k <- inc$zone[i]; z$demand[k] <- z$demand[k] + inc$demand[i]; z$urgency[k] <- z$urgency[k] + inc$severity[i] * inc$demand[i] }
    z[z$demand > 0, ]
  })
  sup <- reactive({
    W <- sum(zones$demand) * input$lg_cap / 100
    list(res = fractional_knapsack(zones$urgency, zones$demand, W), W = W)
  })
  output$lg_k1 <- renderUI(kpi("Urgency served", sprintf("%.0f of %.0f (%.0f%%)", sup()$res$total, sum(zones$urgency), 100 * sup()$res$total / sum(zones$urgency)), "primary"))
  output$lg_k2 <- renderUI(kpi("Supply used", sprintf("%.0f of %.0f units", sup()$res$used, sup()$W)))
  output$lg_sup_plot <- renderPlotly({
    f <- sup()$res$fraction
    plot_ly(x = paste("Zone", zones$zone), y = round(100 * f, 1), type = "bar", marker = list(color = PAL$blue),
            text = sprintf("%.0f%% of %.0f units", 100 * f, zones$demand), hoverinfo = "text") |>
      plot_theme(yaxis = list(title = "% of the zone's demand delivered", range = c(0, 105)), title = "Most urgent zones are served first")
  })

  # ---- upgrades
  up <- reactive(knapsack01(UPGRADES$cost, UPGRADES$benefit, input$lg_budget))
  output$lg_u1 <- renderUI(kpi("Total benefit", up()$value, "primary"))
  output$lg_u2 <- renderUI(kpi("Budget used", sprintf("%d of %d crore", sum(UPGRADES$cost[up()$items]), input$lg_budget)))
  output$lg_up_tbl <- renderTable({ sel <- up()$items; d <- UPGRADES; d$chosen <- ifelse(seq_len(nrow(d)) %in% sel, "yes", ""); names(d) <- c("Upgrade", "Cost", "Benefit", "Chosen"); d })

  # ---- inspection route: exact branch and bound if it finishes, otherwise the 2-approximation
  sites <- reactive({ req(input$lg_k); s <- with_seed(3, sample.int(g$n, 14)); s[seq_len(input$lg_k)] })
  tsp <- reactive({
    s <- sites(); D <- P$D[s, s]
    ap <- tsp_approx(D)
    bb <- tsp_bb(D, "lc", ub = ap$cost + 1e-6, max_nodes = TSP_NODE_CAP)   # the approximation gives the exact search a head start
    exact <- bb$complete && length(bb$tour) == length(s)
    list(res = if (exact) bb else ap, exact = exact, s = s)
  })
  output$lg_t1 <- renderUI(kpi("Tour length", sprintf("%.1f min", tsp()$res$cost / 60), "primary"))
  output$lg_t2 <- renderUI(kpi("Sites to visit", input$lg_k))
  output$lg_tsp_map <- renderLeaflet({
    x <- tsp(); tour <- x$res$tour; s <- x$s
    m <- base_map(mean(g$x[s]), mean(g$y[s]), 13)
    if (length(tour) > 1) {
      o <- c(tour, tour[1]); m <- addPolylines(m, lng = g$x[s[o]], lat = g$y[s[o]], color = PAL$violet, weight = 3, opacity = .9)
    }
    m <- addCircleMarkers(m, lng = g$x[s], lat = g$y[s], radius = 8, color = PAL$blue, fillColor = PAL$blue, fillOpacity = 1,
                          label = paste("Site", seq_along(s)))
    fitBounds(m, min(g$x[s]) - .002, min(g$y[s]) - .002, max(g$x[s]) + .002, max(g$y[s]) + .002)
  })
  output$lg_tsp_text <- renderUI({
    x <- tsp()
    tags$div(class = "small",
      tags$p("Distances are shortest road travel times between the sites."),
      tags$p(if (x$exact)
        "The exact method was used because this number of sites is small enough for it to finish. This is the shortest possible tour."
      else
        "The fast approximate method was used because this many sites is too many for the exact search to finish. The tour is never more than twice the shortest possible."))
  })
}
