# =============================================================================
# Tab 6 - INTEL (incident logs)
#   Search reports       : KMP string matching on the call log
#   Find duplicate reports: LCS similarity of call-log lines (precomputed best matches)
# =============================================================================

intel_ui <- function() {
  nav_panel("Intel",
    navset_card_pill(
      nav_panel("Search reports",
        layout_sidebar(sidebar = sidebar(width = 310,
          textInput("in_kw", "Keyword or phrase", "flood"),
          selectInput("in_size", "Amount of log to search", c("Short (1,500 characters)" = 1500, "Medium (10,000 characters)" = 10000,
                                                              "Long (40,000 characters)" = 40000), selected = 10000),
          viva_box("find every place where a keyword or phrase appears in the call log",
                   "Knuth-Morris-Pratt (KMP): a table of the pattern's own prefixes tells how far to slide after a mismatch, so no text character is read twice",
                   "O(n + m) to search, O(m) to build the table", "O(m)")),
          layout_columns(col_widths = c(6, 6), uiOutput("in_k1"), uiOutput("in_k2")),
          card(card_header("Matches in context"), uiOutput("in_ctx")))),
      nav_panel("Find duplicate reports",
        layout_sidebar(sidebar = sidebar(width = 310,
          sliderInput("in_thr", "How similar two reports must be to count as duplicates", .5, 1, .9, step = .01),
          viva_box("two callers describe the same incident with different wording",
                   "Longest common subsequence (LCS) of the words, as a share of the longer report; similar enough means duplicate",
                   "O(mn) per pair of reports", "O(min(m, n)) with two rows")),
          layout_columns(col_widths = c(4, 4, 4), uiOutput("in_d1"), uiOutput("in_d2"), uiOutput("in_d3")),
          layout_columns(col_widths = c(6, 6), card(plotlyOutput("in_pr", height = 300)), card(tableOutput("in_dup_tbl"))))))
  )
}

intel_server <- function(input, output, session) {
  I <- PC$intel; calls <- I$calls
  full_log <- paste(calls$text, collapse = "\n")

  # ---- search (KMP)
  search <- reactive({
    kw <- input$in_kw; req(nzchar(kw)); txt <- substr(full_log, 1, as.integer(input$in_size))
    list(pos = kmp_match(str_codes(txt), str_codes(kw))$pos, txt = txt, kw = kw)
  })
  output$in_k1 <- renderUI(kpi("Matches found", length(search()$pos), "primary"))
  output$in_k2 <- renderUI(kpi("Characters searched", format(nchar(search()$txt), big.mark = ",")))
  output$in_ctx <- renderUI({
    s <- search()
    if (!length(s$pos)) return(tags$p("No matches."))
    show <- head(s$pos, 8)
    tagList(lapply(show, function(p) tags$div(class = "small font-monospace border-bottom py-1",
      substr(s$txt, max(1, p - 35), p - 1), tags$mark(substr(s$txt, p, p + nchar(s$kw) - 1)), substr(s$txt, p + nchar(s$kw), p + nchar(s$kw) + 35))),
      if (length(s$pos) > 8) tags$p(class = "small text-muted", sprintf("... and %d more", length(s$pos) - 8)))
  })

  # ---- duplicates: how many flags are right, and how many real duplicates are found, at each threshold
  b <- I$best; truth <- calls$is_duplicate
  pr_at <- function(th) { pr <- b$sim >= th; tp <- sum(pr & truth); c(correct = tp / max(1, sum(pr)), found = tp / sum(truth), flagged = sum(pr)) }
  output$in_d1 <- renderUI(kpi("Reports flagged as duplicates", pr_at(input$in_thr)[["flagged"]], "primary"))
  output$in_d2 <- renderUI(kpi("Flags that are correct", sprintf("%.0f%%", 100 * pr_at(input$in_thr)[["correct"]])))
  output$in_d3 <- renderUI(kpi("Real duplicates found", sprintf("%.0f%%", 100 * pr_at(input$in_thr)[["found"]])))
  output$in_pr <- renderPlotly({
    th <- seq(.5, 1, by = .01); m <- t(vapply(th, pr_at, c(0, 0, 0)))
    plot_ly(x = th) |> add_lines(y = m[, 1], name = "Flags that are correct", line = list(color = PAL$blue)) |>
      add_lines(y = m[, 2], name = "Real duplicates found", line = list(color = PAL$orange)) |>
      add_markers(x = input$in_thr, y = pr_at(input$in_thr)[["correct"]], showlegend = FALSE, marker = list(color = PAL$blue, size = 9)) |>
      add_markers(x = input$in_thr, y = pr_at(input$in_thr)[["found"]], showlegend = FALSE, marker = list(color = PAL$orange, size = 9)) |>
      plot_theme(xaxis = list(title = "Similarity needed"), yaxis = list(title = "", range = c(0, 1.02)), legend = list(orientation = "h"))
  })
  output$in_dup_tbl <- renderTable({
    flag <- which(b$sim >= input$in_thr & b$best_i > 0)
    flag <- flag[order(-b$sim[flag])]; flag <- head(flag, 6)
    if (!length(flag)) return(data.frame(Note = "Nothing flagged"))
    data.frame(Report = substr(calls$text[flag], 1, 70), `Similar earlier report` = substr(calls$text[b$best_i[flag]], 1, 70), Similarity = round(b$sim[flag], 2),
               `Really a duplicate` = ifelse(truth[flag], "yes", "no"), check.names = FALSE)
  })
}
