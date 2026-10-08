# =============================================================================
# CityShield dashboard  -  run with:  shiny::runApp("app")   (from cityshield/)
#
# Design: every expensive result is PRECOMPUTED by scripts/05_precompute.R and
# only read here; cheap algorithms (matching, routes, knapsacks, small TSP, Karger,
# KMP search ...) run LIVE when you move a control. Each tab has a map or chart,
# a few summary cards and a short plain-English result. Where a problem has both an
# exact and an approximate solver, the app uses the exact one whenever it can finish.
# Code layout: one file per tab (tab_*.R) defining <tab>_ui() and <tab>_server().
# =============================================================================
library(shiny); library(bslib); library(leaflet); library(plotly)

# ---- locate the project and load the algorithm library ----
root <- Sys.getenv("CITYSHIELD_ROOT", unset = "")
if (!nzchar(root)) { d <- getwd(); while (!file.exists(file.path(d, "R", "core", "load_all.R")) && dirname(d) != d) d <- dirname(d); root <- d }
Sys.setenv(CITYSHIELD_ROOT = root)
setwd(root)
source(file.path(root, "R", "core", "load_all.R"))
for (f in c("helpers.R", "tab_dispatch.R", "tab_evacuation.R", "tab_resilience.R", "tab_placement.R",
            "tab_logistics.R", "tab_intel.R"))
  source(file.path(root, "app", f), local = TRUE)   # local: tab code sees G, SYN, PC defined below

# ---- data: the road graph, synthetic city data, precomputed results ---------------------
G   <- load_chennai_graph()
SYN <- load_synthetic()
pc_path <- function(n) data_path("precomputed", paste0(n, ".rds"))
need <- c("dispatch", "evacuation", "resilience", "placement", "intel")
missing <- need[!file.exists(vapply(need, pc_path, ""))]
if (length(missing)) stop("Missing precomputed files: ", paste(missing, collapse = ", "),
                          "\nRun  Rscript scripts/05_precompute.R  first.")
PC <- list(dispatch = readRDS(pc_path("dispatch")), evac = readRDS(pc_path("evacuation")),
           resil = readRDS(pc_path("resilience")), place = readRDS(pc_path("placement")),
           intel = readRDS(pc_path("intel")))
CENTER <- c(lon = mean(PC$dispatch$incidents$lon), lat = mean(PC$dispatch$incidents$lat))

theme <- bs_theme(version = 5, bg = "#fcfcfb", fg = "#0b0b0b", primary = "#2a78d6", secondary = "#52514e",
                  "navbar-bg" = "#0b0b0b", base_font = "system-ui")

ui <- page_navbar(
  title = "CityShield", theme = theme, fillable = FALSE, inverse = TRUE,
  header = tags$head(tags$style(HTML(".cs-kpi .value-box-value{font-size:1.4rem}.cs-kpi{margin-bottom:.5rem}
    .leaflet-container{border-radius:6px} .card-header{font-weight:600} table{font-size:.85rem}"))),
  dispatch_ui(), evacuation_ui(), resilience_ui(), placement_ui(),
  logistics_ui(), intel_ui()
)

server <- function(input, output, session) {
  dispatch_server(input, output, session)
  evacuation_server(input, output, session)
  resilience_server(input, output, session)
  placement_server(input, output, session)
  logistics_server(input, output, session)
  intel_server(input, output, session)
}

shinyApp(ui, server)
