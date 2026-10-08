# Headless smoke test of the dashboard: sets the controls of every tab and renders
# every output through shiny::testServer(), reporting any error.
#
#   Rscript scripts/06_app_smoke_test.R        (CITYSHIELD_SYNTH=1 for the dev grid)

suppressPackageStartupMessages({ library(shiny); library(testthat) })
app <- shiny::shinyAppDir("app")
fails <- character(0)
try_out <- function(output, ids) {
  for (id in ids) {
    r <- tryCatch({ v <- output[[id]]; if (is.null(v)) stop("NULL output"); "ok" }, error = function(e) conditionMessage(e))
    if (!identical(r, "ok")) { fails <<- c(fails, paste0(id, ": ", r)); cat("FAIL", id, "-", r, "\n") } else cat("ok  ", id, "\n")
  }
}

shiny::testServer(app, {
  # ---- Dispatch
  session$setInputs(dp_inc = "1", dp_limit = 25)
  try_out(output, c("dp_k1", "dp_k2", "dp_map", "dp_text", "dp_table"))
  session$setInputs(dp_inc = "4", dp_limit = 10); try_out(output, c("dp_k1", "dp_k2", "dp_map", "dp_text"))
  # ---- Evacuation
  session$setInputs(ev_scen = "1", ev_layers = c("flow", "cross", "hull"))
  try_out(output, c("ev_k1", "ev_k2", "ev_k3", "ev_k4", "ev_map", "ev_text"))
  session$setInputs(ev_scen = "2", ev_layers = "cross"); try_out(output, c("ev_k2", "ev_map", "ev_text"))
  # ---- Resilience (no inputs: trials and seed are constants)
  try_out(output, c("rs_k1", "rs_text", "rs_map", "rs_load", "rs_load_text"))
  # ---- Placement (radius 1 uses the approximation, radius 2 and 3 the exact search)
  for (r in c("1", "2", "3")) { session$setInputs(pl_sc_rad = r); try_out(output, c("pl_k1", "pl_k2", "pl_sc_map", "pl_sc_text")) }
  try_out(output, c("pl_vk1", "pl_vk2", "pl_vc_map", "pl_vc_text"))
  # ---- Logistics
  session$setInputs(lg_cap = 40, lg_budget = 45, lg_k = 7)
  try_out(output, c("lg_k1", "lg_k2", "lg_sup_plot", "lg_u1", "lg_u2", "lg_up_tbl", "lg_t1", "lg_t2", "lg_tsp_map", "lg_tsp_text"))
  for (k in c(4, 10, 11, 14)) { session$setInputs(lg_k = k); try_out(output, c("lg_t1", "lg_tsp_map", "lg_tsp_text")) }
  # ---- Intel
  session$setInputs(in_kw = "flood", in_size = "10000", in_thr = .9)
  try_out(output, c("in_k1", "in_k2", "in_ctx", "in_d1", "in_d2", "in_d3", "in_pr", "in_dup_tbl"))
  session$setInputs(in_size = "1500", in_kw = "no such phrase", in_thr = 1); try_out(output, c("in_k1", "in_ctx", "in_dup_tbl"))
})
cat(sprintf("\n%d failures\n", length(fails)))
if (length(fails)) { cat(paste(fails, collapse = "\n"), "\n"); quit(status = 1) }
