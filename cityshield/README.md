# CityShield

Emergency-response and infrastructure-resilience planner on the **real Chennai road network
(OpenStreetMap)**: where to place resources, how to dispatch responders, how to evacuate,
and how fragile the network is. A Design and Analysis of Algorithms project: every algorithm the
dashboard uses is implemented **from scratch in R** (explicit loops, hand-written queue / stack /
heap / treap / union-find), tested against a reference, and demonstrated in a Shiny dashboard.
The code contains only those algorithms.

* Report: `report/cityshield_report.Rmd` (render with `rmarkdown::render`) - viva notes: `report/viva_notes.md`
* Dashboard: `shiny::runApp("app")` (from this folder)

## Reproduce everything

Run from this folder. `Rscript` = `C:\Program Files\R\R-4.6.1\bin\Rscript.exe` (R 4.6.1).

```
Rscript -e "renv::init(bare = TRUE, restart = FALSE)"   # once: creates renv/ + .Rprofile
Rscript scripts/01_setup_renv.R      # install packages -> renv.lock (cache/library outside OneDrive, see .Renviron)

# road network - ONE download, cached with saveRDS():
Rscript scripts/02_download_osm.R    # Overpass tiles (resumable)  ... or, if Overpass is overloaded:
curl -L -C - -o C:/cityshield-raw/southern-zone.osm.pbf https://download.geofabrik.de/asia/india/southern-zone-latest.osm.pbf
Rscript scripts/02b_extract_osm_geofabrik.R    # GDAL reads the .pbf -> data/osm/chennai_roads_raw.rds

Rscript scripts/03_build_graph.R     # ways -> contracted CSR graph (largest strongly connected component)
Rscript scripts/04_generate_data.R   # synthetic incidents, call logs, stations, hospitals (seeded)
Rscript tests/testthat.R             # test suite (reference / igraph / sf validation + banned-function scan)
Rscript scripts/05_precompute.R      # results the dashboard reads (data/precomputed/*.rds)
Rscript benchmarks/run_all.R         # the 5 comparisons (05-09) -> benchmarks/results/*.csv + *.png + notes.rds
Rscript benchmarks/algorithm_registry.R   # one benchmark entry per algorithm
Rscript scripts/07_test_summary.R    # per-file test counts for the report
Rscript scripts/06_app_smoke_test.R  # headless test of every dashboard output
```

Development mode: `CITYSHIELD_SYNTH=1` swaps the road network for a synthetic grid and uses `data_synth/`,
so nothing in `data/` is touched.

## Dashboard

Six tabs, in this order: **Dispatch, Evacuation, Resilience, Placement, Logistics, Intel.**

| Tab | Algorithms |
|---|---|
| Dispatch | push-relabel max flow (ambulance matching), Dijkstra with early exit (routes), quicksort (ranks incidents) |
| Evacuation | push-relabel max flow, sweep line / Bentley-Ottmann (roads crossing the flood boundary), Graham scan (affected-area hull) |
| Resilience | Karger min cut (100 trials, fixed seed), divide-and-conquer maximum subarray (worst load period) |
| Placement | set cover (stations) and vertex cover (sensors): exact branch and bound when it finishes, otherwise greedy / maximal-matching 2-approximation |
| Logistics | fractional knapsack (supplies), 0/1 knapsack DP (upgrades), TSP exact branch and bound or MST 2-approximation (inspection route) |
| Intel | KMP (search reports), LCS (find duplicate reports) |

## Layout

```
R/core/         graph.R (CSR + residual flow graph)  structures.R (queue, stack, binary heap, 2-key heap)
                sort_utils.R (merge sort)  geometry_types.R  timer.R (time_it)  study_area.R  data_io.R  load_all.R
R/m1_greedy_dc/ fractional knapsack, maximum subarray (divide and conquer)
R/m2_dp_bt_bb/  LCS, 0/1 knapsack, generic branch & bound + TSP B&B
R/m3_strings/   KMP
R/m4_graphs/    Dijkstra, push-relabel max flow, bipartite matching, evacuation / dispatch glue
R/m5_geometry/  segment tests, treap, sweep line (Bentley-Ottmann), Graham scan
R/m6_randomized/ quicksort, Karger (+ union-find)
R/m7_complexity/ vertex cover, set cover (exact + approximate), TSP 2-approximation
tests/testthat/ one file per module + helper graphs + banned-function compliance test
benchmarks/     common.R, exp_*.R (experiments 05-09), notes.R (text from CSVs), algorithm_registry.R, results/
app/            Shiny dashboard: app.R + one file per tab
scripts/        pipeline scripts 01-07
report/         R Markdown report, viva notes, screenshots
```

## Rules enforced

`tests/testthat/test-banned-functions.R` fails if any file under `R/` *calls* `outer`, `which`,
`order`, `rank`, `apply`-family, `Map`/`Reduce`/`Filter`, or `sort` (the latter allowed only in
`*_ref.R`). `igraph` and `sf` are used only to load data and to validate results in tests (and, for igraph, in the Karger benchmark).
