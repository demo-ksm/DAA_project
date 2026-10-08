# =============================================================================
# Compliance test: graded algorithms must be written with explicit loops.
#
# Scans every R source file under R/ (comments and string literals are removed
# first, so documentation can still MENTION these names) and fails if any
# banned shortcut is CALLED.
#
# Exempt: core/data_io.R, core/timer.R, core/load_all.R  (loading/benchmark
#         plumbing, not graded algorithms).
# sort() is additionally allowed in files named *_ref.R: reference baselines.
# =============================================================================

banned <- c("outer", "which", "which.max", "which.min", "order", "rank",
            "apply", "sapply", "lapply", "vapply", "tapply", "mapply", "rapply",
            "Map", "Reduce", "Filter", "Position", "Find", "Vectorize",
            "sort", "sort.int", "sort.list")

strip_code <- function(lines) {
  lines <- gsub('"([^"\\\\]|\\\\.)*"', '""', lines)   # drop "string literals"
  lines <- gsub("'([^'\\\\]|\\\\.)*'", "''", lines)   # drop 'string literals'
  gsub("#.*$", "", lines)                              # drop comments
}

find_violations <- function(path, banned) {
  code <- strip_code(readLines(path, warn = FALSE))
  hits <- character(0)
  for (fn in banned) {
    # a CALL: name preceded by a non-identifier char and followed by "("
    pat <- paste0("(^|[^A-Za-z0-9_.$@])", gsub(".", "\\.", fn, fixed = TRUE), "\\s*\\(")
    for (i in seq_along(code)) {
      if (grepl(pat, code[i], perl = TRUE))
        hits <- c(hits, sprintf("%s:%d  %s()", basename(path), i, fn))
    }
  }
  hits
}

test_that("the scanner itself detects banned calls (self-test)", {
  tmp <- tempfile(fileext = ".R")
  writeLines(c("x <- sapply(1:3, f)",
               "# which(x) in a comment is fine",
               "msg <- 'order(x) in a string is fine'",
               "y <- sort(x)",
               "z <- my_order(x)",          # different identifier: fine
               "w <- obj$which(1)"),        # method on an object: fine
             tmp)
  v <- find_violations(tmp, banned)
  expect_length(v, 2)
  expect_true(any(grepl("sapply", v)))
  expect_true(any(grepl("sort", v)))
})

test_that("no graded algorithm uses a banned shortcut", {
  files <- list.files(cs_path("R"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
  exempt <- c("data_io.R", "timer.R", "load_all.R")
  files <- files[!(basename(files) %in% exempt)]
  all_hits <- character(0)
  for (f in files) {
    b <- if (grepl("_ref\\.R$", f)) setdiff(banned, c("sort", "sort.int", "sort.list")) else banned
    all_hits <- c(all_hits, find_violations(f, b))
  }
  expect_equal(all_hits, character(0), info = paste(all_hits, collapse = "\n"))
})
