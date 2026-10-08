# Run the whole suite from the project root:   Rscript tests/testthat.R
library(testthat)
test_dir("tests/testthat", reporter = "progress", stop_on_failure = TRUE)
