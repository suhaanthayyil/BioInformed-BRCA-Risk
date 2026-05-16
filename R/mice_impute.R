#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("mice is not installed; cannot run multiple imputation")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript R/mice_impute.R <input_csv> <output_csv>")
}

input_path <- args[[1]]
output_path <- args[[2]]

df <- read.csv(input_path, check.names = FALSE)
set.seed(42)
imp <- mice::mice(df, m = 5, seed = 42, printFlag = FALSE)
completed <- mice::complete(imp, action = "long", include = TRUE)
write.csv(completed, output_path, row.names = FALSE)
